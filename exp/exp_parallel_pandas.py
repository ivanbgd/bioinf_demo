#!/usr/bin/env python
"""
File:       exp/exp_parallel_pandas.py
Author:     Ivan Lazarević
Brief:      Script for experimenting with parallelization of the main logic using third party libraries.

Details:    Trying out Dask and Modin in combination with Pandas.
            Haven't implemented Dask.
            Have implemented Modin. Works correctly.
            Speed is about the same as sequential.
            Input file is not large enough, so parallelization overhead negates the benefits.
            Modin is **very** easy to translate to from Pandas.
            Modin preserves the order of records in the output file.
            Modin can be used in a single computer (shared-memory model), or in a cluster (distributed-memory model).
            It can read and process larger-than-memory files.
"""
# Standard library imports
import csv
import gzip
import os
import sys
from itertools import zip_longest
from pathlib import Path
from typing import List, Set, Tuple, Union

os.environ["MODIN_CPUS"] = "4"

# Third party library imports
import dask
import dask.dataframe as ddf
import modin
import modin.pandas as pd
import numpy as np
import pandas
import pgzip
from distributed import Client

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import ALPHABET, NEWLINE, ADAPTER_LEN, POLY_LEN, POLY_END, ADAPTER_END, OPEN_PARAMS
from src.config import PANDAS_SEPARATOR, PANDAS_COLUMNS
from src.config import INPUT_FASTQ, INPUT_ADAPTER
from src.config import OUTPUT_FASTQ_GZ, TEST_OUT_FASTQ_SIZE_REF, OUTPUT_STATISTICS, OUTPUT_TEXT_FILE
from src.config import OUTPUT_FASTQ_SMALL_GZ
from src.config import TEST_OUT_FASTQ_GZ_REFERENCE, TEST_OUT_STAT_REFERENCE
from src.config import TEST_INP_FASTQ_SMALL_ORIGINAL
from src.utils import exit_program, time_it
from src.trie import build_trie, trie_matching
from src.type_aliases import Trie, AdaptersNaive
from tools.make_gzip import compress_file


def _generate_all_polyx_patterns() -> Set[str]:
    """ Generate all *Poly-X* patterns and return them in a set.

        All *Poly-X* patterns are `POLY_LEN` long.
        They include the four original *Poly-X* patterns, without a mutation,
        and also additional 180 patterns with exactly one mutation.
        This can then be used for exact pattern matching, instead of approximate pattern matching.
    """
    all_polys = []

    for x in range(len(ALPHABET)):
        for i in range(POLY_LEN):
            current_letter = ALPHABET[x]
            poly_x = current_letter * POLY_LEN
            for letter in ALPHABET:
                pattern = poly_x[:i] + poly_x[i:].replace(current_letter, letter, 1)
                all_polys.append(pattern)

    all_polys = set(all_polys)
    return all_polys


def _read_adapters(input_adapter: Path, *, use_set: bool = False) -> AdaptersNaive:
    """ Read adapters from a file and return them as list of str or set of str, as determined by `use_set`.

        Sets have faster lookup than lists, and thus may be preferred over them.
    """
    adapter_list: List[str] = []
    adapter_set: Set[str] = set()
    with open(input_adapter, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as adapter_handle:
        for line in adapter_handle:
            adapter = line.strip()
            if len(adapter) <= ADAPTER_LEN:
                adapter_list.append(adapter)
                adapter_set.add(adapter)
    adapters = adapter_set if use_set else adapter_list
    return adapters


def _filter_out_by_poly_x_naive(sequence: str, poly_patterns: Set[str]) -> bool:
    """Return True if the record should be filtered out (discarded), otherwise False. Naive implementation."""
    for pattern in poly_patterns:
        if pattern in sequence:
            return True
    return False


def _filter_out_by_adapters_naive(sequence: str, adapters: AdaptersNaive) -> bool:
    """Return True if the record should be filtered out (discarded), otherwise False. Naive implementation."""
    for adapter in adapters:
        if adapter in sequence:
            return True
    return False


def _filter_out_by_poly_x_trie(sequence: str, poly_patterns: Trie) -> bool:
    """Return True if the record should be filtered out (discarded), otherwise False. Uses Trie, `poly_patterns`."""
    return trie_matching(sequence, poly_patterns)


def _filter_out_by_adapters_trie(sequence: str, adapters: Trie) -> bool:
    """Return True if the record should be filtered out (discarded), otherwise False. Uses Trie, `adapters`."""
    return trie_matching(sequence, adapters)


@time_it  # 12 s. CORRECT, but don't use this. Modin: +Inf s
def pandas_iter_naive_gzip_counter(input_fastq: Path,
                                   output_fastq: Path,
                                   output_stat: Path,
                                   poly_patterns: Set[str],
                                   adapters: AdaptersNaive) -> None:
    """
    Low-level implementation of the main filtering logic.
    Uses Pandas DataFrame iteration: iterrows or itertuples. This is not recommended. This is antipattern in Pandas.
    Uses naive algorithms for pattern matching, for polyX and adapters.
    Uses simple loop counter.
    Can use sequential or Modin - see code.
    Uses Pandas for reading the input file, which means Modin can parallelize that part.
    Does NOT use Pandas for writing the output file, which means Modin cannot parallelize that part.
    Can use pgzip instead of gzip.
    """
    num_filtered_out_by_poly_x = 0
    num_filtered_out_by_adapters = 0

    # Sequential
    input_frame = pandas.DataFrame(
        pandas.read_csv(
            input_fastq, sep=PANDAS_SEPARATOR, header=None, encoding="ascii", encoding_errors="strict"
        ).values.reshape(-1, 4), columns=PANDAS_COLUMNS
    )

    # Modin:
    # input_frame = pd.DataFrame(
    #     pd.read_csv(
    #         input_fastq, sep=PANDAS_SEPARATOR, header=None, encoding="ascii", encoding_errors="strict"
    #         # input_fastq, header=None, encoding="ascii", encoding_errors="strict", delim_whitespace=True
    #     ).values.reshape(-1, 4), columns=PANDAS_COLUMNS
    # )

    input_frame = input_frame.reset_index(drop=True)
    # print(input_frame.info())
    # print(f"Pandas input frame shape = {input_frame.shape}")  # (100000, 4)
    # print(f"Pandas input frame size = {input_frame.size}")  # 400000
    # print("=" * 50)
    # print(input_frame.head())
    # print("=" * 50)
    # print(input_frame.head().to_string())
    # print("=" * 50)
    # exit_program("Testing...")

    output_list = []  # np.empty(), perhaps?

    # for row in input_frame.iterrows():  # 16 s
    for row in input_frame.itertuples():  # 10 s
        # print(len(row), row)
        # line_no = row[0]

        # read = row[1]
        read = row[1:]
        title, sequence, plus, quality = read

        # print(title, sequence, plus, quality)
        # print(row[0])
        # print(row[1])
        # print(row[1:])
        # print(row[1][0], row[1][1])
        # print(row.Index, row.read_id)
        # exit_program("Testing...")

        # Step 1: Filter by *poly-X*.
        is_filtered_out_by_poly_x = _filter_out_by_poly_x_naive(sequence, poly_patterns)
        if is_filtered_out_by_poly_x:
            num_filtered_out_by_poly_x += 1
            continue

        # Step 2: Filter by *adapters*.
        is_filtered_out_by_adapters = _filter_out_by_adapters_naive(sequence, adapters)
        if is_filtered_out_by_adapters:
            num_filtered_out_by_adapters += 1
            continue

        # Step 3: Write the record to the output file if not filtered out.
        output_list.append(f"{title}{NEWLINE}{sequence}{NEWLINE}+{NEWLINE}{quality}{NEWLINE}")
        # output_list.append(f"{title}{NEWLINE}{sequence}{NEWLINE}+{NEWLINE}{quality}")
        # output_list.append(f"{title}{sequence}+{quality}")
        # output_list.append(f"{title}{NEWLINE}")
        # output_list.append(f"{sequence}{NEWLINE}")
        # output_list.append(f"+{NEWLINE}")
        # output_list.append(f"{quality}{NEWLINE}")

    # print("\n\n\n")
    # print(len(output_list))  # 70058 or 280232, both correct, but use 70058
    # print("*" * 50)
    # print(output_list[0])
    # print("@" * 50)
    output_frame = pandas.DataFrame(output_list)#.values.reshape(-1, 4)
    # output_frame = pandas.DataFrame(output_list, columns=PANDAS_COLUMNS)#.values.reshape(-1, 4)
    output_frame.columns = ["read"]
    # print(output_frame.head())
    # print("#" * 50)
    # print(output_frame.head().to_string())
    # print("$" * 50)
    # print(output_frame.head().to_string().split("\n"))#[:5])
    # print("%" * 50)
    # print(output_frame.explode(["read"]))
    # print(output_frame.assign(read=output_frame["read"].str.split(NEWLINE)).explode(["read"]))
    # # output_frame.assign(var1=output_frame['var1'].str.split(',')).explode('var1')
    # print(output_frame.shape, output_frame.size, output_frame.ndim, output_frame.dtypes)  # (70058, 1) 70058 2 read    object
    # print(output_frame.info())
    # print("*" * 50)
    # print(output_frame)
    # print("*" * 50)
    # print(output_frame[0])
    # print("*" * 50)
    # print(output_frame[0][0])  # One record, the first one
    # output_frame.to_csv(output_fastq, header=False, index=False, sep=NEWLINE, quoting=csv.QUOTE_NONE, line_terminator=NEWLINE, escapechar=NEWLINE)
    # output_frame.to_csv(output_fastq, header=False, index=False, sep=NEWLINE, quoting=csv.QUOTE_ALL, line_terminator=NEWLINE, quotechar=NEWLINE)#, escapechar=NEWLINE)
    # output_frame.to_csv(output_fastq, header=False, index=False, sep=NEWLINE, quoting=csv.QUOTE_ALL, line_terminator=NEWLINE, quotechar="\\")
    # output_frame.to_csv(output_fastq, header=False, index=False, sep=NEWLINE, quoting=csv.QUOTE_NONE, line_terminator=NEWLINE, escapechar="\\")
    # print(f"Pandas output frame shape = {output_frame.shape}")  # (70058, 1)
    # print(f"Pandas output frame size = {output_frame.size}")  # 70058

    ## with open(OUTPUT_TEXT_FILE, "wt", newline=NEWLINE) as output_handle:
    ##     output_handle.write("".join(output_list))
    ## compress_file(OUTPUT_TEXT_FILE, output_fastq)

    with gzip.open(output_fastq, "wt", **OPEN_PARAMS) as dst_handle:  # Correct.
        contents = "".join(output_list)
        dst_handle.write(contents)

    # Step 4: Store the number of records filtered out by *poly-X* and by *adapters*, respectively.
    print(num_filtered_out_by_poly_x, num_filtered_out_by_adapters)  # 9567 20375, which is correct
    stats = f"filterByPolyX:\t{num_filtered_out_by_poly_x}{NEWLINE}" \
            f"filterByAdapter:\t{num_filtered_out_by_adapters}{NEWLINE}"
    with open(output_stat, "wt", newline=NEWLINE) as stat_handle:
        stat_handle.write(stats)


@time_it  # 12 s almost correct, don't use this (for now), rework the writing to file part
def pandas_iter_naive_gzip_counter2(input_fastq: Path,
                                    output_fastq: Path,
                                    output_stat: Path,
                                    poly_patterns: Set[str],
                                    adapters: AdaptersNaive) -> None:
    """
    Low-level implementation of the main filtering logic.
    Uses Numpy array iteration.
    Uses naive algorithms for pattern matching, for polyX and adapters.
    Uses simple loop counter.
    Can use sequential or Modin - see code.
    Uses Pandas for reading the input file, which means Modin can parallelize that part.
    Does NOT use Pandas for writing the output file, which means Modin cannot parallelize that part.
    Can use pgzip instead of gzip.
    """
    # TODO: Try with Modin; Try with Dask, separately of Modin. We have numpy.ndarray here.
    num_filtered_out_by_poly_x = 0
    num_filtered_out_by_adapters = 0

    input_text = pandas.read_csv(
        input_fastq, sep=PANDAS_SEPARATOR, header=None, encoding="ascii", encoding_errors="strict"
    ).values.reshape(-1, 4)
    # print(input_text)
    # print(type(input_text), input_text.shape, input_text.size, input_text.dtype, input_text.ndim)  # <class 'numpy.ndarray'> (100000, 4) 400000 object 2

    output_list = []

    for read_no, read in enumerate(input_text, 1):
        # print(read_no, read)
        title, sequence, plus, quality = read
        # print(title, sequence, plus, quality)

        # Step 1: Filter by *poly-X*.
        is_filtered_out_by_poly_x = _filter_out_by_poly_x_naive(sequence, poly_patterns)
        if is_filtered_out_by_poly_x:
            num_filtered_out_by_poly_x += 1
            continue

        # Step 2: Filter by *adapters*.
        is_filtered_out_by_adapters = _filter_out_by_adapters_naive(sequence, adapters)
        if is_filtered_out_by_adapters:
            num_filtered_out_by_adapters += 1
            continue

        # Step 3: Write the record to the output file if not filtered out.
        output_list.append(f"{title}{NEWLINE}{sequence}{NEWLINE}+{NEWLINE}{quality}{NEWLINE}")

    # TODO: Rework writing to file
    print(len(output_list))  # 70058, correct
    # print(output_list)
    # output_frame = pandas.DataFrame(output_list, columns=PANDAS_COLUMNS)#.values.reshape(-1, 4)
    output_frame = pandas.DataFrame(output_list)#.values.reshape(-1, 4)
    # print(output_frame)
    # print(output_frame[0])
    # print(output_frame[0][0])  # One record, the first one
    output_frame.to_csv(output_fastq, sep=NEWLINE)
    print(f"Pandas output frame shape = {output_frame.shape}")  # (70058, 1)
    print(f"Pandas output frame size = {output_frame.size}")  # 70058

    # Step 4: Store the number of records filtered out by *poly-X* and by *adapters*, respectively.
    print(num_filtered_out_by_poly_x, num_filtered_out_by_adapters)  # 9567 20375, which is correct
    stats = f"filterByPolyX:\t{num_filtered_out_by_poly_x}{NEWLINE}" \
            f"filterByAdapter:\t{num_filtered_out_by_adapters}{NEWLINE}"
    with open(output_stat, "wt", newline=NEWLINE) as stat_handle:
        stat_handle.write(stats)


@time_it  # CORRECT. Seq: 14 s, Modin: 16 s. Can use this.
def pandas_apply_naive(input_fastq: Path,
                       output_fastq: Path,
                       output_stat: Path,
                       poly_patterns: Set[str],
                       adapters: AdaptersNaive) -> None:
    """
    Low-level implementation of the main filtering logic.
    Uses Pandas Series apply. This is recommended in Pandas.
    Uses naive algorithms for pattern matching, for polyX and adapters.
    Can use sequential or Modin - see code.
    Uses Pandas for reading the input file, which means Modin can parallelize that part.
    Uses Pandas for writing the output file, which means Modin can parallelize that part.
    """
    # input_frame_seq = pandas.DataFrame(
    #     pandas.read_csv(
    #         input_fastq, sep=PANDAS_SEPARATOR, header=None
    #     ).values.reshape(-1, 4), columns=PANDAS_COLUMNS
    # )

    input_frame_modin = pd.DataFrame(
        pd.read_csv(
            input_fastq, sep=PANDAS_SEPARATOR, header=None
        ).values.reshape(-1, 4), columns=PANDAS_COLUMNS
    )

    input_frame = input_frame_modin

    # print(input_frame.info())
    # print(input_frame.columns)
    # print(input_frame.seq[:5])
    # print(f"Pandas input frame shape = {input_frame.shape}")  # (100000, 4)
    # print(f"Pandas input frame size = {input_frame.size}")  # 400000
    # print("\n\n")

    # Step 1: Filter by *poly-X*.
    are_filtered_out_by_poly_x = input_frame["seq"].apply(lambda seq: _filter_out_by_poly_x_naive(seq, poly_patterns))
    filtered_out_by_poly_x = are_filtered_out_by_poly_x[are_filtered_out_by_poly_x].index  # We need this.
    # kept_by_poly_x = are_filtered_out_by_poly_x[~are_filtered_out_by_poly_x].index  # Inverse of what we need.
    num_filtered_out_by_poly_x = are_filtered_out_by_poly_x.value_counts()[True]
    # print(input_frame.drop(filtered_out_by_poly_x))  # [90433 rows x 4 columns]
    # print()
    # print(are_filtered_out_by_poly_x)
    # print()
    # print(are_filtered_out_by_poly_x.info())  # Length: 100000
    # print()
    # print(are_filtered_out_by_poly_x.count())  # 100000
    # print()
    # print(are_filtered_out_by_poly_x.values)
    # print()
    # print(are_filtered_out_by_poly_x.value_counts())  # False 90433, True 9567
    # print()
    # print(are_filtered_out_by_poly_x.value_counts()[True])  # 9567
    # print(are_filtered_out_by_poly_x.value_counts()[False])  # 90433
    # print()
    # print(type(are_filtered_out_by_poly_x))  # <class 'pandas.core.series.Series'>
    # print()
    # print(filtered_out_by_poly_x)  # length=9567
    # print()
    # print(input_frame.loc[filtered_out_by_poly_x])  # [9567 rows x 4 columns]
    # print(input_frame.loc[kept_by_poly_x])  # [90433 rows x 4 columns]
    # exit_program("Testing...")

    # We can work with one data frame only, *input_frame*, modifying it, but we'll create a new data frame.
    output_frame = input_frame.drop(filtered_out_by_poly_x)

    # Step 2: Filter by *adapters*.
    are_filtered_out_by_adapters = output_frame["seq"].apply(lambda seq: _filter_out_by_adapters_naive(seq, adapters))
    filtered_out_by_adapters = are_filtered_out_by_adapters[are_filtered_out_by_adapters].index
    num_filtered_out_by_adapters = are_filtered_out_by_adapters.value_counts()[True]
    # print(are_filtered_out_by_adapters)
    # print()
    # print(are_filtered_out_by_adapters.info())  # Length: 90433
    # print()
    # print(are_filtered_out_by_adapters.count())  # 90433
    # print()
    # print(are_filtered_out_by_adapters.values)
    # print()
    # print(are_filtered_out_by_adapters.value_counts())  # False 70058, True 20375
    # print()
    # print(are_filtered_out_by_adapters.value_counts()[True])  # 20375
    # print(are_filtered_out_by_adapters.value_counts()[False])  # 70058
    # print()
    # print(type(are_filtered_out_by_adapters))  # <class 'pandas.core.series.Series'>
    # print()
    # print(filtered_out_by_adapters)  # length=20375
    # print()
    # print(output_frame.loc[filtered_out_by_adapters])  # [20375 rows x 4 columns]
    # exit_program("Testing...")

    output_frame.drop(filtered_out_by_adapters, inplace=True)

    # print(output_frame.head(6).to_string())
    # print(output_frame.tail().to_string())
    # exit_program("Testing...")

    # Step 3: Write the record to the output file if not filtered out.
    output_frame.to_csv(
        output_fastq, header=False, index=False, sep=NEWLINE,
        quoting=csv.QUOTE_NONE, line_terminator=NEWLINE, escapechar=NEWLINE
    )

    # print(f"Pandas output frame shape = {output_frame.shape}")  # (70058, 4)
    # print(f"Pandas output frame size = {output_frame.size}")  # 280232

    # Step 4: Store the number of records filtered out by *poly-X* and by *adapters*, respectively.
    print(num_filtered_out_by_poly_x, num_filtered_out_by_adapters)  # 9567 20375, which is correct
    stats = f"filterByPolyX:\t{num_filtered_out_by_poly_x}{NEWLINE}" \
            f"filterByAdapter:\t{num_filtered_out_by_adapters}{NEWLINE}"
    with open(output_stat, "wt", newline=NEWLINE) as stat_handle:
        stat_handle.write(stats)


@time_it  # CORRECT. Seq: 45 s, Modin: 30 s. Can use this.
def pandas_apply_trie(input_fastq: Path,
                      output_fastq: Path,
                      output_stat: Path,
                      poly_patterns: Trie,
                      adapters: Trie) -> None:
    """
    Low-level implementation of the main filtering logic.
    Uses Pandas Series apply. This is recommended in Pandas.
    Uses Trie algorithms for pattern matching, for polyX and adapters.
    Can use sequential or Modin - see code.
    Uses Pandas for reading the input file, which means Modin can parallelize that part.
    Uses Pandas for writing the output file, which means Modin can parallelize that part.
    """
    # input_frame_seq = pandas.DataFrame(
    #     pandas.read_csv(
    #         input_fastq, sep=PANDAS_SEPARATOR, header=None
    #     ).values.reshape(-1, 4), columns=PANDAS_COLUMNS
    # )

    input_frame_modin = pd.DataFrame(
        pd.read_csv(
            input_fastq, sep=PANDAS_SEPARATOR, header=None
        ).values.reshape(-1, 4), columns=PANDAS_COLUMNS
    )

    input_frame = input_frame_modin

    # Step 1: Filter by *poly-X*.
    are_filtered_out_by_poly_x = input_frame["seq"].apply(lambda seq: _filter_out_by_poly_x_trie(seq, poly_patterns))
    filtered_out_by_poly_x = are_filtered_out_by_poly_x[are_filtered_out_by_poly_x].index  # We need this.
    num_filtered_out_by_poly_x = are_filtered_out_by_poly_x.value_counts()[True]

    # We can work with one data frame only, *input_frame*, modifying it, but we'll create a new data frame.
    output_frame = input_frame.drop(filtered_out_by_poly_x)

    # Step 2: Filter by *adapters*.
    are_filtered_out_by_adapters = output_frame["seq"].apply(lambda seq: _filter_out_by_adapters_trie(seq, adapters))
    filtered_out_by_adapters = are_filtered_out_by_adapters[are_filtered_out_by_adapters].index  # We need this.
    num_filtered_out_by_adapters = are_filtered_out_by_adapters.value_counts()[True]

    output_frame.drop(filtered_out_by_adapters, inplace=True)

    # Step 3: Write the record to the output file if not filtered out.
    output_frame.to_csv(
        output_fastq, header=False, index=False, sep=NEWLINE,
        quoting=csv.QUOTE_NONE, line_terminator=NEWLINE, escapechar=NEWLINE
    )

    # Step 4: Store the number of records filtered out by *poly-X* and by *adapters*, respectively.
    print(num_filtered_out_by_poly_x, num_filtered_out_by_adapters)  # 9567 20375, which is correct
    stats = f"filterByPolyX:\t{num_filtered_out_by_poly_x}{NEWLINE}" \
            f"filterByAdapter:\t{num_filtered_out_by_adapters}{NEWLINE}"
    with open(output_stat, "wt", newline=NEWLINE) as stat_handle:
        stat_handle.write(stats)


@time_it
def main_logic_naive(input_fastq: Path, input_adapter: Path, output_fastq: Path, output_stat: Path) -> None:
    """High-level implementation of the main filtering logic. Naive implementation. Sequential or parallel."""
    all_polyx_patterns = _generate_all_polyx_patterns()
    adapters = _read_adapters(input_adapter, use_set=True)
    # pandas_iter_naive_gzip_counter(input_fastq, output_fastq, output_stat, all_polyx_patterns, adapters)  # 12 s
    # pandas_iter_naive_gzip_counter2(input_fastq, output_fastq, output_stat, all_polyx_patterns, adapters)  # 15 s
    pandas_apply_naive(input_fastq, output_fastq, output_stat, all_polyx_patterns, adapters)  # 14/16 s


@time_it
def main_logic_trie(input_fastq: Path, input_adapter: Path, output_fastq: Path, output_stat: Path) -> None:
    """High-level implementation of the main filtering logic. Separate Tries. Sequential or parallel."""
    all_polyx_patterns = _generate_all_polyx_patterns()
    adapters = _read_adapters(input_adapter, use_set=True)
    polyx_patterns_trie = build_trie(all_polyx_patterns)
    adapters_trie = build_trie(adapters)
    pandas_apply_trie(input_fastq, output_fastq, output_stat, polyx_patterns_trie, adapters_trie)  # 45/30 s


@time_it
def _validate_filtering() -> None:
    """
    Verify that reference and this solution's "out.fq.gz" files are the same.
    Verify that reference and this solution's "out.stat.txt" files are the same.
    These validations are the same as unit tests in "test_main.py".
    """
    print("\n\n VALIDATION \n")

    print("Comparing \"out.fq.gz\" files...")
    with gzip.open(OUTPUT_FASTQ_GZ, "rt", **OPEN_PARAMS) as solution_handle:
        solution_contents = solution_handle.read()
    with gzip.open(TEST_OUT_FASTQ_GZ_REFERENCE, "rt", **OPEN_PARAMS) as reference_handle:
        reference_contents = reference_handle.read()
    print(len(reference_contents), len(solution_contents))
    assert len(reference_contents) == len(solution_contents)
    assert reference_contents == solution_contents

    print("\nComparing \"out.stat.txt\" files...")
    with open(OUTPUT_STATISTICS, "rt", **OPEN_PARAMS) as solution_handle:
        solution_contents = solution_handle.read()
    with open(TEST_OUT_STAT_REFERENCE, "rt", **OPEN_PARAMS) as reference_handle:
        reference_contents = reference_handle.read()
    print(len(reference_contents), len(solution_contents))
    print(solution_contents)
    assert len(reference_contents) == len(solution_contents)
    assert reference_contents == solution_contents


def _validate_pgzip_decompresses_output_file(output_fastq: Path) -> None:
    with pgzip.open(output_fastq, "rt", thread=0) as solution_handle:
        solution_length = len(solution_handle.read())
        print(f"'{output_fastq}' file length is: {solution_length} bytes.")
        assert TEST_OUT_FASTQ_SIZE_REF == solution_length


@time_it  # 1 s
def try_out_pandas_1():
    # Pandas do read gzip files.
    frame = pandas.DataFrame(
        pandas.read_csv(
            INPUT_FASTQ, sep="\r", header=None
        ).values.reshape(-1, 4), columns=["read_id", "seq", "+", "qual"]
    )  # Use "\r" as separator on Windows, and "\n" elsewhere.
    print(f"Pandas frame shape = {frame.shape}")  # (100000, 4)
    print(f"Pandas frame size = {frame.size}")  # 400000
    print(frame.head())
    print("=" * 120)
    print("=" * 120, "\n")

    # # Dask doesn't read gzip files.
    # with pgzip.open(INPUT_FASTQ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
    #     frame = pandas.DataFrame(
    #         pandas.read_csv(input_handle, sep="\r", header=None).values.reshape(-1, 4),
    #         columns=["read_id", "seq", "+", "qual"]
    #     )  # Use "\r" as separator on Windows, and "\n" elsewhere.
    # print(frame.head())
    #
    # print("=" * 120)
    # print("=" * 120)
    # print("=" * 120, "\n")


@time_it  # 3 s
def try_out_pandas_modin_1():
    # Pandas do read gzip files.
    frame = pd.DataFrame(
        pd.read_csv(
            INPUT_FASTQ, sep=PANDAS_SEPARATOR, header=None
        ).values.reshape(-1, 4), columns=["read_id", "seq", "+", "qual"]
    )  # Use "\r" as separator on Windows, and "\n" elsewhere.
    print(frame.head())
    print("=" * 120)
    print("=" * 120, "\n")


def _initialize_dask():
    client = Client()


@time_it
def try_out_dask_1():
    frame = ddf.read_csv(TEST_INP_FASTQ_SMALL_ORIGINAL, blocksize=6400000, lineterminator="\r", header=None)
    print(frame.head())
    print("="*120, "\n")
    print(frame.loc[0, :])
    print("="*120, "\n")
    for row in frame.loc[0, :]:
        print(row)
    print("="*120)
    print("="*120, "\n")

    frame = ddf.read_csv(TEST_INP_FASTQ_SMALL_ORIGINAL, blocksize=None, lineterminator="\r", dtype=np.int8)
    print(frame.head(5))
    print(frame.index)
    print(frame.shape)
    print(frame.size)
    print(frame.info())
    print("=" * 120)
    print("="*120, "\n")

    # frame = dask.delayed(pd.read_csv)(INPUT_FASTQ).compute()
    # print(frame.head(5))
    # # for item in frame:
    # #     print(item)
    # print("=" * 120)
    # print("=" * 120, "\n")

    # with gzip.open(INPUT_FASTQ, "rt", **OPEN_PARAMS) as input_handle:
    # # with gzip.open(INPUT_FASTQ, "rb") as input_handle:
    #     frame = ddf.read_csv(input_handle, blocksize=6400000, lineterminator="\r")
    # print(frame.head(5))

    print("=" * 120)
    print("=" * 120)
    print("=" * 120, "\n")


if __name__ == "__main__":
    _initialize_dask()

    modin.config.Engine.put("Dask")  # noqa
    print(f"Modin num partitions = {modin.config.NPartitions.get()}")  # noqa

    # main_logic_naive(INPUT_FASTQ, INPUT_ADAPTER, OUTPUT_FASTQ_GZ, OUTPUT_STATISTICS)
    main_logic_trie(INPUT_FASTQ, INPUT_ADAPTER, OUTPUT_FASTQ_GZ, OUTPUT_STATISTICS)
    _validate_filtering()
    _validate_pgzip_decompresses_output_file(OUTPUT_FASTQ_GZ)

    # try_out_pandas_1()
    # try_out_pandas_modin_1()
    # try_out_dask_1()
