#!/usr/bin/env python
"""
File:       exp/exp_parallel.py
Author:     Ivan Lazarević
Brief:      Script for experimenting with parallelization of the main logic using standard library.
"""
# TODO: Nothing implemented

# Standard library imports
import concurrent.futures
import gzip
import multiprocessing
import os
import sys
from itertools import zip_longest
from pathlib import Path
from typing import List, Set, Tuple, Union

# Third party library imports
import numpy as np
import pgzip

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import ALPHABET, NEWLINE, ADAPTER_LEN, POLY_LEN, POLY_END, ADAPTER_END, OPEN_PARAMS
from src.config import INPUT_FASTQ, INPUT_ADAPTER
from src.config import OUTPUT_FASTQ_GZ, TEST_OUT_FASTQ_SIZE_REF, OUTPUT_STATISTICS
from src.config import TEST_OUT_FASTQ_GZ_REFERENCE, TEST_OUT_STAT_REFERENCE
from src.config import TEST_INP_FASTQ_SMALL_ORIGINAL
from src.utils import exit_program, time_it
from src.trie import build_trie, trie_matching, trie_matching_combined
from src.type_aliases import Trie, AdaptersNaive


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


@time_it  # 11 s, gzip or pgzip
def _worker_seq_naive_pgzip_counter(input_fastq: Path,
                                    output_fastq: Path,
                                    output_stat: Path,
                                    poly_patterns: Set[str],
                                    adapters: AdaptersNaive) -> None:
    """Low-level implementation of the main filtering logic. Sequential. No *Biopython* at all. Simple loop counter."""
    num_filtered_out_by_poly_x = 0
    num_filtered_out_by_adapters = 0

    with gzip.open(input_fastq, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with gzip.open(output_fastq, "wt", encoding="ascii", errors="strict", newline=NEWLINE) as output_handle:
            keep_record = True
            for line_no, line in enumerate(input_handle, 1):
                if line_no % 4 == 1:
                    title = line[:-1]

                elif line_no % 4 == 2:
                    sequence = line[:-1]

                    # Step 1: Filter by *poly-X*.
                    is_filtered_out_by_poly_x = _filter_out_by_poly_x_naive(sequence, poly_patterns)
                    if is_filtered_out_by_poly_x:
                        num_filtered_out_by_poly_x += 1
                        keep_record = False
                        continue

                    # Step 2: Filter by *adapters*.
                    is_filtered_out_by_adapters = _filter_out_by_adapters_naive(sequence, adapters)
                    if is_filtered_out_by_adapters:
                        num_filtered_out_by_adapters += 1
                        keep_record = False
                        continue

                # Step 3: Write the record to the output file if not filtered out.
                elif line_no % 4 == 0:
                    if keep_record:
                        quality = line[:-1]
                        output_handle.write(f"{title}{NEWLINE}{sequence}{NEWLINE}+{NEWLINE}{quality}{NEWLINE}")
                    keep_record = True

    # Step 4: Store the number of records filtered out by *poly-X* and by *adapters*, respectively.
    print(num_filtered_out_by_poly_x, num_filtered_out_by_adapters)
    stats = f"filterByPolyX:\t{num_filtered_out_by_poly_x}{NEWLINE}" \
            f"filterByAdapter:\t{num_filtered_out_by_adapters}{NEWLINE}"
    with open(output_stat, "wt", newline=NEWLINE) as stat_handle:
        stat_handle.write(stats)


@time_it  # 11 s
def _worker_seq_naive_pgzip_zip(input_fastq: Path,
                                output_fastq: Path,
                                output_stat: Path,
                                poly_patterns: Set[str],
                                adapters: AdaptersNaive) -> None:
    """Low-level implementation of the main filtering logic. Sequential. No *Biopython* at all. Uses *itertools*."""
    num_filtered_out_by_poly_x = 0
    num_filtered_out_by_adapters = 0

    with gzip.open(input_fastq, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with gzip.open(output_fastq, "wt", encoding="ascii", errors="strict", newline=NEWLINE) as output_handle:
            fastq_iterator = (line[:-1] for line in input_handle)
            for record in zip_longest(*[fastq_iterator] * 4):
                title, sequence, quality = record[0], record[1], record[3]

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

                output_handle.write(f"{title}{NEWLINE}{sequence}{NEWLINE}+{NEWLINE}{quality}{NEWLINE}")

    # Step 4: Store the number of records filtered out by *poly-X* and by *adapters*, respectively.
    print(num_filtered_out_by_poly_x, num_filtered_out_by_adapters)
    stats = f"filterByPolyX:\t{num_filtered_out_by_poly_x}{NEWLINE}" \
            f"filterByAdapter:\t{num_filtered_out_by_adapters}{NEWLINE}"
    with open(output_stat, "wt", newline=NEWLINE) as stat_handle:
        stat_handle.write(stats)


@time_it
def main_logic_seq_naive(input_fastq: Path, input_adapter: Path, output_fastq: Path, output_stat: Path) -> None:
    """High-level implementation of the main filtering logic. Naive implementation. Sequential."""
    all_polyx_patterns = _generate_all_polyx_patterns()
    adapters = _read_adapters(input_adapter, use_set=True)
    _worker_seq_naive_pgzip_zip(input_fastq, output_fastq, output_stat, all_polyx_patterns, adapters)  # 11 s


def mp_pool():
    with multiprocessing.Pool() as pool:
        raise NotImplementedError
        # pool.map(_worker_seq_naive_pgzip_zip)


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


if __name__ == "__main__":
    # main_logic_par_naive(INPUT_FASTQ, INPUT_ADAPTER, OUTPUT_FASTQ_GZ, OUTPUT_STATISTICS)

    num_cpus = multiprocessing.cpu_count()
    print(num_cpus)

    mp_pool()
    _validate_filtering()
    _validate_pgzip_decompresses_output_file(OUTPUT_FASTQ_GZ)
