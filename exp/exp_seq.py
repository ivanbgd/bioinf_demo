#!/usr/bin/env python
"""
File:       exp/exp_seq.py
Author:     Ivan Lazarević
Brief:      Script for experimenting with reading, parsing and writing gzipped FASTQ files,
            filtering by the assignment text (the task).

Details:    We vary algorithms for `filter_out_by_poly_x_*` and for `filter_out_by_adapters_*`,
            starting with naive algorithms and moving on to more advanced implementations, such as Trie.

            We vary `worker_*` by changing the way we read input data, process it, and write it.
            We do this sequentially in the beginning, and later in parallel, in "exp/exp_parallel.py".
"""
# Standard library imports
import gzip
import os
import sys
from itertools import zip_longest
from pathlib import Path
from typing import List, Set, Tuple, Union

# Third party library imports
import pgzip
from Bio import Seq, SeqIO, bgzf
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import ALPHABET, NEWLINE, ADAPTER_LEN, POLY_LEN, POLY_END, ADAPTER_END, OPEN_PARAMS
from src.config import INPUT_FASTQ, INPUT_ADAPTER
from src.config import OUTPUT_FASTQ_GZ, TEST_OUT_FASTQ_SIZE_REF, OUTPUT_STATISTICS
from src.config import TEST_OUT_FASTQ_GZ_REFERENCE, TEST_OUT_STAT_REFERENCE
from src.trie import build_trie, trie_matching, trie_matching_combined
from src.trie import build_trie_improved, trie_matching_improved
from src.type_aliases import Trie, AdaptersNaive
from src.utils import exit_program, time_it


num_filtered_out_by_poly_x_global = 0
num_filtered_out_by_adapters_global = 0


# Unused.
def _check_record_naive(sequence: str, polys: Set[str], adapters: AdaptersNaive) -> Tuple[bool, bool]:
    """Check a single record"""
    # Step 1: Filter by *poly-X*.
    is_filtered_out_by_poly_x = _filter_out_by_poly_x_naive(sequence, polys)

    # Step 2: Filter by *adapters*.
    is_filtered_out_by_adapters = _filter_out_by_adapters_naive(sequence, adapters)

    return is_filtered_out_by_poly_x, is_filtered_out_by_adapters


def _test_generating_poly_patterns(all_polys: Set[str]) -> None:
    test_patterns = ["CAAAAAAAAAAAAAA", "GAAAAAAAAAAAAAA",
                     "TAAAAAAAAAAAAAA", "CTTTTTTTTTTTTTT",
                     "AAAAAAAAAAAAAAG", "TTTTTTTTTTTTTTA"]
    expected_len = len(ALPHABET) * (1 + POLY_LEN * (len(ALPHABET) - 1))
    print("Testing _test_generating_poly_patterns...")
    sorted_patterns = sorted(all_polys)
    for pattern in sorted_patterns:
        print(pattern, len(pattern))
        assert POLY_LEN == len(pattern)
    print(len(all_polys), expected_len)
    for pattern in test_patterns:
        print(pattern in all_polys)
    assert expected_len == len(all_polys)
    exit_program("_test_generating_poly_patterns")


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
            # idx = poly_x.find(current_letter, i)
            # print(x, current_letter, i, poly_x, idx)
            for letter in ALPHABET:
                pattern = poly_x[:i] + poly_x[i:].replace(current_letter, letter, 1)
                all_polys.append(pattern)

    all_polys = set(all_polys)
    # _test_generating_poly_patterns(all_polys)
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


def _build_combined_trie(input_adapter: Path) -> Trie:
    """ Build a combined trie from both *poly-X* and *adapter* sequences and return it

        The trie looks like:
        {0: {'C': 1, 'A': 47, 'G': 63, 'T': 148}, 1: {'C': 2, 'G': 103, 'T': 278, 'A': 473}, ...
        It contains POLY_END ADAPTER_END nodes. Leaves are empty nodes, as before.
        Its size is 1712 nodes.
    """
    adapters = _read_adapters(input_adapter, use_set=True)
    adapters_new = set()
    for adapter in adapters:
        adapters_new.add(adapter + ADAPTER_END)

    all_polyx_patterns = _generate_all_polyx_patterns()
    polys_new = set()
    for poly in all_polyx_patterns:
        polys_new.add(poly + POLY_END)

    combined_patterns = polys_new.union(adapters_new)
    combined_trie = build_trie(combined_patterns)
    print(len(combined_trie), combined_trie)
    return combined_trie


def _filter_out_by_poly_x_naive(sequence: str, poly_patterns: Set[str]) -> bool:
    """Return True if the record should be filtered out (discarded), otherwise False. Naive implementation."""
    for pattern in poly_patterns:
        if pattern in sequence:
            return True
    return False


def _filter_out_by_poly_x_trie(sequence: str, poly_patterns: Trie) -> bool:
    """Return True if the record should be filtered out (discarded), otherwise False. Uses Trie, `poly_patterns`."""
    return trie_matching(sequence, poly_patterns)


def _filter_out_by_poly_x_trie_improved(sequence: str, poly_patterns: Trie) -> bool:
    """
    Return True if the record should be filtered out (discarded), otherwise False. Uses Trie, `poly_patterns`.
    Can be used in case a pattern can be a prefix of another pattern.
    """
    return trie_matching_improved(sequence, poly_patterns)


def _filter_out_by_adapters_naive(sequence: str, adapters: AdaptersNaive) -> bool:
    """Return True if the record should be filtered out (discarded), otherwise False. Naive implementation."""
    for adapter in adapters:
        if adapter in sequence:
            return True
    return False


def _filter_out_by_adapters_trie(sequence: str, adapters: Trie) -> bool:
    """Return True if the record should be filtered out (discarded), otherwise False. Uses Trie, `adapters`."""
    return trie_matching(sequence, adapters)


def _filter_out_by_adapters_trie_improved(sequence: str, adapters: Trie) -> bool:
    """
    Return True if the record should be filtered out (discarded), otherwise False. Uses Trie, `adapters`.
    Can be used in case a pattern can be a prefix of another pattern.
    """
    return trie_matching_improved(sequence, adapters)


@time_it  # 40 s
def _worker_seq_naive_bp(input_fastq: Path,
                         # input_adapter: Path,
                         output_fastq: Path,
                         output_stat: Path,
                         poly_patterns: Set[str],
                         adapters: AdaptersNaive) -> None:
    """Low-level implementation of the main filtering logic. Sequential. Uses *Biopython*."""
    num_filtered_out_by_poly_x = 0
    num_filtered_out_by_adapters = 0

    with gzip.open(input_fastq, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with bgzf.BgzfWriter(output_fastq, "wb") as output_handle:
            for title, sequence, quality in FastqGeneralIterator(input_handle):
                record = SeqIO.SeqRecord(
                    id=title,
                    seq=Seq.Seq(sequence),
                    description=title,
                    letter_annotations={"phred_quality": [ord(q) - 33 for q in quality]},
                )

                # Step 1: Filter by *poly-X*.
                is_filtered_out_by_poly_x = _filter_out_by_poly_x_naive(record.seq, poly_patterns)
                if is_filtered_out_by_poly_x:
                    num_filtered_out_by_poly_x += 1
                    # break
                    continue

                # Step 2: Filter by *adapters*.
                is_filtered_out_by_adapters = _filter_out_by_adapters_naive(record.seq, adapters)
                if is_filtered_out_by_adapters:
                    num_filtered_out_by_adapters += 1
                    # break
                    continue

                # Step 3: Write the record to the output file if not filtered out.
                SeqIO.write(sequences=record, handle=output_handle, format="fastq")

    # Step 4: Store the number of records filtered out by *poly-X* and by *adapters*, respectively.
    print(num_filtered_out_by_poly_x, num_filtered_out_by_adapters)
    stats = f"filterByPolyX:\t{num_filtered_out_by_poly_x}{NEWLINE}" \
            f"filterByAdapter:\t{num_filtered_out_by_adapters}{NEWLINE}"
    with open(output_stat, "wt", newline=NEWLINE) as stat_handle:
        stat_handle.write(stats)


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


@time_it  # 80 s
def _worker_seq_trie_bp(input_fastq: Path,
                        output_fastq: Path,
                        output_stat: Path,
                        poly_patterns: Trie,
                        adapters: Trie) -> None:
    """Low-level implementation of the main filtering logic. Sequential. Uses *Biopython*."""
    num_filtered_out_by_poly_x = 0
    num_filtered_out_by_adapters = 0

    with gzip.open(input_fastq, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with bgzf.BgzfWriter(output_fastq, "wb") as output_handle:
            for title, sequence, quality in FastqGeneralIterator(input_handle):
                record = SeqIO.SeqRecord(
                    id=title,
                    seq=Seq.Seq(sequence),
                    description=title,
                    letter_annotations={"phred_quality": [ord(q) - 33 for q in quality]},
                )

                # Step 1: Filter by *poly-X*.
                is_filtered_out_by_poly_x = _filter_out_by_poly_x_trie(record.seq, poly_patterns)  # 80 s
                # is_filtered_out_by_poly_x = trie_matching(record.seq, poly_patterns)  # 80 s
                if is_filtered_out_by_poly_x:
                    num_filtered_out_by_poly_x += 1
                    # break
                    continue

                # Step 2: Filter by *adapters*.
                is_filtered_out_by_adapters = _filter_out_by_adapters_trie(record.seq, adapters)  # 80 s
                # is_filtered_out_by_adapters = trie_matching(record.seq, adapters)  # 80 s
                if is_filtered_out_by_adapters:
                    num_filtered_out_by_adapters += 1
                    # break
                    continue

                # Step 3: Write the record to the output file if not filtered out.
                SeqIO.write(sequences=record, handle=output_handle, format="fastq")

    # Step 4: Store the number of records filtered out by *poly-X* and by *adapters*, respectively.
    print(num_filtered_out_by_poly_x, num_filtered_out_by_adapters)
    stats = f"filterByPolyX:\t{num_filtered_out_by_poly_x}{NEWLINE}" \
            f"filterByAdapter:\t{num_filtered_out_by_adapters}{NEWLINE}"
    with open(output_stat, "wt", newline=NEWLINE) as stat_handle:
        stat_handle.write(stats)


def _filtering_logic_generator(records, poly_patterns, adapters):
    global num_filtered_out_by_poly_x_global, num_filtered_out_by_adapters_global
    num_filtered_out_by_poly_x_global = 0
    num_filtered_out_by_adapters_global = 0

    for record in records:
        # Step 1: Filter by *poly-X*.
        is_filtered_out_by_poly_x = _filter_out_by_poly_x_trie(record.seq, poly_patterns)
        if is_filtered_out_by_poly_x:
            num_filtered_out_by_poly_x_global += 1
            continue

        # Step 2: Filter by *adapters*.
        is_filtered_out_by_adapters = _filter_out_by_adapters_trie(record.seq, adapters)
        if is_filtered_out_by_adapters:
            num_filtered_out_by_adapters_global += 1
            continue

        yield record


@time_it  # 80 s
def _worker_seq_trie_bp_generator(input_fastq: Path,
                                  output_fastq: Path,
                                  output_stat: Path,
                                  poly_patterns: Trie,
                                  adapters: Trie) -> None:
    """Low-level implementation of the main filtering logic. Sequential. Uses Biopython. Writes all records at once."""
    with gzip.open(input_fastq, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with bgzf.BgzfWriter(output_fastq, "wb") as output_handle:
            # Steps 1 & 2: Filter by *poly-X* or by *adapters*.
            fastq_parser = SeqIO.parse(input_handle, format="fastq")

            # Step 3: Write all preserved records at once to the output file.
            SeqIO.write(sequences=_filtering_logic_generator(fastq_parser, poly_patterns, adapters),
                        handle=output_handle,
                        format="fastq")

    # Step 4: Store the number of records filtered out by *poly-X* and by *adapters*, respectively.
    print(num_filtered_out_by_poly_x_global, num_filtered_out_by_adapters_global)
    stats = f"filterByPolyX:\t{num_filtered_out_by_poly_x_global}{NEWLINE}" \
            f"filterByAdapter:\t{num_filtered_out_by_adapters_global}{NEWLINE}"
    with open(output_stat, "wt", newline=NEWLINE) as stat_handle:
        stat_handle.write(stats)


@time_it  # 55 s
def _worker_seq_trie_combined_bp(input_fastq: Path,
                                 output_fastq: Path,
                                 output_stat: Path,
                                 combined_trie: Trie) -> None:
    """Low-level implementation of the main filtering logic. Sequential. Uses *Biopython*."""
    # TODO: Almost fully correct. Namely, the output file contents are perfectly correct.
    # TODO: Only the two counters are off by 1, but their total is correct.
    num_filtered_out_by_poly_x = 0
    num_filtered_out_by_adapters = 0

    with gzip.open(input_fastq, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with bgzf.BgzfWriter(output_fastq, "wb") as output_handle:
            for title, sequence, quality in FastqGeneralIterator(input_handle):
                record = SeqIO.SeqRecord(
                    id=title,
                    seq=Seq.Seq(sequence),
                    description=title,
                    letter_annotations={"phred_quality": [ord(q) - 33 for q in quality]},
                )

                # Steps 1 & 2: Filter by *poly-X* and by *adapters* at the same time.
                result = trie_matching_combined(record.seq, combined_trie)

                # Step 1: Filter by *poly-X*.
                if result == POLY_END:
                    num_filtered_out_by_poly_x += 1
                # Step 2: Filter by *adapters*.
                elif result == ADAPTER_END:
                    num_filtered_out_by_adapters += 1
                # Step 3: Write the record to the output file if not filtered out.
                else:
                    SeqIO.write(sequences=record, handle=output_handle, format="fastq")

    # Step 4: Store the number of records filtered out by *poly-X* and by *adapters*, respectively.
    print(num_filtered_out_by_poly_x, num_filtered_out_by_adapters)
    stats = f"filterByPolyX:\t{num_filtered_out_by_poly_x}{NEWLINE}" \
            f"filterByAdapter:\t{num_filtered_out_by_adapters}{NEWLINE}"
    with open(output_stat, "wt", newline=NEWLINE) as stat_handle:
        stat_handle.write(stats)


@time_it  # 35 s, gzip or pgzip
def _worker_seq_trie_pgzip_counter(input_fastq: Path,
                                   output_fastq: Path,
                                   output_stat: Path,
                                   poly_patterns: Trie,
                                   adapters: Trie) -> None:
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
                    is_filtered_out_by_poly_x = _filter_out_by_poly_x_trie(sequence, poly_patterns)
                    if is_filtered_out_by_poly_x:
                        num_filtered_out_by_poly_x += 1
                        keep_record = False
                        continue

                    # Step 2: Filter by *adapters*.
                    is_filtered_out_by_adapters = _filter_out_by_adapters_trie(sequence, adapters)
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


@time_it  # 38 s, gzip or pgzip
def _worker_improved_seq_trie_pgzip_counter(input_fastq: Path,
                                            output_fastq: Path,
                                            output_stat: Path,
                                            poly_patterns: Trie,
                                            adapters: Trie) -> None:
    """
    Low-level implementation of the main filtering logic. Sequential. No *Biopython* at all. Simple loop counter.
    Can be used in case a pattern can be a prefix of another pattern.
    """
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
                    is_filtered_out_by_poly_x = _filter_out_by_poly_x_trie_improved(sequence, poly_patterns)
                    if is_filtered_out_by_poly_x:
                        num_filtered_out_by_poly_x += 1
                        keep_record = False
                        continue

                    # Step 2: Filter by *adapters*.
                    is_filtered_out_by_adapters = _filter_out_by_adapters_trie_improved(sequence, adapters)
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
    print(num_filtered_out_by_poly_x, num_filtered_out_by_adapters)  # 9567 0
    stats = f"filterByPolyX:\t{num_filtered_out_by_poly_x}{NEWLINE}" \
            f"filterByAdapter:\t{num_filtered_out_by_adapters}{NEWLINE}"
    with open(output_stat, "wt", newline=NEWLINE) as stat_handle:
        stat_handle.write(stats)


@time_it  # 35 s, gzip or pgzip
def _worker_seq_trie_pgzip_zip(input_fastq: Path,
                               output_fastq: Path,
                               output_stat: Path,
                               poly_patterns: Trie,
                               adapters: Trie) -> None:
    """Low-level implementation of the main filtering logic. Sequential. No *Biopython* at all. Uses *itertools*."""
    num_filtered_out_by_poly_x = 0
    num_filtered_out_by_adapters = 0

    with gzip.open(input_fastq, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with gzip.open(output_fastq, "wt", encoding="ascii", errors="strict", newline=NEWLINE) as output_handle:
            fastq_iterator = (line[:-1] for line in input_handle)
            for record in zip_longest(*[fastq_iterator] * 4):
                title, sequence, quality = record[0], record[1], record[3]

                # Step 1: Filter by *poly-X*.
                is_filtered_out_by_poly_x = _filter_out_by_poly_x_trie(sequence, poly_patterns)
                if is_filtered_out_by_poly_x:
                    num_filtered_out_by_poly_x += 1
                    continue

                # Step 2: Filter by *adapters*.
                is_filtered_out_by_adapters = _filter_out_by_adapters_trie(sequence, adapters)
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


@time_it  # 40 s
def _worker_improved_seq_trie_pgzip_zip(input_fastq: Path,
                                        output_fastq: Path,
                                        output_stat: Path,
                                        poly_patterns: Trie,
                                        adapters: Trie) -> None:
    """
    Low-level implementation of the main filtering logic. Sequential. No *Biopython* at all. Uses *itertools*.
    Can be used in case a pattern can be a prefix of another pattern.
    """
    num_filtered_out_by_poly_x = 0
    num_filtered_out_by_adapters = 0

    with gzip.open(input_fastq, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with gzip.open(output_fastq, "wt", encoding="ascii", errors="strict", newline=NEWLINE) as output_handle:
            fastq_iterator = (line[:-1] for line in input_handle)
            for record in zip_longest(*[fastq_iterator] * 4):
                title, sequence, quality = record[0], record[1], record[3]

                # Step 1: Filter by *poly-X*.
                is_filtered_out_by_poly_x = _filter_out_by_poly_x_trie_improved(sequence, poly_patterns)
                if is_filtered_out_by_poly_x:
                    num_filtered_out_by_poly_x += 1
                    continue

                # Step 2: Filter by *adapters*.
                is_filtered_out_by_adapters = _filter_out_by_adapters_trie_improved(sequence, adapters)
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
    # _worker_seq_naive_bp(input_fastq, output_fastq, output_stat, all_polyx_patterns, adapters)  # 40 s
    # _worker_seq_naive_pgzip_counter(input_fastq, output_fastq, output_stat, all_polyx_patterns, adapters)  # 11 s
    _worker_seq_naive_pgzip_zip(input_fastq, output_fastq, output_stat, all_polyx_patterns, adapters)  # 11 s


@time_it
def main_logic_seq_trie(input_fastq: Path, input_adapter: Path, output_fastq: Path, output_stat: Path) -> None:
    """High-level implementation of the main filtering logic. Separate Tries. Sequential."""
    all_polyx_patterns = _generate_all_polyx_patterns()
    adapters = _read_adapters(input_adapter, use_set=True)
    polyx_patterns_trie = build_trie(all_polyx_patterns)
    adapters_trie = build_trie(adapters)
    # print(len(polyx_patterns_trie), polyx_patterns_trie)  # 1477  {0: {'C': 1, 'T': 16, 'A': 63, 'G': 78}, ...
    # print(len(adapters_trie), adapters_trie)  # 57 {0: {'A': 1, 'C': 31}, 1: {'T': 2}, 2: {'A': 3}, 3: {'A': 4}, ...
    # _worker_seq_trie_bp(input_fastq, output_fastq, output_stat, polyx_patterns_trie, adapters_trie)  # 80 s
    # _worker_seq_trie_bp_generator(input_fastq, output_fastq, output_stat, polyx_patterns_trie, adapters_trie)  # 80 s
    # _worker_seq_trie_pgzip_counter(input_fastq, output_fastq, output_stat, polyx_patterns_trie, adapters_trie)  # 35 s
    _worker_seq_trie_pgzip_zip(input_fastq, output_fastq, output_stat, polyx_patterns_trie, adapters_trie)  # 35 s


@time_it
def main_logic_seq_trie_improved(input_fastq: Path, input_adapter: Path, output_fastq: Path, output_stat: Path) -> None:
    """High-level implementation of the main filtering logic. Separate Tries. Sequential."""
    all_polyx_patterns = _generate_all_polyx_patterns()
    adapters = _read_adapters(input_adapter, use_set=True)
    polyx_patterns_trie = build_trie_improved(all_polyx_patterns)
    adapters_trie = build_trie_improved(adapters)
    # _worker_improved_seq_trie_pgzip_counter(
    #     input_fastq, output_fastq, output_stat, polyx_patterns_trie, adapters_trie
    # )  # 35 s
    _worker_improved_seq_trie_pgzip_zip(
        input_fastq, output_fastq, output_stat, polyx_patterns_trie, adapters_trie
    )  # 35 s


@time_it
def main_logic_seq_trie_combined(input_fastq: Path, input_adapter: Path, output_fastq: Path, output_stat: Path) -> None:
    """High-level implementation of the main filtering logic. A combined Trie. Sequential."""
    combined_trie = _build_combined_trie(input_adapter)
    _worker_seq_trie_combined_bp(input_fastq, output_fastq, output_stat, combined_trie)  # 55 s


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
    main_logic_seq_naive(INPUT_FASTQ, INPUT_ADAPTER, OUTPUT_FASTQ_GZ, OUTPUT_STATISTICS)
    _validate_filtering()
    _validate_pgzip_decompresses_output_file(OUTPUT_FASTQ_GZ)
