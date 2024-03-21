#!/usr/bin/env python
"""
File:       exp/exp_keep_ids.py
Author:     Ivan Lazarević
Brief:      Script for experimenting with reading, parsing and writing gzipped FASTQ files, filtering by sequence IDs.
"""
# Standard library imports
import gzip
import os
import sys

# Third party library imports
import pgzip
from Bio import Seq, SeqIO, bgzf
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import NEWLINE, INPUT_FASTQ, OUTPUT_FASTQ_GZ, OUTPUT_TEXT_FILE
from src.utils import time_it


keep_ids_lst = ["xx", "V350044321L1C001R0040217861/1", "V350044321L1C001R0040000003/1", "V350044321L1C001R0040000011/1"]

# https://www.biostars.org/p/10353/#10364
# Read your read names into list1 and change list to set.
# Set is hashable, so checking for presence of an element is much faster than in a list.
keep_ids = set(keep_ids_lst)


@time_it  # 2.65 s
def parse_keep_ids_1():
    """Print to *stdout* only, in a lazy fashion. Uses Biopython."""
    with gzip.open(INPUT_FASTQ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        filtered = (record for record in SeqIO.parse(handle=input_handle, format="fastq") if record.id in keep_ids)
        # for record in filtered:
        #     print(record)
        num_rec_stdout = SeqIO.write(filtered, sys.stdout, "fastq")
        print("\n", num_rec_stdout)


@time_it  # 2.45 s
def parse_keep_ids_2():
    """Print to *stdout* and a file, but eagerly. Uses Biopython."""
    with gzip.open(INPUT_FASTQ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with bgzf.BgzfWriter(OUTPUT_FASTQ_GZ, "wb") as output_handle:
            filtered = [record for record in SeqIO.parse(handle=input_handle, format="fastq") if record.id in keep_ids]
            num_rec_stdout = SeqIO.write(filtered, sys.stdout, "fastq")
            num_rec_file = SeqIO.write(sequences=filtered, handle=output_handle, format="fastq")
            print("\n", num_rec_stdout, num_rec_file)


@time_it  # 2.45 s
def parse_keep_ids_3():
    """Print to *stdout* and a file, lazily. Uses Biopython."""
    with gzip.open(INPUT_FASTQ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with bgzf.BgzfWriter(OUTPUT_FASTQ_GZ, "wb") as output_handle:
            for record in SeqIO.parse(handle=input_handle, format="fastq"):
                if record.id in keep_ids:
                    num_rec_stdout = SeqIO.write(record, sys.stdout, "fastq")
                    num_rec_file = SeqIO.write(sequences=record, handle=output_handle, format="fastq")
                    # print("**", num_rec_stdout, num_rec_file)
                    # print(record)
                    # print(str(record))
                    # print(repr(record))
                    # print(record.letter_annotations)
                    # print(list(record.letter_annotations))
                    # print(record.letter_annotations["phred_quality"])


@time_it  # 0.65 s
def parse_keep_ids_4():
    """Print to *stdout* and a file, lazily. Uses Biopython. Multiple times faster than *parse3*."""
    # http://maq.sourceforge.net/fastq.shtml
    with gzip.open(INPUT_FASTQ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with bgzf.BgzfWriter(OUTPUT_FASTQ_GZ, "wb") as output_handle:
            for title, sequence, quality in FastqGeneralIterator(input_handle):
                if title.split(None, 1)[0] in keep_ids:
                    print(f"@{title}\n{sequence}\n+\n{quality}")
                    # print(f"@{title}{NEWLINE}{sequence}{NEWLINE}+{NEWLINE}{quality}")
                    record = SeqIO.SeqRecord(
                        id=title,
                        seq=Seq.Seq(sequence),
                        # name=title,
                        description=title,
                        # dbxrefs=[],
                        # features=[0],
                        letter_annotations={"phred_quality": [ord(q) - 33 for q in quality]},
                    )
                    num_rec_file = SeqIO.write(sequences=record, handle=output_handle, format="fastq")
                    # print("**", num_rec_file)


@time_it  # 0.45 s
def parse_keep_ids_5():
    """Print to *stdout* and a file, lazily. Uses Biopython for writing only. A little faster than *parse4*."""
    _keep_ids = set("@" + id_ + NEWLINE for id_ in keep_ids_lst)
    # print(_keep_ids)
    len_newline = len(NEWLINE)
    # print(len_newline)
    flag = False
    with gzip.open(INPUT_FASTQ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with bgzf.BgzfWriter(OUTPUT_FASTQ_GZ, "wt") as output_handle:
            for line_no, line in enumerate(input_handle, 1):
                if line_no % 4 == 1:
                    flag = line in _keep_ids
                    title = line.strip()[1:]
                    # title = line[1:-len_newline]
                # print("!!!", line_no, flag, line)
                if flag:
                    print(line, end="\b")
                    if line_no % 4 == 2:
                        sequence = line.strip()
                        # sequence = line[:-len_newline]
                        # print("@@@", sequence)
                    if line_no % 4 == 0:
                        quality = line.strip()
                        # quality = line[:-len_newline]
                        # print("???", title, sequence, quality)
                        record = SeqIO.SeqRecord(
                            id=title,
                            seq=Seq.Seq(sequence),
                            description=title,
                            letter_annotations={"phred_quality": [ord(q) - 33 for q in quality]},
                        )
                        num_rec_file = SeqIO.write(sequences=record, handle=output_handle, format="fastq")
                        # print("**", num_rec_file, "\n")
                # if line_no == 3 * 4:
                #     return


@time_it  # 0.4 s
def parse_keep_ids_6_text():
    """Print to *stdout* and a file, lazily. GZIP. No Biopython at all. Text mode. A little faster than *parse4*."""
    _keep_ids = set("@" + id_ + NEWLINE for id_ in keep_ids_lst)
    flag = False
    with gzip.open(INPUT_FASTQ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with open(OUTPUT_TEXT_FILE, "wt", encoding="ascii", errors="strict", newline=NEWLINE) as output_handle:
            for line_no, line in enumerate(input_handle, 1):
                if line_no % 4 == 1:
                    flag = line in _keep_ids
                    title = line[1:-1]
                if flag:
                    print(line, end="\b")
                    if line_no % 4 == 2:
                        sequence = line[:-1]
                    if line_no % 4 == 0:
                        quality = line[:-1]
                        output_handle.write(f"@{title}{NEWLINE}{sequence}{NEWLINE}+{NEWLINE}{quality}{NEWLINE}")

    with open(OUTPUT_TEXT_FILE, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        contents_to_compress = input_handle.read()
    with gzip.open(OUTPUT_FASTQ_GZ, "wt", encoding="ascii", errors="strict", newline=NEWLINE) as output_handle:
        output_handle.write(contents_to_compress)


@time_it  # 0.4 s
def parse_keep_ids_6_binary():
    """Print to *stdout* and a file, lazily. GZIP. No Biopython at all. Binary mode. A little faster than *parse4*."""
    _keep_ids = set("@" + id_ + NEWLINE for id_ in keep_ids_lst)
    flag = False
    # Both these files must be opened in text mode.
    with gzip.open(INPUT_FASTQ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with open(OUTPUT_TEXT_FILE, "wt", encoding="ascii", errors="strict", newline=NEWLINE) as output_handle:
            for line_no, line in enumerate(input_handle, 1):
                if line_no % 4 == 1:
                    flag = line in _keep_ids
                    title = line[1:-1]
                if flag:
                    print(line, end="\b")
                    if line_no % 4 == 2:
                        sequence = line[:-1]
                    if line_no % 4 == 0:
                        quality = line[:-1]
                        output_handle.write(f"@{title}{NEWLINE}{sequence}{NEWLINE}+{NEWLINE}{quality}{NEWLINE}")

    with open(OUTPUT_TEXT_FILE, "rb") as input_handle:
        contents_to_compress = input_handle.read()
    with gzip.open(OUTPUT_FASTQ_GZ, "wb") as output_handle:
        output_handle.write(contents_to_compress)


# Can use this.
@time_it  # 0.4 s
def parse_keep_ids_6():
    """Print to *stdout* and a file, lazily. GZIP. No Biopython at all. Text mode. A little faster than *parse4*."""
    _keep_ids = set("@" + id_ + NEWLINE for id_ in keep_ids_lst)
    flag = False
    with gzip.open(INPUT_FASTQ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with gzip.open(OUTPUT_FASTQ_GZ, "wt", encoding="ascii", errors="strict", newline=NEWLINE) as output_handle:
            for line_no, line in enumerate(input_handle, 1):
                if line_no % 4 == 1:
                    flag = line in _keep_ids
                    title = line[:-1]
                if flag:
                    print(line, end="\b")
                    if line_no % 4 == 2:
                        sequence = line[:-1]
                    if line_no % 4 == 0:
                        quality = line[:-1]
                        output_handle.write(f"{title}{NEWLINE}{sequence}{NEWLINE}+{NEWLINE}{quality}{NEWLINE}")


# Can use this.
@time_it  # 0.4 s
def parse_keep_ids_7():
    """Print to *stdout* and a file, lazily. Uses *pgzip* instead of *gzip*. A little faster than *parse4*."""
    _keep_ids = set("@" + id_ + NEWLINE for id_ in keep_ids_lst)
    flag = False
    with pgzip.open(INPUT_FASTQ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with pgzip.open(OUTPUT_FASTQ_GZ, "wt", encoding="ascii", errors="strict", newline=NEWLINE) as output_handle:
            for line_no, line in enumerate(input_handle, 1):
                if line_no % 4 == 1:
                    flag = line in _keep_ids
                    title = line[:-1]
                if flag:
                    print(line, end="\b")
                    if line_no % 4 == 2:
                        sequence = line[:-1]
                    if line_no % 4 == 0:
                        quality = line[:-1]
                        output_handle.write(f"{title}{NEWLINE}{sequence}{NEWLINE}+{NEWLINE}{quality}{NEWLINE}")


@time_it
def parse_keep_ids_8():
    """Print to *stdout* and a file, lazily. Using Biopython's "index". A little faster than *parse4* ???."""
    # TODO: Try Biopython with "index"; fix docstring w.r.t. speed.
    raise NotImplementedError


@time_it  # 0.55 s
def parse_keep_ids_9():
    """Print to *stdout* and a file, lazily. No Biopython at all. A little faster than *parse4*."""
    _keep_ids = set("@" + id_ + NEWLINE for id_ in keep_ids_lst)
    flag = False
    eof = False
    with gzip.open(INPUT_FASTQ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as input_handle:
        with gzip.open(OUTPUT_FASTQ_GZ, "wt", encoding="ascii", errors="strict", newline=NEWLINE) as output_handle:
            line_no = 1
            while not eof:
                line = input_handle.readline()
                eof = line == ""
                if line_no % 4 == 1:
                    flag = line in _keep_ids
                    title = line[:-1]
                if flag:
                    print(line, end="\b")
                    if line_no % 4 == 2:
                        sequence = line[:-1]
                    if line_no % 4 == 0:
                        quality = line[:-1]
                        output_handle.write(f"{title}{NEWLINE}{sequence}{NEWLINE}+{NEWLINE}{quality}{NEWLINE}")
                    line_no += 1
                else:
                    line_no += 4


def visually_validate_parsing_gzip():
    """*gzip* wants to decompress gzip files written by Biopython."""
    print("\n\n VALIDATION (gzip) \n\n")
    with gzip.open(OUTPUT_FASTQ_GZ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as handle:
        contents = handle.read()
    print(contents)


def visually_validate_parsing_pgzip():
    """*pgzip* doesn't want to decompress gzip files written by Biopython."""
    print("\n\n VALIDATION (pgzip) \n\n")
    with pgzip.open(OUTPUT_FASTQ_GZ, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as handle:
        contents = handle.read()
    print(contents)


if __name__ == "__main__":
    parse_keep_ids_6()
    visually_validate_parsing_gzip()
    visually_validate_parsing_pgzip()
