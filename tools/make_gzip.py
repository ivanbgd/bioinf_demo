#!/usr/bin/env python3
"""
File:       tools/make_gzip.py
Author:     Ivan Lazarević
Brief:      Script for making test *gzip* files.

Details:
            It is enough to run the script once from CLI.
            Alternatively, the `compress_file function can be called programmatically.
"""
# Standard library imports
import gzip
import os
import sys
from pathlib import Path

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import OPEN_PARAMS
from src.config import TEST_INP_FASTQ_SMALL_ORIGINAL, TEST_INP_FASTQ_SMALL_GZ
from src.config import TEST_OUT_FASTQ_SMALL_REFERENCE_ORIGINAL, TEST_OUT_FASTQ_SMALL_GZ_REFERENCE


def compress_file(source: Path, dst: Path) -> None:
    with open(source, "rt", **OPEN_PARAMS) as source_handle:
        with gzip.open(dst, "wt", **OPEN_PARAMS) as dst_handle:
            contents = source_handle.read()
            dst_handle.write(contents)


if __name__ == "__main__":
    compress_file(TEST_INP_FASTQ_SMALL_ORIGINAL, TEST_INP_FASTQ_SMALL_GZ)
    compress_file(TEST_OUT_FASTQ_SMALL_REFERENCE_ORIGINAL, TEST_OUT_FASTQ_SMALL_GZ_REFERENCE)
