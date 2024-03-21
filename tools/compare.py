#!/usr/bin/env python3
"""
File:       tools/compare.py
Author:     Ivan Lazarević
Brief:      Script for comparing speed of various algorithms and data structures.
"""
# Standard library imports
import os
import sys

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import INPUT_FASTQ  # noqa
from src.fastq_reader import compare_fastq_readers  # noqa


def _compare_fastq_readers():
    print("\nComparing FASTQ readers...")
    compare_fastq_readers(INPUT_FASTQ)


def _compare_everything() -> None:
    _compare_fastq_readers()


if __name__ == "__main__":
    _compare_everything()
