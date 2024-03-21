"""
File:       src/config.py
Author:     Ivan Lazarević
Brief:      Program configuration and constants.
"""
import platform
from pathlib import Path
from typing import Dict, Final

DEBUG: Final[bool] = True

USE_MODIN: Final[bool] = True

NUM_CPUS: Final[int] = 4
MODIN_CPUS: Final[str] = str(NUM_CPUS)

ALPHABET: Final[str] = "ACGT"
NEWLINE: Final[str] = "\n"  # Data files' line ending character(s). FASTQ uses "\n".
POLY_LEN: Final[int] = 15
ADAPTER_LEN: Final[int] = 32
POLY_END: Final[str] = "$"
ADAPTER_END: Final[str] = "#"
TRIE_IMPROVED_END: Final[str] = "@"

OPEN_PARAMS: Final[Dict[str, str]] = {
    "encoding": "ascii",
    "errors": "strict",
    "newline": NEWLINE,
}

# Use "\r" as separator on Windows, and "\n" elsewhere.
PANDAS_SEPARATOR = "\r" if platform.system() == "Windows" else "\n"
PANDAS_COLUMNS = ["read_id", "seq", "plus", "qual"]


# Data files
_INPUT_DATA_DIR = Path("input_data")
INPUT_ADAPTER = Path(_INPUT_DATA_DIR / "adapter.list")
INPUT_FASTQ = Path(_INPUT_DATA_DIR / "input.fq.gz")

_OUTPUT_DATA_DIR = Path("output_data")
OUTPUT_FASTQ_GZ = Path(_OUTPUT_DATA_DIR / "out.fq.gz")
OUTPUT_STATISTICS = Path(_OUTPUT_DATA_DIR / "out.stat.txt")
OUTPUT_TEXT_FILE = Path(_OUTPUT_DATA_DIR / "out.fq")

OUTPUT_FASTQ_SMALL_GZ = Path(_OUTPUT_DATA_DIR / "out_small.fq.gz")
OUTPUT_STATISTICS_SMALL = Path(_OUTPUT_DATA_DIR / "out_small.stat.txt")


# Test files
_TEST_INPUT_DATA_DIR = Path("tests/data/inp")
TEST_INP_ADAPTER = Path(_TEST_INPUT_DATA_DIR / "adapter.list")

TEST_INP_FASTQ = INPUT_FASTQ
TEST_INP_FASTQ_SMALL_ORIGINAL = Path(_TEST_INPUT_DATA_DIR / "input_small.fq")
TEST_INP_FASTQ_SMALL_GZ = Path(_TEST_INPUT_DATA_DIR / "input_small.fq.gz")

_TEST_OUTPUT_DATA_DIR = Path("tests/data/out")

TEST_OUT_FASTQ_GZ_REFERENCE = Path(_TEST_OUTPUT_DATA_DIR / "ref_out.fq.gz")
TEST_OUT_FASTQ_SIZE_REF = 16463630
TEST_OUT_STAT_REFERENCE = Path(_TEST_OUTPUT_DATA_DIR / "ref_out.stat.txt")

TEST_OUT_FASTQ_SMALL_REFERENCE_ORIGINAL = Path(_TEST_OUTPUT_DATA_DIR / "ref_out_small.fq")
TEST_OUT_FASTQ_SMALL_GZ_REFERENCE = Path(_TEST_OUTPUT_DATA_DIR / "ref_out_small.fq.gz")
TEST_OUT_FASTQ_SMALL_SIZE_REF = 940
TEST_OUT_SMALL_STAT_REFERENCE = Path(_TEST_OUTPUT_DATA_DIR / "ref_out_small.stat.txt")
