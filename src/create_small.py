"""
File:       src/create_small.py
Author:     Ivan Lazarević
Brief:      Helper for creating small output files for testing.
"""
# Standard library imports
from pathlib import Path

# Local modules imports
from src.data import generate_all_polyx_patterns, read_adapters
from src.filters import FilterByPolyXNaive


def create_small_output_files(input_fastq: Path, input_adapter: Path, output_fastq: Path, output_stat: Path) -> None:
    """Use this to create `OUTPUT_FASTQ_SMALL` and `OUTPUT_STATISTICS_SMALL` for testing."""
    all_polyx_patterns = generate_all_polyx_patterns()
    adapters = read_adapters(input_adapter, use_set=True)
    filter_ = FilterByPolyXNaive()
    filter_.worker(input_fastq, output_fastq, output_stat, all_polyx_patterns, adapters)
