"""
File:       src/__main__.py
Author:     Ivan Lazarević
Brief:      The high-level main business logic file.
"""
# Standard library imports
import os
import sys
from pathlib import Path

# Third party library imports

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import INPUT_FASTQ, INPUT_ADAPTER
from src.config import OUTPUT_FASTQ_GZ, OUTPUT_STATISTICS
from src.data import generate_all_polyx_patterns, read_adapters
from src.filters import FilterByPolyXTrie
from src.utils import time_it
from src.trie import build_trie


@time_it
def main_logic(input_fastq: Path, input_adapter: Path, output_fastq: Path, output_stat: Path) -> None:
    """High-level implementation of the main filtering logic."""
    all_polyx_patterns = generate_all_polyx_patterns()
    adapters = read_adapters(input_adapter, use_set=True)
    polyx_patterns_trie = build_trie(all_polyx_patterns)
    adapters_trie = build_trie(adapters)
    filter_ = FilterByPolyXTrie()
    filter_.worker(input_fastq, output_fastq, output_stat, polyx_patterns_trie, adapters_trie)


def main() -> None:
    """Local main"""
    main_logic(INPUT_FASTQ, INPUT_ADAPTER, OUTPUT_FASTQ_GZ, OUTPUT_STATISTICS)
