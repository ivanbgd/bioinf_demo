"""
File:       src/filters.py
Author:     Ivan Lazarević
Brief:      Filters for PolyX and Adapters.
"""
# Standard library imports
import csv
from abc import ABC, abstractmethod
from pathlib import Path

# Third party library imports
import modin.pandas as pd

# Local modules imports
from src.config import NEWLINE
from src.config import PANDAS_SEPARATOR, PANDAS_COLUMNS
from src.trie import trie_matching
from src.type_aliases import Adapters, AdaptersNaive, PolyPatterns, Trie
from src.utils import time_it


class Filter(ABC):
    """ Abstract Base Class for Filters

        Can be used for filtering out records by *PolyX* or by *Adapters*.
        Can be used to compare different implementations of the filtering algorithm.
    """

    @abstractmethod
    def filter_out_by_poly_x(self, record: str, patterns: PolyPatterns) -> bool:
        """Searches for a pattern from *patterns* in the *record* and if it finds one, returns True, otherwise False."""
        pass

    @abstractmethod
    def filter_out_by_adapters(self, record: str, patterns: Adapters) -> bool:
        """Searches for a pattern from *patterns* in the *record* and if it finds one, returns True, otherwise False."""
        pass

    @time_it
    def worker(self,
               input_fastq: Path,
               output_fastq: Path,
               output_stat: Path,
               poly_patterns: PolyPatterns,
               adapters: Adapters) -> None:
        """
        Low-level implementation of the main filtering logic.
        Uses Pandas Series apply. This is recommended in Pandas.
        Uses Pandas for reading the input file, which means Modin can parallelize that part.
        Uses Pandas for writing the output file, which means Modin can parallelize that part.
        """
        input_frame_modin = pd.DataFrame(
            pd.read_csv(
                input_fastq, sep=PANDAS_SEPARATOR, header=None
            ).values.reshape(-1, 4), columns=PANDAS_COLUMNS
        )

        input_frame = input_frame_modin

        # Step 1: Filter by *poly-X*.
        are_filtered_out_by_poly_x = input_frame["seq"].apply(
            lambda seq: self.filter_out_by_poly_x(seq, poly_patterns))
        filtered_out_by_poly_x = are_filtered_out_by_poly_x[are_filtered_out_by_poly_x].index
        num_filtered_out_by_poly_x = are_filtered_out_by_poly_x.value_counts()[True]

        # We can work with one data frame only, *input_frame*, modifying it, but we'll create a new data frame.
        output_frame = input_frame.drop(filtered_out_by_poly_x)

        # Step 2: Filter by *adapters*.
        are_filtered_out_by_adapters = output_frame["seq"].apply(
            lambda seq: self.filter_out_by_adapters(seq, adapters))
        filtered_out_by_adapters = are_filtered_out_by_adapters[are_filtered_out_by_adapters].index
        num_filtered_out_by_adapters = are_filtered_out_by_adapters.value_counts()[True]

        output_frame.drop(filtered_out_by_adapters, inplace=True)

        # Step 3: Write the record to the output file if not filtered out.
        output_frame.to_csv(
            output_fastq, header=False, index=False, sep=NEWLINE,
            quoting=csv.QUOTE_NONE, line_terminator=NEWLINE, escapechar=NEWLINE
        )

        # Step 4: Store the number of records filtered out by *poly-X* and by *adapters*, respectively.
        print(num_filtered_out_by_poly_x, num_filtered_out_by_adapters)
        stats = f"filterByPolyX:\t{num_filtered_out_by_poly_x}{NEWLINE}" \
                f"filterByAdapter:\t{num_filtered_out_by_adapters}{NEWLINE}"
        with open(output_stat, "wt", newline=NEWLINE) as stat_handle:
            stat_handle.write(stats)


class FilterByPolyXNaive(Filter):
    """Naive implementations of filtering"""

    def filter_out_by_poly_x(self, record: str, patterns: AdaptersNaive) -> bool:
        for pattern in patterns:
            if pattern in record:
                return True
        return False

    def filter_out_by_adapters(self, record: str, patterns: AdaptersNaive) -> bool:
        for adapter in patterns:
            if adapter in record:
                return True
        return False


class FilterByPolyXTrie(Filter):
    """Trie implementations of filtering"""

    def filter_out_by_poly_x(self, record: str, patterns: Trie) -> bool:
        return trie_matching(record, patterns)

    def filter_out_by_adapters(self, record: str, patterns: Trie) -> bool:
        return trie_matching(record, patterns)
