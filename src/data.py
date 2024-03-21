"""
File:       src/data.py
Author:     Ivan Lazarević
Brief:      Data-reading or data-creating functions.
"""
# Standard library imports
import os
import sys
from pathlib import Path
from typing import List, Set, Union

# Third party library imports

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import ALPHABET, NEWLINE, ADAPTER_LEN, POLY_LEN  # noqa


def generate_all_polyx_patterns() -> Set[str]:
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


def read_adapters(input_adapter: Path, *, use_set: bool = False) -> Union[List[str], Set[str]]:
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
