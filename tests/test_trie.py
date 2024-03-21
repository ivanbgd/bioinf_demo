"""
File:       tests/test_trie.py
Author:     Ivan Lazarević
Brief:      Unit tests for the Trie data structure.
"""
# Standard library imports
import os
import sys
import unittest
from typing import Set

# Third party library imports

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.trie import build_trie, trie_matching_positions, trie_matching
from src.trie import build_trie_improved, trie_matching_improved


class TestTrie(unittest.TestCase):
    """Class for automated testing of the Trie data structure"""

    @staticmethod
    def _create_result(trie) -> Set[str]:
        result = set()
        for node in trie:
            for c in trie[node]:
                result.add(f"{node}->{trie[node][c]}:{c}")
        return result

    def test_build_trie_0(self):
        patterns = ["A", "T", "AC"]
        expected = {"0->1:A", "0->2:T", "1->3:C"}
        trie = build_trie(patterns)
        result = self._create_result(trie)
        self.assertSetEqual(expected, result)

    def test_build_trie_1(self):
        patterns = ["ATA"]
        expected = {"0->1:A", "2->3:A", "1->2:T"}
        trie = build_trie(patterns)
        result = self._create_result(trie)
        self.assertSetEqual(expected, result)

    def test_build_trie_2(self):
        patterns = ["AT", "AG", "AC"]
        expected = {"0->1:A", "1->4:C", "1->3:G", "1->2:T"}
        trie = build_trie(patterns)
        result = self._create_result(trie)
        self.assertSetEqual(expected, result)

    def test_build_trie_3(self):
        patterns = ["ATAGA", "ATC", "GAT"]
        expected = {"0->1:A", "1->2:T", "2->3:A", "3->4:G", "4->5:A", "2->6:C", "0->7:G", "7->8:A", "8->9:T"}
        trie = build_trie(patterns)
        result = self._create_result(trie)
        self.assertSetEqual(expected, result)

    def test_trie_matching_positions_1(self):
        text = "AAA"
        patterns = ["AA"]
        expected = [0, 1]
        trie = build_trie(patterns)
        result = trie_matching_positions(text, trie)
        self.assertListEqual(expected, result)

    def test_trie_matching_positions_2(self):
        text = "AA"
        patterns = ["T"]
        expected = []
        trie = build_trie(patterns)
        result = trie_matching_positions(text, trie)
        self.assertListEqual(expected, result)

    def test_trie_matching_positions_3(self):
        text = "AATCGGGTTCAATCGGGGT"
        patterns = ["ATCG", "GGGT"]
        expected = [1, 4, 11, 15]
        trie = build_trie(patterns)
        result = trie_matching_positions(text, trie)
        self.assertListEqual(expected, result)

    def test_trie_matching_1(self):
        text = "AAA"
        patterns = ["AA"]
        trie = build_trie(patterns)
        result = trie_matching(text, trie)
        self.assertTrue(result)

    def test_trie_matching_2(self):
        text = "AA"
        patterns = ["T"]
        trie = build_trie(patterns)
        result = trie_matching(text, trie)
        self.assertFalse(result)

    def test_trie_matching_3(self):
        text = "AATCGGGTTCAATCGGGGT"
        patterns = ["ATCG", "GGGT"]
        trie = build_trie(patterns)
        result = trie_matching(text, trie)
        self.assertTrue(result)

    def test_trie_matching_improved_1(self):
        text = "AAA"
        patterns = ["AA"]
        trie = build_trie_improved(patterns)
        result = trie_matching_improved(text, trie)
        self.assertTrue(result)

    def test_trie_matching_improved_2(self):
        text = "AA"
        patterns = ["T"]
        trie = build_trie_improved(patterns)
        result = trie_matching_improved(text, trie)
        self.assertFalse(result)

    def test_trie_matching_improved_3(self):
        text = "AATCGGGTTCAATCGGGGT"
        patterns = ["ATCG", "GGGT"]
        trie = build_trie_improved(patterns)
        result = trie_matching_improved(text, trie)
        self.assertTrue(result)

    def test_trie_matching_improved_4(self):
        text = "ACATA"
        patterns = ["AT", "A", "AG"]
        trie = build_trie_improved(patterns)
        result = trie_matching_improved(text, trie)
        self.assertTrue(result)


if __name__ == "__main__":
    unittest.main(argv=[""], verbosity=2, exit=False)
