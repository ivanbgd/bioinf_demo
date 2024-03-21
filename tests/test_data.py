"""
File:       tests/test_data.py
Author:     Ivan Lazarević
Brief:      Unit tests for data-reading or data-creating functions.
"""
# Standard library imports
import os
import sys
import unittest

# Third party library imports

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import ALPHABET, POLY_LEN
from src.data import generate_all_polyx_patterns


class TestData(unittest.TestCase):
    """Class for automated testing of data fetching"""

    def test_generate_all_polyx_patterns(self):
        test_patterns = ["CAAAAAAAAAAAAAA", "GAAAAAAAAAAAAAA",
                         "TAAAAAAAAAAAAAA", "CTTTTTTTTTTTTTT",
                         "AAAAAAAAAAAAAAG", "TTTTTTTTTTTTTTA"]
        expected_len = len(ALPHABET) * (1 + POLY_LEN * (len(ALPHABET) - 1))
        all_polys = generate_all_polyx_patterns()
        sorted_patterns = sorted(all_polys)
        for pattern in sorted_patterns:
            self.assertEqual(POLY_LEN, len(pattern))
        for pattern in test_patterns:
            self.assertIn(pattern, all_polys)
        self.assertEqual(expected_len, len(all_polys))


if __name__ == "__main__":
    unittest.main(argv=[""], verbosity=2, exit=False)
