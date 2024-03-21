"""
File:       tests/test_fastq_reader.py
Author:     Ivan Lazarević
Brief:      Automated unit tests for the FASTQ Reader.
"""
# Standard library imports
import gzip
import os
import sys
import unittest

# Third party library imports

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import NEWLINE, TEST_OUT_FASTQ_GZ_REFERENCE
from src.fastq_reader import read_fastq


class TestFastqReader(unittest.TestCase):
    """Class for automated testing of the FASTQ Reader"""

    def test_read_fastq(self):
        """Validate that this solution properly reads the reference "out.fq.gz" file"""
        solution_contents = "".join(read_fastq(TEST_OUT_FASTQ_GZ_REFERENCE))
        with gzip.open(TEST_OUT_FASTQ_GZ_REFERENCE, "rt", encoding="ascii", errors="strict", newline=NEWLINE) \
                as reference_handle:
            reference_contents = reference_handle.read()
        self.assertEqual(reference_contents, solution_contents)


if __name__ == "__main__":
    unittest.main(argv=[""], verbosity=2, exit=False)
