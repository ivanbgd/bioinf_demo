"""
File:       tests/test_main_logic.py
Author:     Ivan Lazarević
Brief:      Automated unit tests for the main logic.
"""
# Standard library imports
import gzip
import os
import sys
import unittest
from pathlib import Path

# Third party library imports
import pgzip

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import DEBUG, OPEN_PARAMS
from src.config import OUTPUT_FASTQ_GZ, OUTPUT_STATISTICS
from src.config import OUTPUT_FASTQ_SMALL_GZ, OUTPUT_STATISTICS_SMALL
from src.config import TEST_INP_ADAPTER
from src.config import TEST_INP_FASTQ_SMALL_ORIGINAL, TEST_INP_FASTQ_SMALL_GZ
from src.config import TEST_OUT_FASTQ_GZ_REFERENCE, TEST_OUT_FASTQ_SIZE_REF, TEST_OUT_STAT_REFERENCE
from src.config import TEST_OUT_FASTQ_SMALL_REFERENCE_ORIGINAL, TEST_OUT_FASTQ_SMALL_GZ_REFERENCE
from src.config import TEST_OUT_FASTQ_SMALL_SIZE_REF, TEST_OUT_SMALL_STAT_REFERENCE
from src.create_small import create_small_output_files
from tools.make_gzip import compress_file


class TestMain(unittest.TestCase):
    """Class for automated testing of the main business logic"""

    @classmethod
    def setUpClass(cls) -> None:
        if not TEST_INP_FASTQ_SMALL_GZ.exists():
            print(f"'{TEST_INP_FASTQ_SMALL_GZ}' doesn't exist. Creating it now...")
            compress_file(TEST_INP_FASTQ_SMALL_ORIGINAL, TEST_INP_FASTQ_SMALL_GZ)
        if not TEST_OUT_FASTQ_SMALL_GZ_REFERENCE.exists():
            print(f"'{TEST_OUT_FASTQ_SMALL_GZ_REFERENCE}' doesn't exist. Creating it now...")
            compress_file(TEST_OUT_FASTQ_SMALL_REFERENCE_ORIGINAL, TEST_OUT_FASTQ_SMALL_GZ_REFERENCE)

        OUTPUT_FASTQ_SMALL_GZ.unlink(missing_ok=True)
        OUTPUT_STATISTICS_SMALL.unlink(missing_ok=True)
        create_small_output_files(
            TEST_INP_FASTQ_SMALL_GZ, TEST_INP_ADAPTER, OUTPUT_FASTQ_SMALL_GZ, OUTPUT_STATISTICS_SMALL
        )

    @classmethod
    def tearDownClass(cls) -> None:
        if not DEBUG:
            print(f"Deleting '{OUTPUT_FASTQ_SMALL_GZ}' and '{OUTPUT_STATISTICS_SMALL}'...")
            OUTPUT_FASTQ_SMALL_GZ.unlink(missing_ok=True)
            OUTPUT_STATISTICS_SMALL.unlink(missing_ok=True)

    def _validate_output_fastq(self, output_fastq: Path, test_out_fastq_reference: Path) -> None:
        with gzip.open(output_fastq, "rt", **OPEN_PARAMS) as solution_handle:
            solution_contents = solution_handle.read()
        with gzip.open(test_out_fastq_reference, "rt", **OPEN_PARAMS) as reference_handle:
            reference_contents = reference_handle.read()
        self.assertEqual(len(reference_contents), len(solution_contents))
        self.assertEqual(reference_contents, solution_contents)

    def _validate_output_stat(self, output_fastq: Path, test_out_fastq_reference: Path) -> None:
        with open(output_fastq, "rt", **OPEN_PARAMS) as solution_handle:
            solution_contents = solution_handle.read()
        with open(test_out_fastq_reference, "rt", **OPEN_PARAMS) as reference_handle:
            reference_contents = reference_handle.read()
        self.assertEqual(len(reference_contents), len(solution_contents))
        self.assertEqual(reference_contents, solution_contents)

    def test_out_fq_gz(self):
        """Verify that reference and this solution's "out.fq.gz" files are the same"""
        self._validate_output_fastq(OUTPUT_FASTQ_GZ, TEST_OUT_FASTQ_GZ_REFERENCE)

    def test_out_stat_txt(self):
        """Verify that reference and this solution's "out.stat.txt" files are the same"""
        self._validate_output_stat(OUTPUT_STATISTICS, TEST_OUT_STAT_REFERENCE)

    def test_gzip_decompresses_out_fq_gz(self) -> None:
        """Verify by *gzip* that the output *gzip* file was written correctly"""
        with gzip.open(OUTPUT_FASTQ_GZ, "rb") as solution_handle:
            self.assertEqual(TEST_OUT_FASTQ_SIZE_REF, len(solution_handle.read()))

    def test_pgzip_decompresses_out_fq_gz(self) -> None:
        """Biopython doesn't write *gzip* files properly according to *pgzip*, so validate by *pgzip*"""
        with pgzip.open(OUTPUT_FASTQ_GZ, "rb", thread=None) as solution_handle:
            self.assertEqual(TEST_OUT_FASTQ_SIZE_REF, len(solution_handle.read()))

    def test_out_small_fq_gz(self):
        """Verify that reference and this solution's "out_small.fq.gz" files are the same"""
        self._validate_output_fastq(OUTPUT_FASTQ_SMALL_GZ, TEST_OUT_FASTQ_SMALL_GZ_REFERENCE)

    def test_out_small_stat_txt(self):
        """Verify that reference and this solution's "out_small.stat.txt" files are the same"""
        self._validate_output_stat(OUTPUT_STATISTICS_SMALL, TEST_OUT_SMALL_STAT_REFERENCE)

    def test_gzip_decompresses_out_small_fq_gz(self) -> None:
        """Verify by *gzip* that the output small *gzip* file was written correctly"""
        with gzip.open(OUTPUT_FASTQ_SMALL_GZ, "rt", **OPEN_PARAMS) as solution_handle:
            self.assertEqual(TEST_OUT_FASTQ_SMALL_SIZE_REF, len(solution_handle.read()))

    def test_pgzip_decompresses_out_small_fq_gz(self) -> None:
        """Biopython doesn't write *gzip* files properly according to *pgzip*, so validate by *pgzip*"""
        with pgzip.open(OUTPUT_FASTQ_SMALL_GZ, "rt", **OPEN_PARAMS, thread=None) \
                as solution_handle:
            self.assertEqual(TEST_OUT_FASTQ_SMALL_SIZE_REF, len(solution_handle.read()))


if __name__ == "__main__":
    unittest.main(argv=[""], verbosity=2, exit=False)
