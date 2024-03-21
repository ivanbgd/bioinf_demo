"""
File:       src/fastq_reader.py
Author:     Ivan Lazarević
Brief:      FASTQ input facilities.
"""
# Standard library imports
import gzip
import os
import sys
from abc import ABC, abstractmethod
from pathlib import Path
from typing import cast, List, NoReturn, Union

# Third party library imports
from Bio import SeqIO  # noqa
from Bio.SeqRecord import SeqRecord  # noqa

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import NEWLINE  # noqa
from src.utils import time_it  # noqa


class _FastqReader(ABC):
    """Abstract Base Class for FASTQ readers"""

    @abstractmethod
    def read(self, path: Path) -> str:
        """Read eagerly (fully) a FASTQ file in *gzip* format and return its contents as string"""
        pass

    @abstractmethod
    def read_list(self, path: Path) -> Union[List[str], List[SeqRecord]]:
        """
        Read lazily, element by element, a FASTQ file in *gzip* format
        and return its contents as list of `str` or `SeqRecord`
        """
        pass


class _FastqReaderSequential(_FastqReader):
    """A sequential implementation of a FASTQ reader which uses the *gzip* module"""

    @time_it
    def read(self, path: Path) -> str:
        with gzip.open(path, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as handle:
            contents = handle.read()
        return contents

    @time_it
    def read_list(self, path: Path) -> List[str]:
        with gzip.open(path, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as handle:
            contents = handle.readlines()
        return contents


class _FastqReaderPyfastx(_FastqReader):
    """A *pyfastx* implementation of a FASTQ reader"""

    @time_it
    def read(self, path: Path) -> NoReturn:
        raise NotImplementedError("*pyfastx* won't install on python 3.10")

    @time_it
    def read_list(self, path: Path) -> NoReturn:
        raise NotImplementedError("*pyfastx* won't install on python 3.10")


class _FastqReaderBiopython(_FastqReader):
    """A *biopython* implementation of a FASTQ reader"""

    @time_it
    def read(self, path: Path) -> str:
        raise NotImplementedError("Phred qualities missing")
        records = self.read_list(path)  # noqa
        contents = NEWLINE.join(f"{rec.id}\n{rec.seq}\n+\n" for rec in records)
        return contents

    @time_it
    def read_list(self, path: Path) -> List[SeqRecord]:
        with gzip.open(path, "rt", encoding="ascii", errors="strict", newline=NEWLINE) as handle:
            contents = [r for r in SeqIO.parse(handle, "fastq")]
        return contents


def compare_fastq_readers(input_path: Path) -> None:
    """Compare speed of various `_FastqReader` implementations reading the same FASTQ file"""
    readers: List[_FastqReader] = [
        cast(_FastqReader, _FastqReaderSequential()),
        cast(_FastqReader, _FastqReaderBiopython()),
    ]
    reader: _FastqReader
    for reader in readers:
        reader.read_list(input_path)


def read_fastq(input_path: Path) -> Union[str, List[str], List[SeqRecord]]:
    """Read a FASTQ file and return its contents by using the fastest `_FastqReader` implementation"""
    reader: _FastqReader = _FastqReaderSequential()
    contents = reader.read_list(input_path)
    return contents


# def is_file_plain_fastq(read_file: Path) -> bool:
#     pass
