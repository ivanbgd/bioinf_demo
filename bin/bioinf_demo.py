#!/usr/bin/env python3
"""
File:       bin/bioinf_demo.py
Author:     Ivan Lazarević
Brief:      The demo's main "executable" file.
"""
# Standard library imports
import os
import sys

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
import src.__main__  # noqa
from src import main  # noqa
from src.config import USE_MODIN  # noqa
from src.initialize import initialize_dask_and_modin  # noqa


if __name__ == "__main__":
    """Program entry point"""
    if USE_MODIN:
        initialize_dask_and_modin()

    main()
