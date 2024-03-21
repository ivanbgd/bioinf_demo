"""
File:       tests/__init__.py
Author:     Ivan Lazarević
Brief:      The "tests" directory init file.
"""
# Standard library imports
import os
import sys

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import USE_MODIN  # noqa
from src.initialize import initialize_dask  # noqa


if USE_MODIN:
    initialize_dask()
