"""
File:       src/initialize.py
Author:     Ivan Lazarević
Brief:      The demo's main "executable" file.
"""
# Standard library imports
import os

# Local modules imports
from src.config import MODIN_CPUS


def initialize_dask():
    from distributed import Client
    _ = Client()


def _initialize_modin():
    os.environ["MODIN_CPUS"] = MODIN_CPUS
    import modin
    modin.config.Engine.put("Dask")  # noqa
    print(f"Modin num partitions = {modin.config.NPartitions.get()}")  # noqa


def initialize_dask_and_modin():
    initialize_dask()
    _initialize_modin()
