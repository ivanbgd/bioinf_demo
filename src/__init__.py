"""
File:       src/__init__.py
Author:     Ivan Lazarević
Brief:      The "src" directory init file.
"""
# Standard library imports
import os
import sys

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
import src.__main__  # noqa


def main() -> None:
    """Second-level program entry point"""
    return src.__main__.main()
