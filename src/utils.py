"""
File:       src/utils.py
Author:     Ivan Lazarević
Brief:      Utility functions.
"""
# Standard library imports
import datetime
import sys
import timeit
from functools import wraps
from typing import Callable, NoReturn

# Local modules imports
from src.config import DEBUG


def exit_program(msg: str) -> NoReturn:
    """Print message and exit program"""
    print(f"\n{msg}\nExiting program.")
    sys.exit(1)


def time_it(function: Callable) -> Callable:
    if not DEBUG:
        return function

    @wraps(function)
    def inner(*args, **kw):
        start = timeit.default_timer()
        result = function(*args, **kw)
        end = timeit.default_timer()
        diff = end - start
        print(f"\n\t*** TIMING: {function} took {datetime.timedelta(seconds=diff)} or {diff:.3f} s to complete.\n")
        return result
    return inner
