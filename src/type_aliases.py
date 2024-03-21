"""
File:       src/type_aliases.py
Author:     Ivan Lazarević
Brief:      Type aliases for type annotations.
"""
from typing import Dict, List, Set, Union

# Type aliases
Trie = Dict[int, Dict]
AdaptersNaive = Union[List[str], Set[str]]

Adapters = Union[AdaptersNaive, Trie]
PolyPatterns = Union[Set[str], Trie]
