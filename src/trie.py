"""
File:       src/trie.py
Author:     Ivan Lazarević
Brief:      The Trie data structure.
"""
# Standard library imports
import os
import sys
from typing import List, Optional

# Third party library imports

# Local modules imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), r".."))
from src.config import ADAPTER_END, POLY_END, TRIE_IMPROVED_END
from src.type_aliases import Trie


def build_trie(patterns) -> Trie:
    """ Return a Trie built from `patterns`.

        Returns a Trie built from `patterns` in the form of a dictionary of dictionaries,
        e.g. {0:{"A":1,"T":2},1:{"C":3}}, where the key of the external dictionary is the node ID (integer),
        and the internal dictionary contains all the trie edges outgoing from the corresponding node,
        and the keys are the letters on those edges, and the values are the node IDs to which these edges lead.

        External dictionary contains nodes, and the nodes' internal dictionaries contain their outgoing edges.
        Nodes are labeled by integers, uniquely.
        If a node is a leaf, i.e., it doesn't have an outgoing edge, it will contain an empty dictionary.
        So, we create a node by adding a new empty dictionary to it.
        In other words, whenever we create a node, we should create a dictionary and assign it to the node.
        If a node has an outgoing edge, we fill its dictionary. If not, its dictionary remains empty.

        This function executes in a millisecond for our trie.
    """
    trie = {}

    new_node_label = 0
    trie[new_node_label] = {}
    new_node_label += 1
    root = trie[0]
    for pattern in patterns:
        current_node = root
        for current_symbol in pattern:
            try:
                existing_node_label = current_node[current_symbol]
            except KeyError:
                current_node[current_symbol] = new_node_label
                trie[new_node_label] = {}
                current_node = trie[new_node_label]
                new_node_label += 1
            else:
                current_node = trie[existing_node_label]

    return trie


def _prefix_trie_matching(text: str, trie: Trie) -> Optional[str]:
    """The algorithm for matching a `trie` of patterns in `text`"""
    v = trie[0]
    i = 0
    symbol = text[i]
    pattern = []  # Pattern is spelled from root to v.
    while True:
        if not v:  # "v" is a leaf.
            return "".join(pattern)
        elif v.get(symbol, None):
            pattern.append(symbol)
            v = trie[v[symbol]]
            i += 1
            symbol = text[i] if i < len(text) else None
        else:
            return None


def trie_matching_positions(text: str, trie: Trie) -> List[int]:
    """A wrapper that takes `text` and `trie` of patterns and returns a list of all positions of patterns in `text`."""
    positions = []
    for i in range(len(text)):
        result = _prefix_trie_matching(text[i:], trie)
        if result is not None:
            positions.append(i)
    return positions


def trie_matching(text: str, trie: Trie) -> bool:
    """A wrapper that takes `text` and `trie` of patterns and returns whether a pattern is contained in `text`."""
    for i in range(len(text)):
        result = _prefix_trie_matching(text[i:], trie)
        if result is not None:
            return True
    return False


def _prefix_trie_matching_combined(text: str, trie: Trie) -> Optional[str]:
    """The algorithm for matching a `trie` of combined patterns in `text`."""
    v = trie[0]
    i = 0
    symbol = text[i]
    while True:
        if v.get(POLY_END, None):  # "v" is a leaf.
            return POLY_END
        elif v.get(ADAPTER_END, None):  # "v" is a leaf.
            return ADAPTER_END
        elif v.get(symbol, None):
            v = trie[v[symbol]]
            i += 1
            symbol = text[i] if i < len(text) else None
        else:
            return None


def trie_matching_combined(text: str, trie: Trie) -> Optional[str]:
    """
    A wrapper that takes `text` and `trie` of combined patterns and returns whether a pattern is contained in `text`.
    If it isn't contained, then it returns None.
    If it is contained, returns either `POLY_END` or `ADAPTER_END`.
    """
    for i in range(len(text)):
        result = _prefix_trie_matching_combined(text[i:], trie)
        if result is not None:
            return result
    return None


def build_trie_improved(patterns) -> Trie:
    """ Return a Trie built from `patterns` that supports patterns which are prefixes of other patterns.

        Returns a Trie built from `patterns` in the form of a dictionary of dictionaries,
        e.g. {0:{"A":1,"T":2},1:{"C":3}}, where the key of the external dictionary is the node ID (integer),
        and the internal dictionary contains all the trie edges outgoing from the corresponding node,
        and the keys are the letters on those edges, and the values are the node IDs to which these edges lead.

        External dictionary contains nodes, and the nodes' internal dictionaries contain their outgoing edges.
        Nodes are labeled by integers, uniquely.
        If a node is a leaf, i.e., it doesn't have an outgoing edge, it will contain an empty dictionary.
        So, we create a node by adding a new empty dictionary to it.
        In other words, whenever we create a node, we should create a dictionary and assign it to the node.
        If a node has an outgoing edge, we fill its dictionary. If not, its dictionary remains empty.

        Add label `TRIE_IMPROVED_END` to each node to mark it as a node that ends a pattern, if True.

        This function executes in a millisecond for our trie.
    """
    trie = {}

    new_node_label = 0
    trie[new_node_label] = {}
    trie[new_node_label][TRIE_IMPROVED_END] = False
    new_node_label += 1
    root = trie[0]
    for pattern in patterns:
        current_node = root
        for current_symbol in pattern:
            try:
                existing_node_label = current_node[current_symbol]
            except KeyError:
                current_node[current_symbol] = new_node_label
                trie[new_node_label] = {}
                trie[new_node_label][TRIE_IMPROVED_END] = False
                current_node = trie[new_node_label]
                new_node_label += 1
            else:
                current_node = trie[existing_node_label]
        current_node[TRIE_IMPROVED_END] = True

    return trie


def _prefix_trie_matching_improved(text: str, trie: Trie) -> Optional[str]:
    """
    The algorithm for matching a `trie` of patterns in `text`
    that supports cases of patterns being prefixes of other patterns.
    """
    v = trie[0]
    i = 0
    symbol = text[i]
    pattern = []  # Pattern is spelled from root to v.
    while True:
        if v[TRIE_IMPROVED_END]:  # "v" ends a pattern.
            return "".join(pattern)
        elif v.get(symbol, None):
            pattern.append(symbol)
            v = trie[v[symbol]]
            i += 1
            symbol = text[i] if i < len(text) else None
        else:
            return None


def trie_matching_improved(text: str, trie: Trie) -> bool:
    """
    A wrapper that takes `text` and `trie` of patterns and returns whether a pattern is contained in `text`.
    Supports case in which a pattern is a prefix of some other pattern.
    """
    for i in range(len(text)):
        result = _prefix_trie_matching_improved(text[i:], trie)
        if result is not None:
            return True
    return False
