from typing import List

from rosalind.constants import PROTEIN_MASS_MAP
from rosalind._data_structures import BinarySearchTree


def protein_from_mass_spectrum(mass_data: List[float]) -> str:
    inverted_key_values = list(map(lambda x: (x[1], x[0]), PROTEIN_MASS_MAP.items()))
    search_tree = BinarySearchTree(inverted_key_values)
    mass_diff = list(map(lambda x, y: y - x, mass_data[:-1], mass_data[1:]))
    return "".join(map(lambda x: search_tree.search_closest(x).val, mass_diff))
