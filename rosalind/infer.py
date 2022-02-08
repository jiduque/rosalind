from typing import List, Tuple

from rosalind._data_structures import BinarySearchTree
from rosalind.constants import PROTEIN_MASS_MAP
from rosalind.compare import spectral_convolution


def protein_from_mass_spectrum(mass_data: List[float]) -> str:
    inverted_key_values = list(map(lambda x: (x[1], x[0]), PROTEIN_MASS_MAP.items()))
    search_tree = BinarySearchTree(inverted_key_values)
    mass_diff = list(map(lambda x, y: y - x, mass_data[:-1], mass_data[1:]))
    return "".join(map(lambda x: search_tree.search_closest(x).val, mass_diff))


def spectrum_to_protein(mass_data: List[float], proteins: List[str]) -> Tuple[int, str]:
    def complete_spectrum(protein):
        n = len(protein)
        masses = list(map(lambda x: PROTEIN_MASS_MAP[x], protein))
        output = []
        for i in range(n):
            output.append(sum(masses[i:]))
            output.append(sum(masses[:i]))

        return output

    return max(
        [(spectral_convolution(mass_data, complete_spectrum(protein))[0], protein) for protein in proteins],
        key=lambda x: x[0]
    )
