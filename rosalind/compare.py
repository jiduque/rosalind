from typing import Tuple, List, Dict
from itertools import product
from collections import Counter

from suffix_tree import Tree


def transition_transversion_ratio(dna1: str, dna2: str) -> float:
    transition_map = {"A": "G", "G": "A", "C": "T", "T": "C"}
    transitions, transversions = 0, 0

    for x, y in zip(dna1, dna2):
        if transition_map[x] == y:
            transitions += 1
        elif x != y:
            transversions += 1

    return transitions / transversions if transversions else float("inf")


def hamming_distance(dna_1: str, dna_2: str) -> int:
    assert len(dna_1) == len(dna_2)
    return sum(x != y for x, y in zip(dna_1, dna_2))


def overlap(dna1: str, dna2: str, order: int = 3) -> bool:
    return (dna1 != dna2) and dna1.endswith(dna2[:order])


def overlap_graph(dict_of_seq: Dict[str, str], order: int = 3) -> List[Tuple[str, str]]:
    output = []
    for n1 in dict_of_seq.keys():
        for n2 in dict_of_seq.keys():
            if overlap(dict_of_seq[n1], dict_of_seq[n2], order):
                output.append((n1, n2))
    return output


def find_shared_motif(dict_of_seq: Dict[str, str]) -> str:
    k = len(dict_of_seq)
    generalized_suffix_tree = Tree(dict_of_seq)
    lcss = max(filter(lambda x: x[0] == k,
                      generalized_suffix_tree.common_substrings()),
               key=lambda x: x[1])
    if not lcss:
        return ''

    s = lcss[2]

    return ''.join(map(str, s.S[s.start:s.end]))


def distance_matrix(dict_of_seq: Dict[str, str], order: List[str]):
    n = len(order)
    output = [[0.0 for _ in range(n)] for _ in range(n)]

    for i, specimen1 in enumerate(order):
        dna1 = dict_of_seq[specimen1]
        m = len(dna1)

        for j, specimen2 in enumerate(order[i:]):
            if i == j:
                continue

            dna2 = dict_of_seq[specimen2]
            dist = hamming_distance(dna1, dna2) / m
            output[i][j] = dist
            output[j][i] = dist

    return output


def spectral_convolution(spectra1: List[float], spectra2: List[float]) -> Tuple[int, float]:
    minkowski_diff_counter = Counter(map(lambda x: round(x[0] - x[1], 6), product(spectra1, spectra2)))
    solution = max(minkowski_diff_counter, key=lambda x: minkowski_diff_counter[x])
    return minkowski_diff_counter[solution], abs(solution)
