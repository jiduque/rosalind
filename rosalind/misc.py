
from typing import List, Tuple
from itertools import product

from rosalind._data_structures import UnionFind


def lexicographic_strings(ordered_alphabet: List[str], length: int) -> List[str]:
    return list(map(lambda x: "".join(x),
                    product(ordered_alphabet, repeat=length)))


def edges_to_complete_the_tree(number_of_nodes: int, list_of_edges: List[Tuple[int, int]]) -> int:
    union_find = UnionFind(number_of_nodes)
    for p, q in list_of_edges:
        union_find.union(p, q)
    return union_find.n_connected_sets - 1

