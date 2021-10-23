from typing import Tuple, List

from suffix_tree import Tree


def hamming_distance(dna_1: str, dna_2: str) -> int:
    assert len(dna_1) == len(dna_2)
    return sum(x != y for x, y in zip(dna_1, dna_2))


def overlap(dna1: str, dna2: str, order: int = 3) -> bool:
    return (dna1 != dna2) and dna1.endswith(dna2[:order])


def overlap_graph(dict_of_seq: dict, order: int = 3) -> List[Tuple[str, str]]:
    output = []
    for n1 in dict_of_seq.keys():
        for n2 in dict_of_seq.keys():
            if overlap(dict_of_seq[n1], dict_of_seq[n2], order):
                output.append((n1, n2))
    return output


def find_shared_motif(dict_of_seq: dict) -> str:
    k = len(dict_of_seq)
    generalized_suffix_tree = Tree(dict_of_seq)
    lcss = max(filter(lambda x: x[0] == k,
                      generalized_suffix_tree.common_substrings()),
               key=lambda x: x[1])
    if not lcss:
        return ''

    s = lcss[2]

    return ''.join(map(str, s.S[s.start:s.end]))
