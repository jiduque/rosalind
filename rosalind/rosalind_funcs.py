import re
import argparse

from os import path
from typing import List, Tuple
from collections import Counter
from functools import reduce
from itertools import combinations_with_replacement
from difflib import SequenceMatcher

from rosalind.helpers import DNA_COMPLEMENT_MAP, RNA_CODON_MAP, read_fasta_format


def count_nucleotides_in_dna(dna_string: str) -> dict:
    output = {}
    counter = Counter(dna_string)
    for nucleotide in ["A", "C", "G", "T"]:
        output[nucleotide] = counter[nucleotide]
    return output


def transcribe_dna_to_rna(dna_string: str) -> str:
    return "".join(map(lambda x: "U" if x == "T" else x, dna_string))


def get_reverse_complement(dna_string: str) -> str:
    return "".join(map(lambda x: DNA_COMPLEMENT_MAP[x], dna_string.strip()[::-1]))


def number_of_rabbits_with_offspring(number_months: int, number_offspring: int) -> int:
    if 0 < number_months < 3:
        return 1
    i = 3
    r1, r2 = 1, 1
    r3 = r1 + number_offspring * r2

    while i < number_months:
        r1, r2 = r3, r1
        r3 = r1 + number_offspring * r2
        i += 1

    return r3


def gc_content_of_dna(dna_string: str) -> float:
    output = 0
    nucleotide_content = count_nucleotides_in_dna(dna_string)
    sum_nucleotides = sum(nucleotide_content.values())

    if sum_nucleotides:
        output = nucleotide_content["C"] + nucleotide_content["G"]
        output /= sum_nucleotides

    return output


def max_gc_content_in_list_of_dna(dict_of_dna: dict) -> Tuple[str, float]:
    gc_contents = {}

    for specimen, dna in dict_of_dna.items():
        gc_contents[specimen] = gc_content_of_dna(dna)

    return max(gc_contents.items(), key=lambda x: x[1])


def hamming_distance(dna_1: str, dna_2: str) -> int:
    assert len(dna_1) == len(dna_2)
    return sum(x != y for x, y in zip(dna_1, dna_2))


def probability_offspring_has_dominant_allele(homo_dom: int, hetero: int, homo_rec: int) -> float:
    k, m, n = homo_dom, hetero, homo_rec
    N = homo_dom + hetero + homo_rec
    if N == 0:
        return 0

    output = k * (k - 1) + 2 * k * m + 2 * k * n
    output += m * n + 0.75 * m * (m - 1)
    output /= N * (N - 1)

    return output


def translate_rna_to_protein(rna_string: str) -> str:
    n = len(rna_string)
    output, i = "", 0

    while i < n - 3 and RNA_CODON_MAP.get(rna_string[i:i+3]) != "Stop":
        output += RNA_CODON_MAP.get(rna_string[i:i+3])
        i += 3

    return output


def find_motif_in_dna(dna_string: str, motif: str) -> List[int]:
    pattern = rf"(?={motif})"
    return [m.start() for m in re.finditer(pattern, dna_string)]


def get_profile(dict_of_seq: dict) -> dict:
    n = max(len(seq) for seq in dict_of_seq.values())
    output = {}

    for nuc in ["A", "C", "G", "T"]:
        output[nuc] = [0 for _ in range(n)]

    for specimen in dict_of_seq:
        for i in range(n):
            nuc = dict_of_seq[specimen][i]
            output[nuc][i] += 1

    return output


def get_consensus_from_profile(profile_dict: dict) -> str:
    nucleotides = profile_dict.keys()
    n = len(profile_dict["A"])

    output = ""
    for i in range(n):
        max_val = max(
            ((nuc, profile_dict[nuc][i]) for nuc in nucleotides),
            key=lambda x: x[1]
        )
        output += max_val[0]

    return output


def get_consensus_and_profile(dict_of_seq: dict) -> Tuple[str, dict]:
    profile = get_profile(dict_of_seq)
    consensus = get_consensus_from_profile(profile)
    return consensus, profile


def number_of_mortal_rabbits(number_months: int, months_until_death: int) -> int:
    return 3


def overlap(dna1: str, dna2: str, order: int = 3) -> bool:
    return (dna1 != dna2) and (dna1[-order:] == dna2[:order])


# TODO: This solution is not working for some reason
def get_overlap_graph(dict_of_seq: dict, order: int = 3) -> List[Tuple[str, str]]:
    return list(filter(lambda x: overlap(dict_of_seq[x[0]], dict_of_seq[x[1]], order),
                       combinations_with_replacement(dict_of_seq.keys(), 2)
                       )
                )


def expected_offspring(number_individuals: Tuple) -> float:
    probs = (1.0, 1.0, 1.0, 0.75, 0.5, 0)
    return 2 * sum(map(lambda x: x[0] * x[1],
                       zip(number_individuals, probs)
                       )
                   )


# TODO: Might need to return all strings
def _shared_motif(dna1: str, dna2: str) -> str:
    matcher = SequenceMatcher(a=dna1, b=dna2)
    a_index, _, size = matcher.find_longest_match()
    print(a_index, _, size)

    if size == 0:
        return ""

    return dna1[a_index:a_index+size+1]


# TODO: might need to turn this into a back tracking problem
def find_shared_motif(list_of_seq: List[str]) -> str:
    list_of_seq.sort(key=len)
    return reduce(_shared_motif, list_of_seq)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count the nucleotides in a DNA sequence')

    parser.add_argument('file', type=str, nargs=1, help='path to files with DNA sequence')

    args = parser.parse_args()

    if not path.exists(args.file[0]):
        print(f"File does not exists: {args.file[0]}")

    else:
        data = list(read_fasta_format(args.file[0]).values())
        answer = find_shared_motif(data)
        print(answer)
