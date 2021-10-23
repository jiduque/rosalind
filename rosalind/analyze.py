
from math import log
from collections import Counter
from functools import reduce

from typing import List, Tuple

from rosalind.constants import PROTEIN_MASS_MAP, PROTEIN_TO_NUM_RNA


def count_nucleotides_in_dna(dna_string: str) -> dict:
    output = {}
    counter = Counter(dna_string)
    for nucleotide in ["A", "C", "G", "T"]:
        output[nucleotide] = counter[nucleotide]
    return output


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


def profile(dict_of_seq: dict) -> dict:
    n = max(len(seq) for seq in dict_of_seq.values())
    output = {}

    for nuc in ["A", "C", "G", "T"]:
        output[nuc] = [0 for _ in range(n)]

    for specimen in dict_of_seq:
        for i in range(n):
            nuc = dict_of_seq[specimen][i]
            output[nuc][i] += 1

    return output


def consensus_from_profile(profile_dict: dict) -> str:
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


def consensus_and_profile(dict_of_seq: dict) -> Tuple[str, dict]:
    prof = profile(dict_of_seq)
    consensus = consensus_from_profile(prof)
    return consensus, prof


def probability_of_random_string(seq_to_analyze: str, list_of_probs: List[float]) -> List[float]:
    counter = Counter(seq_to_analyze)
    p_at = counter["A"] + counter["T"]
    p_cg = counter["C"] + counter["G"]

    return list(map(lambda x: p_at * log((1-x) / 2, 10) + p_cg * log(x/2, 10), list_of_probs))


def failure_array(str_seq: str) -> List[int]:
    n = len(str_seq)
    output = [0 for _ in range(n)]

    i, last_match_length = 1, 0

    while i < n:
        if str_seq[i] != str_seq[last_match_length] and last_match_length:
            last_match_length = output[last_match_length - 1]
            continue

        if str_seq[i] == str_seq[last_match_length]:
            last_match_length += 1
            output[i] = last_match_length

        i += 1

    return output


def protein_mass(protein_seq: str) -> float:
    return sum(map(lambda x: PROTEIN_MASS_MAP[x], protein_seq))


def number_rna_from_protein(protein_seq: str, mod: int = 1E6) -> int:
    return int(3 * reduce(lambda x, y: x * y % mod,
                          map(lambda x: PROTEIN_TO_NUM_RNA[x], protein_seq)
                          ) % mod
               )


def transition_transversion_ratio(dna1: str, dna2: str) -> float:
    transition_map = {"A": "G", "G": "A", "C": "T", "T": "C"}
    transitions, transversions = 0, 0

    for x, y in zip(dna1, dna2):
        if transition_map[x] == y:
            transitions += 1
        elif x != y:
            transversions += 1

    return transitions / transversions if transversions else float("inf")
