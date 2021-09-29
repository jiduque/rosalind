import re
import argparse

from os import path
from math import log, comb, pow
from typing import List, Tuple
from collections import Counter
from functools import reduce
from itertools import combinations_with_replacement

from rosalind.helpers import BinarySearchTree, read_fasta_format
from rosalind.constants import DNA_COMPLEMENT_MAP, PROTEIN_MASS_MAP, RNA_CODON_MAP, PROTEIN_TO_NUM_RNA

from suffix_tree import Tree


def count_nucleotides_in_dna(dna_string: str) -> dict:
    output = {}
    counter = Counter(dna_string)
    for nucleotide in ["A", "C", "G", "T"]:
        output[nucleotide] = counter[nucleotide]
    return output


def transcribe_dna_to_rna(dna_string: str) -> str:
    return "".join(map(lambda x: "U" if x == "T" else x, dna_string))


def reverse_complement(dna_string: str) -> str:
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
    t = homo_dom + hetero + homo_rec
    if t == 0:
        return 0

    output = k * (k - 1) + 2 * k * m + 2 * k * n
    output += m * n + 0.75 * m * (m - 1)
    output /= t * (t - 1)

    return output


def translate_rna_to_protein(rna_string: str) -> str:
    n = len(rna_string)
    output = ""
    for i in range(0, n, 3):
        protein = RNA_CODON_MAP.get(rna_string[i:i+3])
        if protein in ["Stop", None]:
            break
        output += protein

    return output


def translate_dna_to_protein(dna_string: str) -> str:
    return translate_rna_to_protein(transcribe_dna_to_rna(dna_string))


def find_motif_in_dna(dna_string: str, motif: str) -> List[int]:
    pattern = rf"(?={motif})"
    return [m.start() for m in re.finditer(pattern, dna_string)]


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


def number_of_mortal_rabbits(number_months: int, months_until_death: int) -> int:
    assert months_until_death > 2
    if 0 < number_months < 3:
        return 1

    dp = [0 for _ in range(number_months)]
    dp[0], dp[1] = 1, 1
    dp[months_until_death], dp[months_until_death + 1] = -1, -1

    for i in range(2, number_months):
        dp[i] += dp[i-1] + dp[i-2]
        if i > months_until_death + 1:
            dp[i] -= dp[i - months_until_death - 1]

    return dp[number_months - 1]


def overlap(dna1: str, dna2: str, order: int = 3) -> bool:
    return (dna1 != dna2) and dna1.endswith(dna2[:order])


def overlap_graph(dict_of_seq: dict, order: int = 3) -> List[Tuple[str, str]]:
    output = []
    for n1 in dict_of_seq.keys():
        for n2 in dict_of_seq.keys():
            if overlap(dict_of_seq[n1], dict_of_seq[n2], order):
                output.append((n1, n2))
    return output


def expected_offspring(number_individuals: Tuple) -> float:
    probs = (1.0, 1.0, 1.0, 0.75, 0.5, 0)
    return 2 * sum(map(lambda x: x[0] * x[1],
                       zip(number_individuals, probs)
                       )
                   )


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


def probability_at_least_n_organisms(generations: int, n_individuals: int) -> float:
    total_population = 2 ** generations
    return 1 - sum(
        map(lambda x: comb(total_population, x) * pow(0.25, x) * pow(0.75, total_population - x),
            range(n_individuals)
            )
    )


def number_rna_from_protein(protein_seq: str, mod: int = 1E6) -> int:
    return int(3 * reduce(lambda x, y: x * y % mod,
                          map(lambda x: PROTEIN_TO_NUM_RNA[x], protein_seq)
                          ) % mod
               )


def candidate_proteins_from_dna(dna_str: str) -> List[str]:
    pattern = r'(?=(AUG(?:...)*?)(?=UAG|UGA|UAA))'
    rna_str = transcribe_dna_to_rna(dna_str)
    forward = set(map(translate_rna_to_protein, re.findall(pattern, rna_str)))

    rna_str = transcribe_dna_to_rna(reverse_complement(dna_str))
    backward = set(map(translate_rna_to_protein, re.findall(pattern, rna_str)))

    return list(forward.union(backward))


def lexicographic_strings(ordered_alphabet: List[str], length: int) -> List[str]:
    return list(map(lambda x: "".join(x),
                    combinations_with_replacement(ordered_alphabet, length)))


def protein_from_mass_spectrum(mass_data: List[float]) -> str:
    inverted_key_values = list(map(lambda x: (x[1], x[0]), PROTEIN_MASS_MAP.items()))
    search_tree = BinarySearchTree(inverted_key_values)
    mass_diff = list(map(lambda x, y: y - x, mass_data[:-1], mass_data[1:]))
    return "".join(map(lambda x: search_tree.search_closest(x).val, mass_diff))


def splice_dna(dict_of_seq: dict) -> str:
    specimen, dna = max(dict_of_seq.items(), key=lambda x: len(x[1]))

    for key, val in dict_of_seq.items():
        if key != specimen:
            dna = dna.replace(val, '')

    return translate_dna_to_protein(dna)


def transition_transversion_ratio(dna1: str, dna2: str) -> float:
    transition_map = {"A": "G", "G": "A", "C": "T", "T": "C"}
    transitions, transversions = 0, 0

    for x, y in zip(dna1, dna2):
        if transition_map[x] == y:
            transitions += 1
        elif x != y:
            transversions += 1

    return transitions / transversions if transversions else float("inf")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count the nucleotides in a DNA sequence')

    parser.add_argument('file', type=str, nargs=1, help='path to files with DNA sequence')

    args = parser.parse_args()

    if not path.exists(args.file[0]):
        print(f"File does not exists: {args.file[0]}")

    else:
        data = read_fasta_format(args.file[0])
        answer = transition_transversion_ratio(*tuple(data.values()))
        print(answer)

        """
        with open(args.file[0], 'r') as file:
            data = list(map(float, file.read().split()))
            answer = protein_from_mass_spectrum(data)
            print(f"{answer}")

        """
