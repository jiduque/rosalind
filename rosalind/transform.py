import re
from typing import List

from rosalind.constants import RNA_CODON_MAP, DNA_COMPLEMENT_MAP


def transcribe_dna_to_rna(dna_string: str) -> str:
    return "".join(map(lambda x: "U" if x == "T" else x, dna_string))


def reverse_complement(dna_string: str) -> str:
    return "".join(map(lambda x: DNA_COMPLEMENT_MAP[x], dna_string.strip()[::-1]))


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


def candidate_proteins_from_dna(dna_str: str) -> List[str]:
    pattern = r'(?=(AUG(?:...)*?)(?=UAG|UGA|UAA))'
    rna_str = transcribe_dna_to_rna(dna_str)
    forward = set(map(translate_rna_to_protein, re.findall(pattern, rna_str)))

    rna_str = transcribe_dna_to_rna(reverse_complement(dna_str))
    backward = set(map(translate_rna_to_protein, re.findall(pattern, rna_str)))

    return list(forward.union(backward))


def splice_dna(dict_of_seq: dict) -> str:
    specimen, dna = max(dict_of_seq.items(), key=lambda x: len(x[1]))

    for key, val in dict_of_seq.items():
        if key != specimen:
            dna = dna.replace(val, '')

    return translate_dna_to_protein(dna)
