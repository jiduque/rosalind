import re

from typing import List, Tuple
from concurrent.futures import ThreadPoolExecutor

from rosalind.analyze import gc_content_of_dna
from rosalind.constants import UNIPROT_DB
from rosalind.transform import reverse_complement

import requests


def find_max_gc(dict_of_dna: dict) -> Tuple[str, float]:
    gc_contents = {}

    for specimen, dna in dict_of_dna.items():
        gc_contents[specimen] = gc_content_of_dna(dna)

    return max(gc_contents.items(), key=lambda x: x[1])


def find_motif_in_dna(dna_string: str, motif: str) -> List[int]:
    pattern = rf"(?={motif})"
    return [m.start() for m in re.finditer(pattern, dna_string)]


def n_glycosylation_motif(protein_str: str) -> List:
    pattern = r'(?=(N[^P](S|T)[^P]))'
    url = f"{UNIPROT_DB}/{protein_str}.fasta"
    response = requests.get(url)
    if response.ok:
        protein_sequence = "".join(response.text.split("\n")[1:])
        output = [m.start(0) + 1 for m in re.finditer(pattern, protein_sequence)]
        return output
    return []


def list_n_glycosylation_motif(list_of_proteins: List[str]) -> dict:
    output = {}

    with ThreadPoolExecutor() as executor:
        results = executor.map(n_glycosylation_motif, list_of_proteins)
        for x, y in zip(list_of_proteins, results):
            if y:
                output[x] = y

    return output


# TODO: still need to work out the logic for this one
def find_restriction_sites(dna_str: str) -> List[Tuple[int, int]]:
    """
    http://rosalind.info/problems/revp/
    :param dna_str:
    :return:
    """

    # turn reverse_comp into suffix tree
    reverse_comp = reverse_complement(dna_str)
    output = []
    n = len(dna_str)
    length_tracker = [0 for _ in range(n)]

    for i in range(n):
        if dna_str[i] == reverse_comp[i]:
            length_tracker[i] = 1 + int(i > 0) * length_tracker[i - 1]
    return list(map(lambda x: (x[0] - x[1] + 1, x[1]),
                    filter(lambda x: 3 < x[1] < 13,
                           enumerate(output)
                           )
                    )
                )

