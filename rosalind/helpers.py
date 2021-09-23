import pickle

from pathlib import Path


CURRENT_LOCATION = Path(__file__)
DATA_LOCATION = CURRENT_LOCATION.parents[1] / 'data'

DNA_COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}

with (DATA_LOCATION / 'rna_codon_table.pkl').open('rb') as handle:
    RNA_CODON_MAP = pickle.load(handle)


def read_fasta_format(file_path: str) -> dict:
    output = {}
    with open(file_path, 'r') as file:
        for pair in file.read().split(">")[1:]:
            data = pair.split()
            output[data[0]] = "".join(data[1:])

    return output
