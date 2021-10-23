import pickle

from pathlib import Path


CURRENT_LOCATION = Path(__file__)
DATA_LOCATION = CURRENT_LOCATION.parents[1] / 'data'

DNA_COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}

with (DATA_LOCATION / 'rna_codon_table.pkl').open('rb') as handle:
    RNA_CODON_MAP = pickle.load(handle)
    PROTEIN_TO_NUM_RNA = {}

    for key, value in RNA_CODON_MAP.items():
        if value not in PROTEIN_TO_NUM_RNA:
            PROTEIN_TO_NUM_RNA[value] = 0
        PROTEIN_TO_NUM_RNA[value] += 1

with (DATA_LOCATION / 'protein_mass.pkl').open('rb') as handle:
    PROTEIN_MASS_MAP = pickle.load(handle)

UNIPROT_DB = "http://www.uniprot.org/uniprot"
