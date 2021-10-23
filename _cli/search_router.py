from pathlib import Path
from typing import Tuple, List

import rosalind.search as search
from rosalind.helpers import fasta_to_dict, read_several_sequences


class Search:
    """
    A collection of tools used to search for substrings or the largest/smallest out of a lot
    """

    @staticmethod
    def max_gc(fasta_file: Path) -> Tuple[str, float]:
        """Finds specimen with largest GC content in FASTA file"""
        fasta_dict = fasta_to_dict(str(fasta_file))
        return search.max_gc_content_in_list_of_dna(fasta_dict)

    @staticmethod
    def motif(sequence_file: Path) -> List[int]:
        """Finds location(s) of a motif in DNA. File is two sequences, smaller is motif"""
        sequences = read_several_sequences(str(sequence_file))[:2]
        sequences.sort(key=len)
        return search.find_motif_in_dna(*sequences)

    @staticmethod
    def n_glyco(protein_str: str) -> List:
        """Downloads protein sequence from UNIPROT and finds location(s) of n-glycosylation in sequence """
        return search.n_glycosylation_motif(protein_str)
