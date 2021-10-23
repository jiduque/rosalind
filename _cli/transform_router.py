from pathlib import Path
from typing import List, Union

import rosalind.transform as transform
from rosalind.helpers import read_sequence, fasta_to_dict


class Transform:
    """
    A collection of tools used to edit sequences
    """

    @staticmethod
    def d2r(sequence_file: Path) -> str:
        """Converts a DNA sequence into an RNA sequence"""
        dna_seq = read_sequence(str(sequence_file))
        return transform.transcribe_dna_to_rna(dna_seq)

    @staticmethod
    def rc(sequence_file: Path) -> str:
        """Gets a DNA sequence's reverse compliment"""
        dna_seq = read_sequence(str(sequence_file))
        return transform.reverse_complement(dna_seq)

    @staticmethod
    def r2p(sequence_file: Path) -> str:
        """Converts RNA sequence into protein sequence"""
        rna_seq = read_sequence(str(sequence_file))
        return transform.translate_rna_to_protein(rna_seq)

    @staticmethod
    def d2p(sequence_file: Path, offset: bool = False) -> Union[str, List[str]]:
        """
        Converts DNA sequence into protein sequence.
        If offset is not flagged, it assumes the encoding begins at the start.
        Else it will begin encoding at several locations.
        """
        dna_seq = read_sequence(str(sequence_file))

        if offset:
            return transform.candidate_proteins_from_dna(dna_seq)

        return transform.translate_dna_to_protein(dna_seq)

    @staticmethod
    def splice(fasta_file: Path) -> str:
        """Given a FASTA file, it will remove the smaller specimens from the largest one"""
        fasta_dict = fasta_to_dict(str(fasta_file))
        return transform.splice_dna(fasta_dict)
