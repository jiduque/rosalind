from pathlib import Path
from typing import Dict, List

import rosalind.analyze as analyze

from rosalind.helpers import fasta_to_dict, read_sequence


class Analyze:
    """
    Collection of tools to analyze a sequence or sequences. Not to compare sequences, but to summarize them
    """

    @staticmethod
    def count_nucleo(sequence_file: Path) -> Dict:
        """Counts the nucleotides in a DNA string"""
        dna_string = read_sequence(str(sequence_file))
        return analyze.count_nucleotides_in_dna(dna_string)

    @staticmethod
    def gc_content(sequence_file: Path) -> float:
        """Calculates the GC percentage of a DNA string"""
        dna_string = read_sequence(str(sequence_file))
        return analyze.gc_content_of_dna(dna_string)

    @staticmethod
    def profile(fasta_file: Path) -> Dict:
        """Computes the profile of several specimen in FASTA file"""
        fasta_dict = fasta_to_dict(str(fasta_file))
        return analyze.profile(fasta_dict)

    @staticmethod
    def consensus(fasta_file: Path) -> str:
        """Computes the consensus of several specimen in FASTA file"""
        fasta_dict = fasta_to_dict(str(fasta_file))
        con, _ = analyze.consensus_and_profile(fasta_dict)
        return con

    @staticmethod
    def fail_array(sequence_file: Path) -> List[int]:
        """Computes the failure array of a file a sequence"""
        seq = read_sequence(str(sequence_file))
        return analyze.failure_array(seq)

    @staticmethod
    def protein_mass(sequence_file: Path) -> float:
        """Computes the mass of a protein sequence"""
        protein_str = read_sequence(str(sequence_file))
        return analyze.protein_mass(protein_str)

    @staticmethod
    def num_sources(sequence_file: Path, mod: int = int(1E6)) -> int:
        """Calculates the number of RNA sequences that map to this protein string (modulo an optional number)"""
        protein_seq = read_sequence(str(sequence_file))
        return analyze.number_rna_from_protein(protein_seq, mod)