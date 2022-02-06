from pathlib import Path
from typing import List, Tuple

import rosalind.compare as compare
from rosalind.helpers import read_several_sequences, fasta_to_dict, fasta_to_list


class Compare:
    """
    A collection of tools to compare two or more sequences
    """

    @staticmethod
    def tt_ratio(sequence_file: Path) -> float:
        """Computes the transition-transversion ratio of two DNA sequences in FASTA file"""
        dna_seq = fasta_to_list(str(sequence_file))[:2]
        dna1, dna2 = dna_seq[0][1], dna_seq[1][1]
        return compare.transition_transversion_ratio(dna1, dna2)

    @staticmethod
    def hamming(sequence_file: Path) -> int:
        """Computes the Hamming distance between two sequences of the same length"""
        seqs = read_several_sequences(str(sequence_file))[:2]
        return compare.hamming_distance(*seqs)

    @staticmethod
    def overlap(fasta_file: Path, order: int = 3) -> List[Tuple[str, str]]:
        """Computes the overlap graph (adjacency list) of a giver order of two specimens in FASTA file"""
        fasta_dict = fasta_to_dict(str(fasta_file))
        return compare.overlap_graph(fasta_dict, order)

    @staticmethod
    def shared_motif(fasta_file: Path) -> str:
        """Computes the largest shared motif of several specimens in a FASTA file"""
        fasta_dict = fasta_to_dict(str(fasta_file))
        return compare.find_shared_motif(fasta_dict)

    @staticmethod
    def distance_matrix(fasta_file: Path) -> List[List[float]]:
        """Computest the distance matrix of several sequences in a FASTA file"""
        dna_list = list(map(
            lambda x: x[1],
            fasta_to_list(str(fasta_file))
        ))
        return compare.distance_matrix(dna_list)

    @staticmethod
    def spc_conv(sequence_file: Path) -> Tuple[int, float]:
        """Computes the largest multiplicity of the spectral convolution of two spectra in a file"""
        spc_lines = read_several_sequences(str(sequence_file))[:2]
        spectra = tuple(map(lambda line: list(map(float, line.split())), spc_lines))
        return compare.spectral_convolution(*spectra)
