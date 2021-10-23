from pathlib import Path

import rosalind.infer as infer
from rosalind.helpers import read_several_sequences


class Infer:
    """
    A collection of tools used to do inference on sequences, ie a not exact or unique analysis
    """

    @staticmethod
    def protein_from_mass(prefix_sequence_file: Path) -> str:
        """Computes a protein from a prefix spectrum file (diff mass on each line)"""
        mass_data = read_several_sequences(str(prefix_sequence_file))
        mass_data = list(map(float, mass_data))
        return infer.protein_from_mass_spectrum(mass_data)
