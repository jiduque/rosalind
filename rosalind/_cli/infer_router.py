from pathlib import Path
from typing import Tuple

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

    @staticmethod
    def protein_from_spec(data_file: Path) -> Tuple[int, str]:
        """
        Finds protein from a given list with highest multiplicity with mass spectrum.
        First line of file is the number of proteins, then the protein strings, and finally the spectrum.
        All values should be in a new line. 
        """
        data = read_several_sequences(str(data_file))
        n_proteins = int(data[0])
        proteins = data[1:n_proteins+1]
        mass_data = list(map(float, data[n_proteins+1:]))
        return infer.spectrum_to_protein(mass_data, proteins)
