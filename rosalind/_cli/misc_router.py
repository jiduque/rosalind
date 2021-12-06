from pathlib import Path
from typing import List

import rosalind.misc as misc
from rosalind.helpers import read_several_sequences


class Misc:
    """
    Miscellaneous tools used for solving the problems on rosalind
    """

    @staticmethod
    def lexi_kmer(file_path: Path) -> List[str]:
        symbols, num = read_several_sequences(str(file_path))[:2]

        symbols = symbols.strip().split()
        num = int(num)

        return misc.lexicographic_strings(symbols, num)
