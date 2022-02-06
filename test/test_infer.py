import unittest

from rosalind.infer import *


class TestProteinFromMassSpectrum(unittest.TestCase):
    def test_protein_from_mass_spectrum(self):
        x = [3524.8542, 3710.9335, 3841.974, 3970.0326, 4057.0646]
        y = "WMQS"
        self.assertEqual(protein_from_mass_spectrum(x), y)


if __name__ == '__main__':
    unittest.main()
