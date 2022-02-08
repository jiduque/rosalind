import unittest

from rosalind.infer import *


class TestProteinFromMassSpectrum(unittest.TestCase):
    def test_protein_from_mass_spectrum(self):
        x = [3524.8542, 3710.9335, 3841.974, 3970.0326, 4057.0646]
        y = "WMQS"
        self.assertEqual(protein_from_mass_spectrum(x), y)


class TestSpectrumToProtein(unittest.TestCase):
    def test_valid(self):
        x1 = [445.17838, 115.02694, 186.07931, 314.13789, 317.1198, 215.09061]
        x2 = ["GSDMQS", "VWICN", "IASWMQS", "PVSMGAD"]
        y = [(3, "IASWMQS"), (3, 'GSDMQS')]
        self.assertIn(spectrum_to_protein(x1, x2), y)


if __name__ == '__main__':
    unittest.main()
