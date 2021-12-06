import unittest

from rosalind.compare import *

FLOAT_TOLERANCE = 3

# TODO: Finish adding tests for the compare module


class TestTransitionTransversionRatio(unittest.TestCase):
    def test_valid(self):
        x1 = "GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT"
        x2 = "TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT"
        y = 1.21428571429

        self.assertAlmostEqual(transition_transversion_ratio(x1, x2), y, FLOAT_TOLERANCE)

    def test_empty(self):
        x1, x2 = "", ""
        y = float("inf")
        self.assertEqual(transition_transversion_ratio(x1, x2), y)


class TestHammingDistance(unittest.TestCase):
    def test_valid(self):
        x1, x2 = "GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"
        y = 7
        self.assertEqual(hamming_distance(x1, x2), y)

    def test_empty(self):
        x1, x2 = "", ""
        y = 0
        self.assertEqual(hamming_distance(x1, x2), y)


if __name__ == '__main__':
    unittest.main()
