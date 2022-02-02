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


class TestSpectralConvolution(unittest.TestCase):
    def test_valid(self):
        s1 = [186.07931, 287.12699, 548.20532, 580.18077, 681.22845, 706.27446, 782.27613, 968.35544, 968.35544]
        s2 = [101.04768, 158.06914, 202.09536, 318.09979, 419.14747, 463.17369]
        y = (3, 85.03163)

        output = spectral_convolution(s1, s2)

        self.assertEqual(output[0], y[0])
        self.assertAlmostEqual(output[1], y[1], FLOAT_TOLERANCE)

    def test_empty(self):
        s1 = []
        s2 = [1, 3, 19.2]

        self.assertRaises(ValueError, spectral_convolution, s1, s2)


if __name__ == '__main__':
    unittest.main()
