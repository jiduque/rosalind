import unittest

from rosalind.compare import *

FLOAT_TOLERANCE = 3


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


class TestOverlaps(unittest.TestCase):
    def test_overlap(self):
        x1, x2, y = "AAATAAA", "AAATTTT", True
        self.assertEqual(overlap(x1, x2), y)

        x1, x2, y = "AAATAAA", "TTTTCCC", False
        self.assertEqual(overlap(x1, x2), y)

    def test_overlap_graph(self):
        x = {"Rosalind_0498": "AAATAAA", "Rosalind_2391": "AAATTTT",
             "Rosalind_2323": "TTTTCCC", "Rosalind_0442": "AAATCCC",
             "Rosalind_5013": "GGGTGGG"}
        y = {("Rosalind_0498", "Rosalind_2391"), ("Rosalind_0498", "Rosalind_0442"),
             ("Rosalind_2391", "Rosalind_2323")}

        f_x = overlap_graph(x)
        self.assertEqual(len(y), len(f_x))

        for edge in f_x:
            self.assertIn(edge, y)


class TestDistanceMatrix(unittest.TestCase):
    def test_valid(self):
        x = ["TTTCCATTTA", "GATTCATTTC", "TTTCCATTTT", "GTTCCATTTA"]
        y = [[0.00000, 0.40000, 0.10000, 0.10000],
             [0.40000, 0.00000, 0.40000, 0.30000],
             [0.10000, 0.40000, 0.00000, 0.20000],
             [0.10000, 0.30000, 0.20000, 0.00000]]

        output = distance_matrix(x)

        for i in range(4):
            for j in range(4):
                self.assertAlmostEqual(output[i][j], y[i][j], FLOAT_TOLERANCE)


class TestFindSharedMotif(unittest.TestCase):
    def test_find_shared_motif(self):
        x = {"1": "GATTACA", "2": "TAGACCA", "3": "ATACA"}
        y = {"AC", "CA", "TA"}
        self.assertIn(find_shared_motif(x), y)


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
