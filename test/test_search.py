import unittest

from rosalind.search import *

FLOAT_TOLERANCE = 3


class TestFindMaxGC(unittest.TestCase):
    def test_valid_input_max_list(self):
        x = {"Rosalind_6404": "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG",
             "Rosalind_5959": "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC",
             "Rosalind_0808": "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"}
        y = ("Rosalind_0808", 0.60919540)

        f_x = find_max_gc(x)

        self.assertEqual(f_x[0], y[0])
        self.assertAlmostEqual(f_x[1], y[1], FLOAT_TOLERANCE)

    def test_empty_max_list(self):
        x = {}
        self.assertRaises(ValueError, find_max_gc, x)


class TestFindMotifInDNA(unittest.TestCase):
    def test_find_motif_in_dna(self):
        x, y = ("GATATATGCATATACTT", "ATAT"), [1, 3, 9]
        self.assertEqual(find_motif_in_dna(*x), y)


class TestNGlycosulationMotif(unittest.TestCase):
    def test_n_glycosylation_motif(self):
        x = "B5ZC00"
        y = [85, 118, 142, 306, 395]
        self.assertEqual(n_glycosylation_motif(x), y)

    def test_list_n_glycosylation_motif(self):
        x = ["A2Z669", "B5ZC00", "P07204_TRBM_HUMAN", "P20840_SAG1_YEAST"]
        y = {"B5ZC00": [85, 118, 142, 306, 395],
             "P07204_TRBM_HUMAN": [47, 115, 116, 382, 409],
             "P20840_SAG1_YEAST": [79, 109, 135, 248, 306, 348, 364, 402, 485, 501, 614]
             }
        f_x = list_n_glycosylation_motif(x)

        for key in y:
            self.assertIn(key, f_x)
            self.assertEqual(f_x[key], y[key])


class TestRestrictionSite(unittest.TestCase):
    # TODO: Implement test and code
    pass


if __name__ == '__main__':
    unittest.main()
