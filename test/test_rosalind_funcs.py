# TODO: Continue to refactor tests
import unittest

from rosalind.compare import *
from rosalind.infer import *
from rosalind.search import *
from rosalind.simulate import *
from rosalind.transform import *
from rosalind.misc import *

FLOAT_TOLERANCE = 3


class MyTestCase(unittest.TestCase):
    def test_transcribe_dna_to_rna(self):
        x = "GATGGAACTTGACTACGTAAATT"
        y = "GAUGGAACUUGACUACGUAAAUU"
        self.assertEqual(transcribe_dna_to_rna(x), y)

    def test_reverse_compliment_of_dna(self):
        x, y = "AAAACCCGGT", "ACCGGGTTTT"
        self.assertEqual(reverse_complement(x), y)

    def test_number_of_rabbits_with_offspring(self):
        x, y = (5, 3), 19
        self.assertEqual(number_of_rabbits_with_offspring(*x), y)

    def test_probability_offspring_has_dominant_allele(self):
        x, y = (2, 2, 2), 0.7833
        self.assertAlmostEqual(probability_offspring_has_dominant_allele(*x), y, FLOAT_TOLERANCE)

    def test_translate_rna_to_protein(self):
        x = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
        y = "MAMAPRTEINSTRING"
        self.assertEqual(translate_rna_to_protein(x), y)

    def test_find_motif_in_dna(self):
        x, y = ("GATATATGCATATACTT", "ATAT"), [1, 3, 9]
        self.assertEqual(find_motif_in_dna(*x), y)

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

    def test_expected_offspring(self):
        x, y = (1, 0, 0, 1, 0, 1), 3.5
        self.assertAlmostEqual(expected_offspring(x), y, FLOAT_TOLERANCE)

    def test_find_shared_motif(self):
        x = {"1": "GATTACA", "2": "TAGACCA", "3": "ATACA"}
        y = {"AC", "CA", "TA"}
        self.assertIn(find_shared_motif(x), y)

    def test_number_of_mortal_rabbits(self):
        x, y = (6, 3), 4
        self.assertEqual(number_of_mortal_rabbits(*x), y)

    def test_probability_at_least_n_organisms(self):
        x, y = (2, 1), 0.684
        self.assertAlmostEqual(probability_at_least_n_organisms(*x), y, FLOAT_TOLERANCE)

    def test_candidate_proteins_from_dna(self):
        x = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
        y = {"MLLGSFRLIPKETLIQVAGSSPCNLS", "M", "MGMTPRLGLESLLE", "MTPRLGLESLLE"}

        vals = set(candidate_proteins_from_dna(x))

        self.assertFalse(y - vals)

    def test_protein_from_mass_spectrum(self):
        x = [3524.8542, 3710.9335, 3841.974, 3970.0326, 4057.0646]
        y = "WMQS"
        self.assertEqual(protein_from_mass_spectrum(x), y)

    def test_splice_dna(self):
        x = {"Rosalind_10": "ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG",
             "Rosalind_12": "ATCGGTCGAA", "Rosalind_15": "ATCGGTCGAGCGTGT"}
        y = "MVYIADKQHVASREAYGHMFKVCA"

        self.assertEqual(splice_dna(x), y)

    def test_edges_to_complete_the_tree(self):
        x1 = 10
        x2 = [(0, 1), (1, 7), (3, 9), (4, 8), (5, 9), (6, 8)]
        y = 3
        self.assertEqual(edges_to_complete_the_tree(x1, x2), y)

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


if __name__ == '__main__':
    unittest.main()
