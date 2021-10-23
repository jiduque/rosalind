import unittest

from rosalind.analyze import *

FLOAT_TOLERANCE = 3


class TestCountNucleotidesInDNA(unittest.TestCase):
    def test_valid_input(self):
        x = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
        y = {"A": 20, "C": 12, "G": 17, "T": 21}
        self.assertEqual(count_nucleotides_in_dna(x), y)

    def test_empty(self):
        x = " "
        y = {"A": 0, "C": 0, "G": 0, "T": 0}
        self.assertEqual(count_nucleotides_in_dna(x), y)


class TestGCContentInDNA(unittest.TestCase):
    def test_valid_input(self):
        x = "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"
        y = 0.60919540
        self.assertAlmostEqual(gc_content_of_dna(x), y, FLOAT_TOLERANCE)

    def test_empty(self):
        x = " "
        y = 0.0
        self.assertEqual(gc_content_of_dna(x), y)

    def test_valid_input_max_list(self):
        x = {"Rosalind_6404": "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG",
             "Rosalind_5959": "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC",
             "Rosalind_0808": "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"}
        y = ("Rosalind_0808", 0.60919540)

        f_x = max_gc_content_in_list_of_dna(x)

        self.assertEqual(f_x[0], y[0])
        self.assertAlmostEqual(f_x[1], y[1], FLOAT_TOLERANCE)

    def test_empty_max_list(self):
        x = {}
        self.assertRaises(ValueError, max_gc_content_in_list_of_dna, x)


class TestConsensusAndProfile(unittest.TestCase):
    def test_profile_valid(self):
        x = {"Rosalind_1": "ATCCAGCT", "Rosalind_2": "GGGCAACT",
             "Rosalind_3": "ATGGATCT", "Rosalind_4": "AAGCAACC",
             "Rosalind_5": "TTGGAACT", "Rosalind_6": "ATGCCATT",
             "Rosalind_7": "ATGGCACT"}
        y = {"A": [5, 1, 0, 0, 5, 5, 0, 0], "C": [0, 0, 1, 4, 2, 0, 6, 1],
             "G": [1, 1, 6, 3, 0, 1, 0, 0], "T": [1, 5, 0, 0, 0, 1, 1, 6]}

        self.assertEqual(profile(x), y)

    def test_consensus_from_profile_valid(self):
        x = {"A": [5, 1, 0, 0, 5, 5, 0, 0], "C": [0, 0, 1, 4, 2, 0, 6, 1],
             "G": [1, 1, 6, 3, 0, 1, 0, 0], "T": [1, 5, 0, 0, 0, 1, 1, 6]}
        y = "ATGCAACT"

        self.assertEqual(consensus_from_profile(x), y)

    def test_empty_profile(self):
        x = {}
        self.assertRaises(ValueError, profile, x)

    def test_empty_consensus(self):
        x = {}
        self.assertRaises(KeyError, consensus_from_profile, x)


class TestProbabilityOfRandomString(unittest.TestCase):
    def test_valid(self):
        x1 = "ACGATACAA"
        x2 = [0.129, 0.287, 0.423, 0.476, 0.641, 0.742, 0.783]
        y = [-5.737, -5.217, -5.263, -5.360, -5.958, -6.628, -7.009]

        for f_x, y in zip(probability_of_random_string(x1, x2), y):
            self.assertAlmostEqual(f_x, y, FLOAT_TOLERANCE)

    def test_empty(self):
        x1, x2 = '', []
        self.assertEqual(probability_of_random_string(x1, x2), x2)


class TestFailureArray(unittest.TestCase):
    def test_valid(self):
        x = "CAGCATGGTATCACAGCAGAG"
        y = [0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 1, 2, 1, 2, 3, 4, 5, 3, 0, 0]
        self.assertEqual(failure_array(x), y)

    def test_empty(self):
        x, y = '', []
        self.assertEqual(failure_array(x), y)


class TestProteinMass(unittest.TestCase):
    def test_valid(self):
        x, y = "SKADYEK", 821.392
        self.assertAlmostEqual(protein_mass(x), y, FLOAT_TOLERANCE)

    def test_empty(self):
        x, y = '', 0
        self.assertEqual(protein_mass(x), y)


class TestNumberRNAFromProtein(unittest.TestCase):
    def test_valid(self):
        x, y = "MA", 12
        self.assertEqual(number_rna_from_protein(x), y)

    def test_empty(self):
        self.assertRaises(TypeError, number_rna_from_protein, '')


class TestTransitionTranscversionRatio(unittest.TestCase):
    def test_valid(self):
        x1 = "GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT"
        x2 = "TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT"
        y = 1.21428571429

        self.assertAlmostEqual(transition_transversion_ratio(x1, x2), y, FLOAT_TOLERANCE)
    
    def test_empty(self):
        x1, x2 = "", ""
        y = float("inf")
        self.assertEqual(transition_transversion_ratio(x1, x2), y)
    

if __name__ == '__main__':
    unittest.main()
