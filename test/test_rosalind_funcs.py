import unittest

from rosalind.rosalind_funcs import *

FLOAT_TOLERANCE = 3


class MyTestCase(unittest.TestCase):
    def test_count_nucleotides_in_dna(self):
        x = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
        y = {"A": 20, "C": 12, "G": 17, "T": 21}
        self.assertEqual(count_nucleotides_in_dna(x), y)

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

    def test_gc_content_in_dna(self):
        x = "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"
        y = 0.60919540
        self.assertAlmostEqual(gc_content_of_dna(x), y, FLOAT_TOLERANCE)

    def test_max_gc_content_in_list_of_dna(self):
        x = {"Rosalind_6404": "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG",
             "Rosalind_5959": "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC",
             "Rosalind_0808": "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"}
        y = ("Rosalind_0808", 0.60919540)

        f_x = max_gc_content_in_list_of_dna(x)

        self.assertEqual(f_x[0], y[0])
        self.assertAlmostEqual(f_x[1], y[1], FLOAT_TOLERANCE)

    def test_hamming_distance(self):
        x1, x2 = "GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"
        y = 7
        self.assertEqual(hamming_distance(x1, x2), y)

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

    def test_profile(self):
        x = {"Rosalind_1": "ATCCAGCT", "Rosalind_2": "GGGCAACT",
             "Rosalind_3": "ATGGATCT", "Rosalind_4": "AAGCAACC",
             "Rosalind_5": "TTGGAACT", "Rosalind_6": "ATGCCATT",
             "Rosalind_7": "ATGGCACT"}
        y = {"A": [5, 1, 0, 0, 5, 5, 0, 0], "C": [0, 0, 1, 4, 2, 0, 6, 1],
             "G": [1, 1, 6, 3, 0, 1, 0, 0], "T": [1, 5, 0, 0, 0, 1, 1, 6]}

        self.assertEqual(profile(x), y)

    def test_consensus_from_profile(self):
        x = {"A": [5, 1, 0, 0, 5, 5, 0, 0], "C": [0, 0, 1, 4, 2, 0, 6, 1],
             "G": [1, 1, 6, 3, 0, 1, 0, 0], "T": [1, 5, 0, 0, 0, 1, 1, 6]}
        y = "ATGCAACT"

        self.assertEqual(consensus_from_profile(x), y)

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
        x = ["GATTACA", "TAGACCA", "ATACA"]
        y = {"AC", "CA", "TA"}
        self.assertIn(find_shared_motif(x), y)

    def test_number_of_mortal_rabbits(self):
        x, y = (6, 3), 4
        self.assertEqual(number_of_mortal_rabbits(*x), y)

    def test_probability_of_random_string(self):
        x1 = "ACGATACAA"
        x2 = [0.129, 0.287, 0.423, 0.476, 0.641, 0.742, 0.783]
        y = [-5.737, -5.217, -5.263, -5.360, -5.958, -6.628, -7.009]

        for f_x, y in zip(probability_of_random_string(x1, x2), y):
            self.assertAlmostEqual(f_x, y, FLOAT_TOLERANCE)

    def test_failure_array(self):
        x = "CAGCATGGTATCACAGCAGAG"
        y = [0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 1, 2, 1, 2, 3, 4, 5, 3, 0, 0]
        self.assertEqual(failure_array(x), y)

    def test_protein_mass(self):
        x, y = "SKADYEK", 821.392
        self.assertAlmostEqual(protein_mass(x), y, FLOAT_TOLERANCE)

    def test_probability_at_least_n_organisms(self):
        x, y = (2, 1), 0.684
        self.assertAlmostEqual(probability_at_least_n_organisms(*x), y, FLOAT_TOLERANCE)

    def test_number_rna_from_protein(self):
        x, y = "MA", 12
        self.assertEqual(number_rna_from_protein(x), y)

    def test_candidate_proteins_from_dna(self):
        x = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
        y = {"MLLGSFRLIPKETLIQVAGSSPCNLS", "M", "MGMTPRLGLESLLE", "MTPRLGLESLLE"}

        vals = set(candidate_proteins_from_dna(x))

        self.assertFalse(y - vals)

    def test_protein_from_mass_spectrum(self):
        x = [3524.8542, 3710.9335, 3841.974, 3970.0326, 4057.0646]
        y = "WMQS"
        self.assertEqual(protein_from_mass_spectrum(x), y)


if __name__ == '__main__':
    unittest.main()
