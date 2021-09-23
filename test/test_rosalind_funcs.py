import unittest

from rosalind.rosalind_funcs import *


class MyTestCase(unittest.TestCase):
    def test_count_nucleotides_in_dna(self):
        x = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
        y = {"A": 20, "C": 12, "G": 17, "T": 21}
        self.assertEqual(count_nucleotides_in_dna(x), y)

    def test_transcribe_dna_to_rna(self):
        x = "GATGGAACTTGACTACGTAAATT"
        y = "GAUGGAACUUGACUACGUAAAUU"
        self.assertEqual(transcribe_dna_to_rna(x), y)

    def test_get_reverse_compliment_of_dna(self):
        x, y = "AAAACCCGGT", "ACCGGGTTTT"
        self.assertEqual(get_reverse_complement(x), y)

    def test_number_of_rabbits_with_offspring(self):
        x, y = (5, 3), 19
        self.assertEqual(number_of_rabbits_with_offspring(*x), y)

    def test_gc_content_in_dna(self):
        x = "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"
        y = 0.60919540
        self.assertAlmostEqual(gc_content_of_dna(x), y, 5)

    def test_max_gc_content_in_list_of_dna(self):
        x = {"Rosalind_6404": "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG",
             "Rosalind_5959": "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC",
             "Rosalind_0808": "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"}
        y = ("Rosalind_0808", 0.60919540)

        f_x = max_gc_content_in_list_of_dna(x)

        self.assertEqual(f_x[0], y[0])
        self.assertAlmostEqual(f_x[1], y[1], 3)

    def test_hamming_distance(self):
        x1, x2 = "GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"
        y = 7
        self.assertEqual(hamming_distance(x1, x2), y)

    def test_probability_offspring_has_dominant_allele(self):
        x, y = (2, 2, 2), 0.7833
        self.assertAlmostEqual(probability_offspring_has_dominant_allele(*x), y, 3)

    def test_translate_rna_to_protein(self):
        x = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
        y = "MAMAPRTEINSTRING"
        self.assertEqual(translate_rna_to_protein(x), y)

    def test_find_motif_in_dna(self):
        x, y = ("GATATATGCATATACTT", "ATAT"), [1, 3, 9]
        self.assertEqual(find_motif_in_dna(*x), y)

    def test_get_profile(self):
        x = {"Rosalind_1": "ATCCAGCT", "Rosalind_2": "GGGCAACT",
             "Rosalind_3": "ATGGATCT", "Rosalind_4": "AAGCAACC",
             "Rosalind_5": "TTGGAACT", "Rosalind_6": "ATGCCATT",
             "Rosalind_7": "ATGGCACT"}
        y = {"A": [5, 1, 0, 0, 5, 5, 0, 0], "C": [0, 0, 1, 4, 2, 0, 6, 1],
             "G": [1, 1, 6, 3, 0, 1, 0, 0], "T": [1, 5, 0, 0, 0, 1, 1, 6]}

        self.assertEqual(get_profile(x), y)

    def test_get_consensus_from_profile(self):
        x = {"A": [5, 1, 0, 0, 5, 5, 0, 0], "C": [0, 0, 1, 4, 2, 0, 6, 1],
             "G": [1, 1, 6, 3, 0, 1, 0, 0], "T": [1, 5, 0, 0, 0, 1, 1, 6]}
        y = "ATGCAACT"

        self.assertEqual(get_consensus_from_profile(x), y)

    def test_overlap(self):
        x1, x2, y = "AAATAAA", "AAATTTT", True
        self.assertEqual(overlap(x1, x2), y)

        x1, x2, y = "AAATAAA", "TTTTCCC", False
        self.assertEqual(overlap(x1, x2), y)

    def test_get_overlap_graph(self):
        x = {"Rosalind_0498": "AAATAAA", "Rosalind_2391": "AAATTTT",
             "Rosalind_2323": "TTTTCCC", "Rosalind_0442": "AAATCCC",
             "Rosalind_5013": "GGGTGGG"}
        answer_list = [("Rosalind_0498", "Rosalind_2391"), ("Rosalind_0498", "Rosalind_0442"),
                       ("Rosalind_2391", "Rosalind_2323")]
        y = set()
        for edge in answer_list:
            y.add(edge)

        for edge in get_overlap_graph(x):
            self.assertIn(edge, y)

    def test_expected_offspring(self):
        x, y = (1, 0, 0, 1, 0, 1), 3.5
        self.assertAlmostEqual(expected_offspring(x), y, 3)

    def test_find_shared_motif(self):
        x = ["GATTACA", "TAGACCA", "ATACA"]
        y = {"AC", "CA", "TA"}
        self.assertIn(find_shared_motif(x), y)


if __name__ == '__main__':
    unittest.main()
