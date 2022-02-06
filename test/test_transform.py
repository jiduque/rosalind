import unittest

from rosalind.transform import *


class TestTranscribe(unittest.TestCase):
    def test_valid(self):
        x = "GATGGAACTTGACTACGTAAATT"
        y = "GAUGGAACUUGACUACGUAAAUU"
        self.assertEqual(transcribe_dna_to_rna(x), y)

    def test_empty(self):
        x, y = "", ""
        self.assertEqual(transcribe_dna_to_rna(x), y)


class TestReverseCompliment(unittest.TestCase):
    def test_reverse_compliment_of_dna(self):
        x, y = "AAAACCCGGT", "ACCGGGTTTT"
        self.assertEqual(reverse_complement(x), y)


class TestTranslateRNAToProtein(unittest.TestCase):
    def test_translate_rna_to_protein(self):
        x = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
        y = "MAMAPRTEINSTRING"
        self.assertEqual(translate_rna_to_protein(x), y)


class TestCandidateProteins(unittest.TestCase):
    def test_candidate_proteins_from_dna(self):
        x = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
        y = {"MLLGSFRLIPKETLIQVAGSSPCNLS", "M", "MGMTPRLGLESLLE", "MTPRLGLESLLE"}

        vals = set(candidate_proteins_from_dna(x))

        self.assertFalse(y - vals)


class TestDNASplice(unittest.TestCase):
    def test_splice_dna(self):
        x = {"Rosalind_10": "ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG",
             "Rosalind_12": "ATCGGTCGAA", "Rosalind_15": "ATCGGTCGAGCGTGT"}
        y = "MVYIADKQHVASREAYGHMFKVCA"

        self.assertEqual(splice_dna(x), y)


if __name__ == '__main__':
    unittest.main()
