import unittest

from rosalind.search import *

FLOAT_TOLERANCE = 3

# TODO: finish adding the test cases from for the search module


class FindMaxGC(unittest.TestCase):
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


if __name__ == '__main__':
    unittest.main()
