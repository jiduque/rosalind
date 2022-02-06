import unittest

from rosalind.simulate import *

FLOAT_TOLERANCE = 3


class TestSimulations(unittest.TestCase):
    def test_number_of_rabbits_with_offspring(self):
        x, y = (5, 3), 19
        self.assertEqual(number_of_rabbits_with_offspring(*x), y)

    def test_probability_offspring_has_dominant_allele(self):
        x, y = (2, 2, 2), 0.7833
        self.assertAlmostEqual(probability_offspring_has_dominant_allele(*x), y, FLOAT_TOLERANCE)

    def test_expected_offspring(self):
        x, y = (1, 0, 0, 1, 0, 1), 3.5
        self.assertAlmostEqual(expected_offspring(x), y, FLOAT_TOLERANCE)

    def test_number_of_mortal_rabbits(self):
        x, y = (6, 3), 4
        self.assertEqual(number_of_mortal_rabbits(*x), y)

    def test_probability_at_least_n_organisms(self):
        x, y = (2, 1), 0.684
        self.assertAlmostEqual(probability_at_least_n_organisms(*x), y, FLOAT_TOLERANCE)


if __name__ == '__main__':
    unittest.main()
