import unittest

from rosalind.misc import *


class TestMiscFuncs(unittest.TestCase):
    def test_edges_to_complete_the_tree(self):
        x1 = 10
        x2 = [(0, 1), (1, 7), (3, 9), (4, 8), (5, 9), (6, 8)]
        y = 3
        self.assertEqual(edges_to_complete_the_tree(x1, x2), y)


if __name__ == '__main__':
    unittest.main()
