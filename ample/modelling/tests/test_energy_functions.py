"""Test functions for python.energy_functions"""

__author__ = "Felix Simkovic"
__date__ = "05 Dec 2016"

import itertools
import unittest
from ample.modelling.energy_functions import DynamicDistances


class TestDynamicDistances(unittest.TestCase):

    def test_cutoff(self):
        amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        for (a1, a2) in itertools.combinations(amino_acids, 2):
            self.assertEqual(
                DynamicDistances._CB_CB_CUTOFF[a1][a2],
                DynamicDistances.cutoff(a1, a2)
            )

    def test_percentile(self):
        amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        for (a1, a2) in itertools.combinations(amino_acids, 2):
            self.assertEqual(
                DynamicDistances._CB_CB_PERCENT[a1][a2],
                DynamicDistances.percentile(a1, a2)
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
