"""Test functions for python.energy_functions"""

__author__ = "Felix Simkovic"
__date__ = "05 Dec 2016"

import itertools
import unittest
from ample.modelling.energy_functions import DynamicDistances, RosettaFunctionConstructs


class TestDynamicDistances(unittest.TestCase):
    def test_cutoff(self):
        amino_acids = [
            "A",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "K",
            "L",
            "M",
            "N",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "V",
            "W",
            "Y",
        ]
        for (a1, a2) in itertools.combinations(amino_acids, 2):
            self.assertEqual(DynamicDistances._CB_CB_CUTOFF[a1][a2], DynamicDistances.cutoff(a1, a2))

    def test_percentile(self):
        amino_acids = [
            "A",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "K",
            "L",
            "M",
            "N",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "V",
            "W",
            "Y",
        ]
        for (a1, a2) in itertools.combinations(amino_acids, 2):
            self.assertEqual(DynamicDistances._CB_CB_PERCENT[a1][a2], DynamicDistances.percentile(a1, a2))


class TestRosettaFunctionConstructs(unittest.TestCase):
    def test_gaussian(self):
        optd = {'atom1': 'CA', 'res1_seq': 1, 'atom2': 'CA', 'res2_seq': 100, 'mean': 6.0, 'stddev': 5.0}
        restraint = RosettaFunctionConstructs().GAUSSIAN.format(**optd)
        ref = 'AtomPair CA    1 CA  100 GAUSSIANFUNC 6.00 5.00 TAG'
        self.assertEqual(restraint, ref, "Incorrrect restraint: {} -> {}".format(restraint, ref))


if __name__ == "__main__":
    unittest.main(verbosity=2)
