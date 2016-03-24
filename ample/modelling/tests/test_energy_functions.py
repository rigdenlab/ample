"""Test functions for python.energy_functions"""

import itertools
import unittest
from ample.modelling import energy_functions

class Test(unittest.TestCase):
    
    def test_cbcbdata(self):
        ########################################################################
        # Test Case 1
        amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        for (a1, a2) in itertools.product(amino_acids, repeat=2):
            self.assertIn(a1, energy_functions.cbcbcutoff)
            self.assertIn(a1, energy_functions.cbcbpercent)
            self.assertIn(a2, energy_functions.cbcbcutoff[a1])
            self.assertIn(a2, energy_functions.cbcbpercent[a2])

    def test_BOUNDED_gremlin(self):
        ########################################################################
        # Test Case 1
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13, 
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 2, 'scalar_score': 0.55555}
        l = energy_functions.BOUNDED_gremlin(contact)
        ref = "AtomPair CB 13 CB 50 SCALARWEIGHTEDFUNC 0.556 BOUNDED 0 11.810 1 0.5"
        self.assertEqual(ref, l)

    def test_FADE(self):
        ########################################################################
        # Test Case 1
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13, 
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 1, 'scalar_score': 0.55555}
        l = energy_functions.FADE(contact)
        ref = "AtomPair CB 13 CB 50 FADE -10 19 10 -15.00 0"
        self.assertEqual(ref, l)
        
        ########################################################################
        # Test Case 2
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13, 
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 2, 'scalar_score': 0.55555}
        l = energy_functions.FADE(contact)
        ref = "AtomPair CB 13 CB 50 FADE -10 19 10 -30.00 0"
        self.assertEqual(ref, l)

    def test_SIGMOID_gremlin(self):
        ########################################################################
        # Test Case 1
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13,
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 1, 'scalar_score': 0.55555}
        l = energy_functions.SIGMOID_gremlin(contact)
        ref = "AtomPair CB 13 CB 50 SCALARWEIGHTEDFUNC 0.556 SUMFUNC 2 SIGMOID 10.250 1.560 CONSTANTFUNC -0.5"
        self.assertEqual(ref, l)

    def test_BOUNDED_default(self):
        ########################################################################
        # Test Case 1
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13,
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 1, 'scalar_score': 0.55555}
        l = energy_functions.BOUNDED_default(contact)
        ref = "AtomPair CB 13 CB 50 BOUNDED 0.00 8.00 1 0.5 #"
        self.assertEqual(ref, l)

    def test_FADE_default(self):
        ########################################################################
        # Test Case 1
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13,
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 1, 'scalar_score': 0.55555}
        l = energy_functions.FADE_default(contact)
        ref = "AtomPair CB 13 CB 50 FADE -10 19 10 -15.00 0"
        self.assertEqual(ref, l)

        ########################################################################
        # Test Case 2
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13,
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 2, 'scalar_score': 0.55555}
        l = energy_functions.FADE_default(contact)
        ref = "AtomPair CB 13 CB 50 FADE -10 19 10 -15.00 0"
        self.assertEqual(ref, l)

    def test_SIGMOID_default(self):
        ########################################################################
        # Test Case 1
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13,
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 1, 'scalar_score': 0.55555}
        l = energy_functions.SIGMOID_default(contact)
        ref = "AtomPair CB 13 CB 50 SIGMOID 8.00 1.00 #ContactMap: 0.33"
        self.assertEqual(ref, l)

    def test_BOUNDED_scalarweighted(self):
        ########################################################################
        # Test Case 1
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13,
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 1, 'scalar_score': 0.55555}
        l = energy_functions.BOUNDED_scalarweighted(contact)
        ref = "AtomPair CB 13 CB 50 SCALARWEIGHTEDFUNC 0.556 BOUNDED 0.00 8.00 1 0.5 #"
        self.assertEqual(ref, l)

    def test_FADE_scalarweighted(self):
        ########################################################################
        # Test Case 1
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13,
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 1, 'scalar_score': 0.55555}
        l = energy_functions.FADE_scalarweighted(contact)
        ref = "AtomPair CB 13 CB 50 SCALARWEIGHTEDFUNC 0.556 FADE -10 19 10 -15.00 0"
        self.assertEqual(ref, l)

    def test_SIGMOID_scalarweighted(self):
        ########################################################################
        # Test Case 1
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13,
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 1, 'scalar_score': 0.55555}
        l = energy_functions.SIGMOID_scalarweighted(contact)
        ref = "AtomPair CB 13 CB 50 SCALARWEIGHTEDFUNC 0.556 SIGMOID 8.00 1.00"
        self.assertEqual(ref, l)

    def test_SIGMOID_bbcontacts(self):
        ########################################################################
        # Test Case 1
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13,
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 1, 'scalar_score': 0.55555}
        l = energy_functions.SIGMOID_bbcontacts(contact)
        ref = "AtomPair CB 13 CB 50 SCALARWEIGHTEDFUNC 0.556 SUMFUNC 2 SIGMOID 10.250 1.560 CONSTANTFUNC -0.5"
        self.assertEqual(ref, l)
        
        ########################################################################
        # Test Case 2
        contact = {'atom1': 'CB', 'atom2': 'CB', 'res1': 'R', 'res2': 'M', 'res1_index': 13,
                   'res2_index': 50, 'raw_score': 0.333, 'weight': 2, 'scalar_score': 0.55555}
        l = energy_functions.SIGMOID_bbcontacts(contact)
        ref = "AtomPair CB 13 CB 50 SCALARWEIGHTEDFUNC 1.111 SUMFUNC 2 SIGMOID 10.250 1.560 CONSTANTFUNC -0.5"
        self.assertEqual(ref, l)

if __name__ == "__main__":
    unittest.main()
