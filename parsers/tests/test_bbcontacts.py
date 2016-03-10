"""Test functions for parsers.bbcontacts"""

import unittest

from ample.parsers import bbcontacts

class Test(unittest.TestCase):
    def setUp(self):
        self.bp = bbcontacts.BBcontactsContactParser()

    def testFirstLast(self):
        self.bp.contacts = [{'res2': 43, 
                             'res1': 37, 
                             'raw_score': -4.767827, 
                             'strand_index': 7, 
                             'file': 'test', 
                             'strand_orientation': 'Antiparallel', 
                             'internal_strand_position': 'first', 
                             'diversity_factor': 1.02, 
                             'identifier': '1ABA', 
                             'method': 'bbcontacts'}]
        
        line_last1 = "1ABA 1.02 Antiparallel -4.767827 7 last 44 36".split()
        line_last2 = "1ABA 1.02 Antiparallel -4.767827 7 internal 44 36".split()
        
        self.assertTrue(self.bp.isFirstLast(line_last1))
        self.assertFalse(self.bp.isFirstLast(line_last2))

if __name__ == "__main__":
    unittest.main()
