"""Test functions for util.tmscore_util"""

import unittest
from ample.util import tmscore_util

class TestTMScore(unittest.TestCase):

    def test_gaps(self):        
        TM = tmscore_util.TMscore("cho")
        
        seq = "AAAA---AA--AA"
        gaps = TM.find_gaps(seq)
        ref_gaps = [5, 6, 7, 10, 11]
        self.assertEqual(ref_gaps, gaps)

        seq = "---AA-AA"
        gaps = TM.find_gaps(seq)
        ref_gaps = [1, 2, 3, 6]
        self.assertEqual(ref_gaps, gaps)

        seq = "-AAA--"
        gaps = TM.find_gaps(seq)
        ref_gaps = [1, 5, 6]
        self.assertEqual(ref_gaps, gaps)

if __name__ == "__main__":
    unittest.main()
