"""Test functions for python.tmscore"""

import unittest

from ample.python import tmscore

class TestTMScore(unittest.TestCase):

    def test_read(self):
        TM = tmscore.TMscorer("foo", "bar", "cho")
        
        seq = "-----AAAA---"
        offset = TM.read_sequence(seq)
        ref_offset = 5
        self.assertEqual(ref_offset, offset)

        seq = "AAA---"
        offset = TM.read_sequence(seq)
        ref_offset = 0
        self.assertEqual(ref_offset, offset)

        seq = "AAAA"
        offset = TM.read_sequence(seq)
        ref_offset = 0
        self.assertEqual(ref_offset, offset)

    def test_gaps(self):        
        TM = tmscore.TMscorer("foo", "bar", "cho")
        
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