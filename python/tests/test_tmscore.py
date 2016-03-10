"""Test functions for python.tmscore"""

import unittest

from ample.python.tmscore import TMscorer

class TestTMScore(unittest.TestCase):
    def setUp(self):
        self.TM = TMscorer("foo", "bar", "cho")

    def testRead(self):
        seq1 = "-----AAAA---"
        ref1_offset = 5

        seq2 = "AAA---"
        ref2_offset = 0

        seq3 = "AAAA"
        ref3_offset = 0

        offset1 = self.TM.read_sequence(seq1)
        offset2 = self.TM.read_sequence(seq2)
        offset3 = self.TM.read_sequence(seq3)

        self.assertEqual(ref1_offset, offset1)
        self.assertEqual(ref3_offset, offset2)
        self.assertEqual(ref3_offset, offset3)

    def testGaps(self):
        seq1 = "AAAA---AA--AA"
        ref_gaps1 = [5, 6, 7, 10, 11]

        seq2 = "---AA-AA"
        ref_gaps2 = [1, 2, 3, 6]

        seq3 = "-AAA--"
        ref_gaps3 = [1, 5, 6]

        gaps1 = self.TM.find_gaps(seq1)
        gaps2 = self.TM.find_gaps(seq2)
        gaps3 = self.TM.find_gaps(seq3)

        self.assertEqual(ref_gaps1, gaps1)
        self.assertEqual(ref_gaps2, gaps2)
        self.assertEqual(ref_gaps3, gaps3)

if __name__ == "__main__":
    unittest.main()