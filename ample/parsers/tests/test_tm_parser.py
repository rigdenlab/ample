"""Test functions for parsers.tm_parser"""
import os
import unittest
from ample import constants
from ample.parsers import tm_parser


class TestTMscore(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share, 'testfiles')

    def test_parse(self):
        logfile = os.path.join(self.testfiles_dir, "tmscore.log")
        TM = tm_parser.TMscoreLogParser()
        TM.parse(logfile)

        self.assertEqual(173, TM.nr_residues_common)
        self.assertEqual(6.654, TM.rmsd)
        self.assertEqual(0.5512, TM.tm)
        self.assertEqual(0.3147, TM.maxsub)
        self.assertEqual(0.4292, TM.gdtts)
        self.assertEqual(0.2283, TM.gdtha)


class TestTMalign(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share, 'testfiles')

    def test_parse(self):
        logfile = os.path.join(self.testfiles_dir, "tmalign.log")
        TM = tm_parser.TMalignLogParser()
        TM.parse(logfile)

        self.assertEqual(143, TM.nr_residues_common)
        self.assertEqual(0.70502, TM.tm)
        self.assertEqual(2.68, TM.rmsd)
        self.assertEqual(0.182, TM.seq_id)


if __name__ == "__main__":
    unittest.main()
