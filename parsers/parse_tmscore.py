#!/usr/bin/env ccp4-python

import os
import unittest

class TMscoreLogParser(object):
    """ Class to mine information from a tmscore log """

    def __init__(self):
        self.tm=0.0
        self.maxsub=0.0
        self.gdtts=0.0
        self.gdtha=0.0
        self.rmsd=0.0
        self.nrResiduesCommon=0
        return

    def parse(self, logfile):

        with open(logfile, 'r') as f:
            for line in iter(f.readline, ''):
                if line.startswith("Number of residues in common"):
                    self.nrResiduesCommon=int(line.split()[5])
                if line.startswith("RMSD"):         self.rmsd  = float(line.split()[5])
                if line.startswith("TM-score"):     self.tm    = float(line.split()[2])
                if line.startswith("MaxSub-score"): self.maxsub= float(line.split()[1])
                if line.startswith("GDT-TS-score"): self.gdtts = float(line.split()[1])
                if line.startswith("GDT-HA-score"): self.gdtha = float(line.split()[1])

        return
##End TMscoreLogParser

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -1 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        return


    def testParse(self):
        logfile = os.path.join(self.testfiles_dir, "tmscore.log")
        TM = TMscoreLogParser()
        TM.parse(logfile)

        self.assertEqual(173, TM.nrResiduesCommon)
        self.assertEqual(6.654, TM.rmsd)
        self.assertEqual(0.5512, TM.tm)
        self.assertEqual(0.3147, TM.maxsub)
        self.assertEqual(0.4292, TM.gdtts)
        self.assertEqual(0.2283, TM.gdtha)
##End Test
