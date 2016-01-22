#!/usr/bin/env ccp4-python


import collections
import os
import unittest


class PsipredSs2Parser(object):
    """Parser for psipred ss2 file"""
    
    def __init__(self, ss2file=None):
        self.residues = None
        if ss2file: self.parse(ss2file)
        
    def parse(self, ss2file):
        
        PSIPredResidueInfo = collections.namedtuple("PSIPredResidueInfo", ["rank", "residue", "ss", "coil", "helix", "strand"])
        residues = []
        
        with open(ss2file, 'r') as fhin:
            for line in iter( fhin.readline, '' ):
                if line[0] == '#' or not line.strip(): continue
                
                line = line.split()
                
                rank = int(line[0])
                residue = line[1]
                ss = line[2]
                coil, helix, strand = map(float, line[3:6])
                
                residues.append(PSIPredResidueInfo(rank=rank, residue=residue, ss=ss, 
                                                   coil=coil, helix=helix, strand=strand))
        
        self.residues = tuple(residues)
        return
        
    def checkContent(self):
        # Calculate secondary structure content in sequence          
        H = len([i for i in self.residues if i.ss=="H"])
        E = len([i for i in self.residues if i.ss=="E"])

        #E_content = float(E) / len(pred) * 100
        #H_content = float(H) / len(pred) * 100
        if H > 0 and E > 0:
            print  'Your protein is predicted to be mixed alpha beta, your chances of success are intermediate'
        if H == 0 and E > 0:
            print  'Your protein is predicted to be all beta, your chances of success are low'
        if H > 0 and E == 0:
            print  'Your protein is predicted to be all alpha, your chances of success are high'
        if  H == 0 and E == 0:
            print  'Your protein is has no predicted secondary structure, your chances of success are low'
        return
    
    def getSecondaryStructure(self):
        return "".join([i.ss for i in self.residues])
##End PsipredSs2Parser()


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
        ss2file = os.path.join(self.testfiles_dir, "1aba_.psipred_ss2")
        PA = PsipredSs2Parser(ss2file)
        
        ref_ss2 = "CEEEEEECCCCCCCCHHHHHHHHHHHCCCCEEEEEECCCCCCCCHHHHHHHHHHHCCCCCCCCCCCEEEEECCEEEECHHHHHHHHC"
        ss2 = "".join([i.ss for i in PA.residues])
        self.assertEqual(ref_ss2, ss2 )
        

