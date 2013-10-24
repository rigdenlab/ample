#!/usr/bin/env python

import unittest

class PhaserPdbParser(object):
    """
    Class to mine information from a phaser pdb file
    """

    def __init__(self,pdbfile):

        self.pdbfile = pdbfile
        self.phaserLLG = None
        self.phaserTFZ = None

        self.parse()

        return

    def parse(self):
        """parse"""

        # print os.path.join(os.getcwd(), logfile)

        for line in open(self.pdbfile, 'r'):
            if line.startswith("REMARK") and "TFZ=" in line:
                llist = line.split()
                llist.reverse()
                for i in llist:
                    if "TFZ==" in i and "*" not in i:
                        self.phaserTFZ = float(i.replace("TFZ==", ""))
                        break
                    if "TFZ=" in i and "TFZ==" not in i and "*" not in i:
                        self.phaserTFZ = float(i.replace("TFZ=", ""))
                        break
        
                for i in llist:
                    if "LLG==" in i:
                        self.phaserLLG = float(i.replace("LLG==", ""))
                        break
                    if "LLG=" in i and "LLG==" not in i:
                        self.phaserLLG = float(i.replace("LLG=", ""))
                        break
        return
    

class PhaserLogParser(object):
    """
    Class to mine information from a phaser log
    """

    def __init__(self,logfile):

        self.logfile = logfile
        self.phaserLLG = None
        self.phaserTFZ = None
        self.phaserTime = None

        self.parse()

        return

    def parse(self):
        """parse"""

        # print os.path.join(os.getcwd(), logfile)
        fh = open(self.logfile, 'r')
        
        #print "Checking logfile ",self.logfile

        for line in reversed(fh.readlines()):
            if "CPU Time" in line:
                self.phaserTime=float( line.split()[-2] )
                break
        fh.close()

        fh = open(self.logfile, 'r')
        CAPTURE = False
        solline = ""
        self.phaserLLG = 0.0
        self.phaserTFZ = 0.0
        line = fh.readline()
        while line:
            if CAPTURE:
                if "SOLU SPAC" in line:
                    CAPTURE = False
                else:
                    solline += line.strip() + " "
            if  "Solution #1 annotation (history):" in line or  "Solution annotation (history):" in line:
                CAPTURE = True
            line = fh.readline()
        fh.close()

        llist = solline.split()
        llist.reverse()
        for i in llist:
            if "TFZ==" in i and "*" not in i:
                self.phaserTFZ = float(i.replace("TFZ==", ""))
                break
            if "TFZ=" in i and "TFZ==" not in i and "*" not in i:
                self.phaserTFZ = float(i.replace("TFZ=", ""))
                break

        for i in llist:
            if "LLG==" in i:
                self.phaserLLG = float(i.replace("LLG==", ""))
                break
            if "LLG=" in i and "LLG==" not in i:
                self.phaserLLG = float(i.replace("LLG=", ""))
                break
        return
    
class TestParsers(unittest.TestCase):
    """
    Unit test
    """
        
    def testPdbParser1(self):
        """foo"""
        
        pdb = "/media/data/shared/coiled-coils/1BYZ/ROSETTA_MR_0/MRBUMP/cluster_1/search_SCWRL_reliable_sidechains_trunc_0.005734_rad_1_mrbump/data/loc0_ALL_SCWRL_reliable_sidechains_trunc_0.005734_rad_1/unmod/mr/phaser/refine/phaser_loc0_ALL_SCWRL_reliable_sidechains_trunc_0.005734_rad_1_UNMOD.1.pdb"
        
        pp = PhaserPdbParser( pdb )
        
        self.assertEqual( pp.phaserLLG, 11)
        self.assertEqual( pp.phaserTFZ, 6.2)
        
        return
    
    def testPdbParser2(self):
        """foo"""
        
        pdb = "/media/data/shared/coiled-coils/1BYZ/ROSETTA_MR_0/MRBUMP/cluster_1/search_All_atom_trunc_0.039428_rad_1_mrbump/data/loc0_ALL_All_atom_trunc_0.039428_rad_1/unmod/mr/phaser/refine/phaser_loc0_ALL_All_atom_trunc_0.039428_rad_1_UNMOD.1.pdb"
        
        pp = PhaserPdbParser( pdb )
        
        self.assertEqual( pp.phaserLLG, 36)
        self.assertEqual( pp.phaserTFZ, 8.6)
        
        return
#
# Run unit tests
if __name__ == "__main__":
    unittest.main()
