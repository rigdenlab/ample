#!/usr/bin/env python

import os
import unittest

class PhaserPdbParser(object):
    """
    Class to mine information from a phaser pdb file
    """

    def __init__(self,pdbfile):

        self.pdbfile = pdbfile
        self.LLG = None
        self.TFZ = None

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
                        self.TFZ = float(i.replace("TFZ==", ""))
                        break
                    if "TFZ=" in i and "TFZ==" not in i and "*" not in i:
                        self.TFZ = float(i.replace("TFZ=", ""))
                        break
        
                for i in llist:
                    if "LLG==" in i:
                        self.LLG = float(i.replace("LLG==", ""))
                        break
                    if "LLG=" in i and "LLG==" not in i:
                        self.LLG = float(i.replace("LLG=", ""))
                        break
        return
    

class PhaserLogParser(object):
    """
    Class to mine information from a phaser log
    """

    def __init__(self,logfile, noLLG=True):


        self.noLLG = noLLG # Whether to skip the LLG/TFZ parsing if we are getting these from the PDB
        self.logfile = logfile
        self.LLG = None
        self.TFZ = None
        self.time = None
        self.killed = False

        self.parse()

        return
    
    
    def parse(self):
        """parse"""

        # print os.path.join(os.getcwd(), logfile)
        #print "Checking logfile ",self.logfile
        
        # First read the file backwards to get the time and see if it was killed
        BUFFSIZE=1024
        with open(self.logfile, 'r') as fh:
            
            # Find out how long the file is
            fh.seek( 0, os.SEEK_END )
            fsize = fh.tell()
            
            # Position pointer BUFFSIZE bytes from end
            fh.seek( max( fsize - BUFFSIZE, 0 ), os.SEEK_SET )
            
            # Now read lines in reverse
            for line in reversed( fh.readlines() ) :
                if "CPU Time" in line:
                    self.time=float( line.split()[-2] )
                    continue
                
                if line.startswith("KILL-TIME ELAPSED ERROR: Job killed for elapsed time exceeding limit"):
                    self.killed=True
                    return
        
        # Return if we're reading the LLG from the PDB
        if self.noLLG:
            return

        fh = open(self.logfile, 'r')
        CAPTURE = False
        solline = ""
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
                self.TFZ = float(i.replace("TFZ==", ""))
                break
            if "TFZ=" in i and "TFZ==" not in i and "*" not in i:
                self.TFZ = float(i.replace("TFZ=", ""))
                break

        for i in llist:
            if "LLG==" in i:
                self.LLG = float(i.replace("LLG==", ""))
                break
            if "LLG=" in i and "LLG==" not in i:
                self.LLG = float(i.replace("LLG=", ""))
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
        
        self.assertEqual( pp.LLG, 11)
        self.assertEqual( pp.TFZ, 6.2)
        
        return
    
    def testPdbParser2(self):
        """foo"""
        
        pdb = "/media/data/shared/coiled-coils/1BYZ/ROSETTA_MR_0/MRBUMP/cluster_1/search_All_atom_trunc_0.039428_rad_1_mrbump/data/loc0_ALL_All_atom_trunc_0.039428_rad_1/unmod/mr/phaser/refine/phaser_loc0_ALL_All_atom_trunc_0.039428_rad_1_UNMOD.1.pdb"
        
        pp = PhaserPdbParser( pdb )
        
        self.assertEqual( pp.LLG, 36)
        self.assertEqual( pp.TFZ, 8.6)
        
        return
    

    def testLogParser1(self):
        """foo"""
        
        log = "/media/data/shared/coiled-coils/1BYZ/ROSETTA_MR_0/MRBUMP/cluster_1/search_SCWRL_reliable_sidechains_trunc_0.005734_rad_1_mrbump/data/loc0_ALL_SCWRL_reliable_sidechains_trunc_0.005734_rad_1/unmod/mr/phaser/phaser_loc0_ALL_SCWRL_reliable_sidechains_trunc_0.005734_rad_1_UNMOD.log"
        
        pp = PhaserLogParser( log )
        
        self.assertEqual( pp.time, 24991.56)
        self.assertEqual( pp.LLG, None)
        self.assertEqual( pp.TFZ, None)
        self.assertEqual( pp.killed, True)
        
        return
    
    def testLogParser2(self):
        """foo"""
        
        log = "/media/data/shared/coiled-coils/1BYZ/ROSETTA_MR_0/MRBUMP/cluster_1/search_All_atom_trunc_0.039428_rad_1_mrbump/data/loc0_ALL_All_atom_trunc_0.039428_rad_1/unmod/mr/phaser/phaser_loc0_ALL_All_atom_trunc_0.039428_rad_1_UNMOD.log"
        
        pp = PhaserLogParser( log, noLLG=False )
        self.assertEqual( pp.time, 9648.5)
        self.assertEqual( pp.LLG, 36)
        self.assertEqual( pp.TFZ, 8.6)
        self.assertEqual( pp.killed, False)
        return
#
# Run unit tests
if __name__ == "__main__":
    unittest.main()
