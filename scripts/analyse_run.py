'''
Created on 19 Jul 2013

@author: jmht

Data we will be collecting:

Target:
length in AA
resolution in A
secondary structure (% and per AA)
?radius of gyration

Rosetta Models:
Score
maxsub cf target

Cluster:

Ensemble:
number of models (and their scores)
truncation level
number residues
side_chain_treatment

Solution
number (and identity? of residues)

shelxe rebuild:
* CC
* av. fragment length
* RMSD to native
* Maxsub to native
* TM to native

remac refined result
* reforigin score to native
* rmsd to native
* maxsub to native
* TM to native

TO THINK ABOUT
* multiple models/chains in native
* multiple chains in solution (e.g. 3PCV)



'''
import os
import re
import sys
import unittest

import mrbump_results
import pdbEd


nativePdb = "/media/data/shared/TM/3U2F/3U2F.pdb"
fasta = "/media/data/shared/TM/3U2F/3U2F.fasta"
workdir = "/home/jmht/Documents/test/3U2F/test"
# For now assume we just use the first model/chain from the native pdb


result = mrbump_results.MrBumpResult()
result.resultDir = "/media/data/shared/TM/3U2F/ROSETTA_MR_0/MRBUMP/cluster_1/search_poly_ala_trunc_0.21093_rad_2_molrep_mrbump/data/loc0_ALL_poly_ala_trunc_0.21093_rad_2/unmod/mr/molrep"
result.program = "molrep"
resultk = "loc0_ALL_poly_ala_trunc_0.21093_rad_2_UNMOD"
result.ensembleName = "poly_ala_trunc_0.21093_rad_2"

os.chdir( workdir )

# Run a pass to find the # chains/models and extract solvent content etc
info = pdbEd.get_info( nativePdb )

# Hard-wire in first native chain for now
nativeChainID = info.models[0].chains[0]

# Extract first model/chain from native
d, f = os.path.split( nativePdb )
name, ext = os.path.splitext( f )
nativePdb_m1c1 = os.path.join( workdir, name+"_m1c1"+ext )
#retcode = pdbEd.to_1_std_chain( nativePdb, nativePdb_m1c1 )


# Get path of refined structure
refinedPdb = os.path.join( result.resultDir, "refine", "refmac_{0}_loc0_ALL_{1}_UNMOD.pdb".format( result.program, result.ensembleName ) )

# Find out how many chains it contains.
info = pdbEd.get_info( refinedPdb )

# Loop over each chain and run a comparison
for chainID in info.models[0].chains:
    
    # Extract the chain from the pdb
    n = os.path.splitext( os.path.basename( refinedPdb ) )[0]
    chainPdb = os.path.join( workdir, n+"_chain{0}.pdb".format(chainID) )
    pdbEd.extract_chain( refinedPdb, chainPdb, chainID=chainID, newChainID=nativeChainID )
    
    # Extract the matching atoms from the nativePdb
    n = os.path.splitext( os.path.basename( nativePdb_m1c1 ) )[0]
    matchingPdb = os.path.join( workdir, n+"_matching_chain{0}.pdb".format(chainID) )
    pdbEd.keep_matching( chainPdb, nativePdb_m1c1, matchingPdb )
    
    n = os.path.splitext( os.path.basename( chainPdb ) )[0]
    outpdb = os.path.join( workdir, n+"_reforigin.pdb" )    
    rms = pdbEd.reforigin_rmsd(matchingPdb, chainPdb, outpdb )
    print rms

sys.exit()


# Extract matching atoms from natvePdb
n = os.path.splitext( os.path.basename( nativePdb_m1c1 ) )[0]
matchingPdb = os.path.join( workdir, n+"_matching.pdb")
pdbEd.keep_matching(refinedPdb, nativePdb_m1c1, matchingPdb)

# CHECK TWO CHAINS!

print nativePdb_m1c1
print refinedPdb
print matchingPdb



class ShelxeLogParser(object):
    """
    Class to mine information from a shelxe log
    """
    
    def __init__(self,logfile):
        
        self.logfile = logfile
        self.CC = None
        self.avgChainLength = None
        
        self.parse()
        
        return
        
    def parse(self):
        """Parse a shelxe log file to get the CC and average Chain length
        """
        
        cycleData = [] # List (CC,avgChainLength) tuples - ordered by cycle
        fh = open( self.logfile, 'r')
        
        line = fh.readline()
        while line:
            
            # find should be quicker then re match
            if line.find("residues left after pruning, divided into chains as follows:") != -1:
                (cc, avgChainLength) = self._parseCycle(fh)
                cycleData.append( (cc, avgChainLength) )
            
            
            if  line.find( "Best trace (cycle" ) != -1:
                # Expecting:
                #  "Best trace (cycle   1 with CC 37.26%) was saved as shelxe-input.pdb"
                cycle = int( re.search("\s\d+\s", line).group(0) )
                cc = float( re.search("\s\d+\.\d+", line).group(0) )
                
                # Check it matches
                if cycleData[ cycle-1 ][0] != cc:
                    raise RuntimeError,"Error getting final CC!"
                
                self.CC =  cycleData[ cycle-1 ][0]
                self.avgChainLength = cycleData[ cycle-1 ][1]

            line = fh.readline()
        #End while
        
        fh.close()
        
        return
        
    def _parseCycle(self, fh):
        """
        Working on assumption each cycle contains something like the below:
<log>
           223 residues left after pruning, divided into chains as follows:
 A:   6   B:   7   C:   6   D:   6   E:   8   F:   7   G:  12   H:  12   I:   5
 J:  10   K:   6   L:   6   M:   6   N:   7   O:   6   P:   7   Q:   8   R:   6
 S:   5   T:   6   U:  10   V:   9   W:  12   X:  11   Y:   8   Z:   6   Z:   6
 Z:   6   Z:   7   Z:   6

 CC for partial structure against native data =  30.92 %
 </log>
 """
        
        lengths = []
        while True:
            
            line = fh.readline().strip()
            line = line.rstrip(os.linesep)
            if not line:
                # End of reading lengths
                break
            
            # Loop through integers & add to list
            for m in re.finditer("\s\d+", line):
                lengths.append( int(m.group(0)) )
                
        # Now calculate average chain length
        if not len( lengths ):
            raise RuntimeError, "Failed to read any fragment lengths"
        
        # Average chain lengths
        avgChainLength = sum(lengths) / int( len(lengths) )        
        
        # Here should have read the  lengths so now just get the CC
        count=0
        while True:
            line = fh.readline().strip()
            if line.startswith("CC for partial structure against native data"):
                break
            else:
                count += 1
                if count > 5:
                    raise RuntimeError,"Error parsing CC score"
            
        cc = float( re.search("\d+\.\d+", line).group(0) )
        
        return ( cc, avgChainLength )


class Test(unittest.TestCase):


    def testShelxeLogParser(self):
        logfile = "/media/data/shared/TM/2BHW/ROSETTA_MR_0/MRBUMP/cluster_1/search_poly_ala_trunc_9.355791_rad_3_molrep_mrbump/" + \
        "data/loc0_ALL_poly_ala_trunc_9.355791_rad_3/unmod/mr/molrep/build/shelxe/shelxe_run.log"
        
        p = ShelxeLogParser( logfile )
        self.assertEqual(37.26, p.CC)
        self.assertEqual(7, p.avgChainLength)
        

#if __name__ == "__main__":
#    #import sys;sys.argv = ['', 'Test.testName']
#    unittest.main()

