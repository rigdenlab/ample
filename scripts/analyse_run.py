#!/usr/bin/env python
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

sys.path.append("/Users/jmht/Documents/AMPLE/ample-dev1/python")
import mrbump_results
import pdbEd


# For now assume we just use the first model/chain from the native pdb

# We need to cout models/chains and then work out what to do accordingly 

def reforigin_rmsd( refinedPdb, nativePdbStd, model=None, chain=None ):
    """Use reforigin to calculate rmsd between native and refined"""
    
    workdir=os.getcwd()
    
    # Now create a PDB with the matching atoms from native that are in refined
    n = os.path.splitext( os.path.basename( nativePdbStd ) )[0]
    nativePdbStdMatch = os.path.join( workdir, n+"_matched.pdb" )
    pdbEd.keep_matching( refpdb=refinedPdb, targetpdb=nativePdbStd, outpdb=nativePdbStdMatch )
    
    # Now get the rmsd
    n = os.path.splitext( os.path.basename( refinedPdb ) )[0]
    reforiginOut = os.path.join( workdir, n+"_model{0}_chain{1}_reforigin.pdb".format( model,chain ) )
    return pdbEd.reforigin_rmsd( refpdb=nativePdbStdMatch, targetpdb=refinedPdb, outpdb=reforiginOut )


result = mrbump_results.MrBumpResult()


# workdir = "/home/jmht/Documents/test/3PCV"
# result.resultDir = "/media/data/shared/TM/3PCV/ROSETTA_MR_0/MRBUMP/cluster_1/search_poly_ala_trunc_2.822761_rad_1_phaser_mrbump/data/loc0_ALL_poly_ala_trunc_2.822761_rad_1/unmod/mr/phaser"
# nativePdb = "/media/data/shared/TM/3PCV/3PCV.pdb"
# fasta = "/media/data/shared/TM/3PCV/3PCV.fasta"
# result.program = "phaser"
# result.ensembleName = "poly_ala_trunc_2.822761_rad_1"

workdir = "/home/jmht/Documents/test/1GU8"
result.resultDir = "/media/data/shared/TM/1GU8/ROSETTA_MR_0/MRBUMP/cluster_1/search_SCWRL_reliable_sidechains_trunc_19.511671_rad_3_phaser_mrbump/data/loc0_ALL_SCWRL_reliable_sidechains_trunc_19.511671_rad_3/unmod/mr/phaser"
nativePdb = "/media/data/shared/TM/1GU8/1GU8.pdb"
fasta = "/media/data/shared/TM/1GU8/1GU8.fasta"
result.program = "phaser"
result.ensembleName = "SCWRL_reliable_sidechains_trunc_19.511671_rad_3"

# workdir = "/home/jmht/Documents/test/3U2F"
# result.resultDir = "/media/data/shared/TM/3U2F/ROSETTA_MR_0/MRBUMP/cluster_1/search_poly_ala_trunc_0.21093_rad_2_molrep_mrbump/data/loc0_ALL_poly_ala_trunc_0.21093_rad_2/unmod/mr/molrep"
# nativePdb = "/media/data/shared/TM/3U2F/3U2F.pdb"
# fasta = "/media/data/shared/TM/3U2F/3U2F.fasta"
# result.program = "molrep"
# result.ensembleName = "poly_ala_trunc_0.21093_rad_2"

os.chdir( workdir )

refinedPdb = os.path.join( result.resultDir, "refine", "refmac_{0}_loc0_ALL_{1}_UNMOD.pdb".format( result.program, result.ensembleName ) )
refinedInfo = pdbEd.get_info( refinedPdb )

# Run a pass to find the # chains/models and extract solvent content etc
nativeInfo = pdbEd.get_info( nativePdb )

#if len(nativeInfo.models) > 1:
#    raise RuntimeError,"More then one model!"

multiModel=False
if len( nativeInfo.models ) > 1:
    multiModel=True
    
for imodel, model in enumerate( nativeInfo.models ):
    
    if multiModel:
        modelNum = model.serial
        # Need to extract the relevant model
        n = os.path.splitext( os.path.basename( nativePdb ) )[0]
        nativePdbModel = os.path.join( workdir, n+"_model{0}.pdb".format( modelNum ) )
        pdbEd.extract_model( nativePdb, nativePdbModel, modelID=modelNum)
    else:
        modelNum = 0
        nativePdbModel = nativePdb
        
    print "processing model ",modelNum
    
    # First standardise the native pdb
    n = os.path.splitext( os.path.basename( nativePdbModel ) )[0]
    nativePdbStd = os.path.join( workdir, n+"_std.pdb" )
    pdbEd.standardise( nativePdbModel, nativePdbStd )
    
    # Make sure the number of chains in the native and refined match
    native_chains = nativeInfo.models[ imodel ].chains
    refined_chains = refinedInfo.models[ 0 ].chains # only ever one model in the refined pdb
    
    if len( native_chains ) == len( refined_chains ):
        
        print "MATCHING CHAINS "
        rmsd = reforigin_rmsd( refinedPdb, nativePdbStd, model=modelNum, chain='A' )
        print rmsd
        
    else:
        
        rmsds = {} # dict of rmsd -> ( chainIDnative, chainIDrefined )
        print "DIFFERENT NUMBERS OF CHAINS ", native_chains, refined_chains
        
        # Match each chain in native against refined and pick the best
        for nativeChainID in native_chains:
            
            #print "native_chain: {0}".format( nativeChainID )
                    
            if len( native_chains ) == 1:
                # Don't need to do owt as we are just using the native as is
                nativeChainPdb = nativePdbStd
            else:
                
                # Extract the chain from the pdb
                n = os.path.splitext( os.path.basename( nativePdbStd ) )[0]
                nativeChainPdb = os.path.join( workdir, n+"_chain{0}.pdb".format( nativeChainID ) ) 
                pdbEd.extract_chain( nativePdbStd, nativeChainPdb, chainID=nativeChainID )
                
                assert os.path.isfile( nativeChainPdb  ), nativeChainPdb
            
            for refinedChainID in refined_chains:
                
                #print "refined_chain: {0}".format( refinedChainID )
                
                assert os.path.isfile( nativeChainPdb  ), nativeChainPdb
                
                # Extract the chain from the pdb
                n = os.path.splitext( os.path.basename( refinedPdb ) )[0]
                refinedChainPdb = os.path.join( workdir, n+"_chain{0}.pdb".format( refinedChainID ) ) 
                pdbEd.extract_chain( refinedPdb, refinedChainPdb, chainID=refinedChainID, newChainID=nativeChainID )
                
                #print "calculating for {0} vs. {1}.".format( refinedChainPdb, nativeChainPdb  )
                
                rmsd = reforigin_rmsd( refinedChainPdb, nativeChainPdb, model=modelNum, chain=refinedChainID )
                
                rmsds[ rmsd ] = ( nativeChainID, refinedChainID )
                
        # End loop over chains
        # Now pick the best...
        print rmsds
        pass
    
    # End loop over models
    

sys.exit()



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