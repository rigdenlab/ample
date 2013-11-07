import os
import shutil

import ample_util
import pdb_edit
import residue_map

class ReforiginRmsd(object):
    """Class to use reforigin to determine how well the model was placed.
    """
    
    def __init__( self, nativePdb, refinedPdb, refModelPdb, workdir=None, cAlphaOnly=True ):
        
        self.workdir = workdir
        if not self.workdir:
            self.workdir = os.getcwd()
        
        self.cAlphaOnly = cAlphaOnly# Whether to only compare c-alpha atoms
        
        self.rmsd = None
        self.bestNativeChain = None
        self.bestRefinedChain = None
        self.bestReforiginPdb = None
        self.refModelPdb = None
        
        self.run( nativePdb, refinedPdb, refModelPdb )
        
    def calc_reforigin_rmsd( self, refpdb=None, targetpdb=None, outpdb=None, DMAX=100 ):
        
        assert refpdb and targetpdb and outpdb
        
        # HACK - REFORIGIN has a limit on the length of the command line, so we need to create copies of inputfile
        # as this has the potentially longest path
        tmptarget = ample_util.tmpFileName()+".pdb"
        shutil.copy(targetpdb, tmptarget)
        
        logfile = outpdb+".log"
        cmd="reforigin xyzin {0} xyzref {1} xyzout {2} DMAX {3}".format( tmptarget, refpdb, outpdb, DMAX ).split()
        
        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False)
        
        if retcode != 0:
            raise RuntimeError, "Error running command: {0}".format( " ".join(cmd) )
        else:
            os.unlink( tmptarget )
        
        # Parse logfile to get RMSD
        rms = None
        for line in open( logfile, 'r' ):
            if line.startswith("RMS deviation:"):
                rms = float( line.split()[-1] )
                break
        
        if not rms:
            raise RuntimeError, "Error extracting RMS from logfile: {0}".format( logfile )
        
        return rms

    def reforigin_rmsd( self, refinedPdb, nativePdb, nativeChainID=None ):
        """Use reforigin to calculate rmsd between native and refined.
        
        NB: Still have a bug with (e.g. 2UUI) where the final residue in the native file only has
        an N for the last residue (ALA 150). If we calculate a map and then strip to C-alpha, there
        is a missing residue.
        """
        
        # Calculate the RefSeqMap - need to do this before we reduce to c-alphas
        resSeqMap = residue_map.residueSequenceMap( nativePdb, self.refModelPdb )
        
        # Find out if there are atoms in the model that we need to remove
        incomparable = resSeqMap.targetIncomparable()
        if len( incomparable ):
            
            refinedPdbCut = ample_util.filename_append( filename=refinedPdb, astr="cut", directory=self.workdir )
            logfile = "{0}.log".format( refinedPdb )
            cmd="pdbcur xyzin {0} xyzout {1}".format( refinedPdb, refinedPdbCut ).split()
            
            # Build up stdin - I'm too thick to work out the selection syntax for a discrete list
            stdin = ""
            for e in incomparable:
                stdin += "delresidue {0}\n".format( e )
            
            retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=self.workdir, dolog=False, stdin=stdin)
            
            if retcode == 0:
                # remove temporary files
                os.unlink(logfile)
            else:
                raise RuntimeError,"Error deleting residues {0}".format( incomparable )
            
            refinedPdb = refinedPdbCut
        
        PE = pdb_edit.PDBEdit()
        if self.cAlphaOnly:
            # If only alpha atoms are required, we create a copy of the model with only alpha atoms
            tmp = ample_util.filename_append( filename=refinedPdb, astr="CA", directory=self.workdir )
            PE.calpha_only( refinedPdb, tmp )
            refinedPdb = tmp
        else:
            # Strip down to backbone atoms
            tmp = ample_util.filename_append( filename=refinedPdb, astr="BB", directory=self.workdir )
            PE.backbone( refinedPdb, tmp  )
            refinedPdb = tmp

        # Now create a PDB with the matching atoms from native that are in refined
        nativePdbMatch = ample_util.filename_append( filename=nativePdb, astr="matched", directory=self.workdir )
        PE.keep_matching( refpdb=refinedPdb, targetpdb=nativePdb, outpdb=nativePdbMatch, resSeqMap=resSeqMap )
        
        # Now get the rmsd
        astr = "chain{0}_reforigin".format( nativeChainID )
        reforiginOut = ample_util.filename_append( filename=refinedPdb, astr=astr, directory=self.workdir )
        rms = self.calc_reforigin_rmsd( refpdb=nativePdbMatch, targetpdb=refinedPdb, outpdb=reforiginOut )
        return ( rms, reforiginOut )
    
    def run( self, nativePdb, refinedPdb, refModelPdb ):
        """For now just save lowest rmsd - can look at collecting more nativeInfo later
        
        Currently we assume we are only given one model and that it has already been standardised.
        """
        
        self.refModelPdb = refModelPdb

        # Run a pass to find the # chains
        pdbedit = pdb_edit.PDBEdit()
        refinedInfo = pdbedit.get_info( refinedPdb )
        nativeInfo = pdbedit.get_info( nativePdb )
        native_chains = nativeInfo.models[ 0 ].chains
        refined_chains = refinedInfo.models[ 0 ].chains # only ever one model in the refined pdb
        
        #print "got native chains ", native_chains
        #print "got refined chains ", refined_chains
            
        rmsds = {} # dict of rmsd -> ( chainIDnative, chainIDrefined, reforiginLogfile )
        
        # Match each chain in native against refined and pick the best
        for nativeChainID in native_chains:
            
            #print "native_chain: {0}".format( nativeChainID )
                    
            if len( native_chains ) == 1:
                # Don't need to do owt as we are just using the native as is
                nativeChainPdb = nativePdb
            else:
                
                # Extract the chain from the pdb
                astr = "chain{0}".format( nativeChainID )
                nativeChainPdb = ample_util.filename_append( filename=nativePdb, astr=astr, directory=self.workdir )
                pdbedit.extract_chain( nativePdb, nativeChainPdb, chainID=nativeChainID )
                
                assert os.path.isfile( nativeChainPdb  ), nativeChainPdb
            
            for refinedChainID in refined_chains:
                
                #print "refined_chain: {0}".format( refinedChainID )
                
                assert os.path.isfile( nativeChainPdb  ), nativeChainPdb
                
                # Extract the chain from the pdb
                astr = "chain{0}".format( refinedChainID )
                refinedChainPdb = ample_util.filename_append( filename=refinedPdb, astr=astr, directory=self.workdir )
                pdbedit.extract_chain( refinedPdb, refinedChainPdb, chainID=refinedChainID, newChainID=nativeChainID, cAlphaOnly=self.cAlphaOnly )
                
                #print "calculating for {0} vs. {1}".format( refinedChainID, nativeChainID  )
                #print "calculating for {0} vs. {1}".format( refinedChainPdb, nativeChainPdb  )
                try:
                    rmsd, refPdb  = self.reforigin_rmsd( refinedChainPdb, nativeChainPdb, nativeChainID=nativeChainID )
                except RuntimeError, e:
                    print "GOT REFORIGIN ERROR for {0},{1},{2}".format( refinedChainPdb, nativeChainPdb, nativeChainID )
                    print e
                    rmsd = 99999
                    refPdb = None
                
                rmsds[ rmsd ] = ( nativeChainID, refinedChainID, refPdb )
                
        # End loop over chains
        # Now pick the best...
        rmsd = sorted( rmsds.keys() )[ 0 ]
        #print "Got rmsds over chains: {0}".format( rmsds )
        
        self.rmsd = rmsd
        self.bestNativeChain = rmsds[ rmsd ][0]
        self.bestRefinedChain = rmsds[ rmsd ][1]
        self.bestReforiginPdb = rmsds[ rmsd ][2]
        #print "best chain rmsd is {0} for nativeChain {1} vs refinedChain {2}".format( self.rmsd, self.bestChains[0], self.bestChains[1] )
            
        return
