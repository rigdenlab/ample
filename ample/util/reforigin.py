import os

from ample.util import ample_util
from ample.util import pdb_edit
from ample.util import residue_map

class ReforiginRmsd(object):
    """Class to use reforigin to determine how well the model was placed.
    """
    
    def __init__( self, workdir=None, cAlphaOnly=True ):
        
        self.workdir = workdir
        if not self.workdir:
            self.workdir = os.getcwd()
        
        self.cAlphaOnly = cAlphaOnly# Whether to only compare c-alpha atoms
        
        self.rmsd = None
        self.bestNativeChain = None
        self.bestPlacedChain = None
        self.bestReforiginPdb = None
        
        return
    
    def calculate( self, refpdb=None, targetpdb=None, outpdb=None, DMAX=100 ):
        
        # Don't need to do this as we process the file before hand
        # HACK - REFORIGIN has a limit on the length of the command line, so we need to create copies of inputfile
        # as this has the potentially longest path
        #tmptarget = ample_util.tmp_file_name()+".pdb"
        #shutil.copy(targetpdb, tmptarget)
        
        logfile = outpdb+".log"
        cmd="reforigin xyzin {0} xyzref {1} xyzout {2} DMAX {3}".format( targetpdb, refpdb, outpdb, DMAX ).split()
        
        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False)
        if retcode != 0:
            raise RuntimeError, "Error running command: {0}".format( " ".join(cmd) )
        
        # Parse logfile to get RMSD
        rms = None
        for line in open( logfile, 'r' ):
            if line.startswith("RMS deviation:"):
                rms = float( line.split()[-1] )
                break
        
        if not rms:
            raise RuntimeError, "Error extracting RMS from logfile: {0}".format( logfile )
        else:
            os.unlink(logfile)
        
        return rms

    def preparePlacedPdb( self, placedPdb=None, placedChainID=None, nativeChainID=None, resSeqMap=None ):
        """
        Use pdbcur to:
        - extract chain to compare
        - strip down to CA/BB
        - remove any atoms that cannot be compared to the native
        """
        # Build up stdin
        
        # Extract the chain to compare
        stdin = "lvchain {0}\n".format( placedChainID )
        
        # Rename it to match the native
        if placedChainID != nativeChainID:
            stdin += "renchain {0} {1}\n".format( placedChainID, nativeChainID )
        
        # Find out if there are atoms in the model that we need to remove
        incomparable = resSeqMap.targetIncomparable( bbMask=not self.cAlphaOnly )
        if len( incomparable ):
            # Build up stdin - I'm too thick to work out the selection syntax for a discrete list
            for e in incomparable:
                stdin += "delresidue {0}\n".format( e )
                
        if self.cAlphaOnly:
            # Strip down to CA
            stdin += 'lvatom "CA[C]:*"\n'
        else:
            # Strip down to backbone atoms
            stdin += 'lvatom "N,CA,C,O,CB[N,C,O]"\n'
            
        # Renumber?
        stdin += "sernum\n"
        
        # Name the output file accordingly
        astr = "chain{0}".format( placedChainID )
        placedChainPdb = ample_util.filename_append( filename=placedPdb, astr=astr, directory=self.workdir )
        
        # Now run pdbcur to do it all
        cmd="pdbcur xyzin {0} xyzout {1}".format( placedPdb, placedChainPdb ).split()
        logfile = "{0}.log".format( placedChainPdb )
        retcode = ample_util.run_command( cmd=cmd, logfile=logfile, directory=self.workdir, dolog=False, stdin=stdin)
        if retcode != 0:
            raise RuntimeError,"Error extracting chain from placed PDB {0} in directory {1}".format( placedPdb, self.workdir )
        else:
            os.unlink(logfile)
        return placedChainPdb
    
    def getRmsd( self, nativePdbInfo=None, placedPdbInfo=None, refModelPdbInfo=None, workdir=None, cAlphaOnly=True  ):
        """For now just save lowest rmsd - can look at collecting more nativeInfo later
        
        Currently we assume we are only given one model and that it has already been standardised.
        """

        if workdir:
            self.workdir = workdir
        if not self.workdir:
            self.workdir = os.getcwd()
        
        self.cAlphaOnly = cAlphaOnly# Whether to only compare c-alpha atoms
        
        # Run a pass to find the # chains
        native_chains = nativePdbInfo.models[ 0 ].chains
        placed_chains = placedPdbInfo.models[ 0 ].chains
        
        #print "got native chains ", native_chains
        #print "got placed chains ", placed_chains
            
        rmsds = {} # dict of rmsd -> ( chainIDnative, chainIDrefined, reforiginLogfile )
        
        # Match each chain in native against refined and pick the best
        for nativeChainID in native_chains:
            
            #print "native_chain: {0}".format( nativeChainID )
                    
            if len( native_chains ) == 1:
                # Don't need to do owt as we are just using the native as is
                nativeChainPdb = nativePdbInfo.pdb
            else:
                # Extract the chain from the pdb
                astr = "chain{0}".format( nativeChainID )
                nativeChainPdb = ample_util.filename_append( filename=nativePdbInfo.pdb, astr=astr, directory=self.workdir )
                pdb_edit.extract_chain( nativePdbInfo.pdb, nativeChainPdb, chainID=nativeChainID )
                
            # Calculate the RefSeqMap - need to do this before we reduce to c-alphas
            # The second chain may be a different composition to the first, so we only generate a traceback if we fail
            # on the first chain. The model only has one chain, so the residueMap has to be the same for all the chains
            try:
                resSeqMap = residue_map.residueSequenceMap()
                resSeqMap.fromInfo( refInfo=nativePdbInfo,
                                    refChainID=nativeChainID,
                                    targetInfo=refModelPdbInfo,
                                    targetChainID='A' # Model only has one chain
                                    )
                
            except RuntimeError:
                if nativeChainID == native_chains[0]:
                    raise
                else:
                    # Only compare the first chain
                    break
                
            for placedChainID in placed_chains:
                
                # Prepare the placed PDB
                placedChainPdb = self.preparePlacedPdb( placedPdb=placedPdbInfo.pdb, placedChainID=placedChainID, nativeChainID=nativeChainID, resSeqMap=resSeqMap )
                
                # Now create a PDB with the matching atoms from native that are in refined
                nativePdbMatch = ample_util.filename_append( filename=nativeChainPdb, astr="matched", directory=self.workdir )
                pdb_edit.keep_matching( refpdb=placedChainPdb, targetpdb=nativeChainPdb, outpdb=nativePdbMatch, resSeqMap=resSeqMap )
                
                # Now get the rmsd
                astr = "chain{0}_reforigin".format( nativeChainID )
                reforiginOut = ample_util.filename_append( filename=placedChainPdb, astr=astr, directory=self.workdir )
                
                try:
                    rms = self.calculate( refpdb=nativePdbMatch, targetpdb=placedChainPdb, outpdb=reforiginOut )
                except RuntimeError, e:
                    print "GOT REFORIGIN ERROR for {0},{1},{2}".format( placedChainPdb, nativeChainPdb, nativeChainID )
                    print e
                    rms = 99999
                 
                rmsds[ rms ] = ( nativeChainID, placedChainID, reforiginOut )
                
                # Clean up
                os.unlink(placedChainPdb)
                os.unlink(nativePdbMatch)
        # End loop over chains
        
        # Clean up
        
        # Now pick the best...
        rmsd = sorted( rmsds.keys() )[ 0 ]
        #print "Got rmsds over chains: {0}".format( rmsds )
        
        self.rmsd = rmsd
        self.bestNativeChain = rmsds[ rmsd ][0]
        self.bestPlacedChain = rmsds[ rmsd ][1]
        self.bestReforiginPdb = rmsds[ rmsd ][2]
        
        # Clean up
        for k in rmsds.keys():
            if k != rmsd:
                try: os.unlink(rmsds[k][2])
                except: pass
        #print "best chain rmsd is {0} for nativeChain {1} vs refinedChain {2}".format( self.rmsd, self.bestChains[0], self.bestChains[1] )
            
        return
