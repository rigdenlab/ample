#!/usr/bin/env python


"""

get residue map
for each chain in native
    for each chain in placed
        run_csymmatch
        # Gives pdb with oriented model
        create copy of csymmatch pdb with residue numbering matching native and rename chain to X
        concatenate the native and reformatted csymmatch pdb
        run ncont to generate contacts
        parse ncont file to generate map & analyse whether placed bits match and what type of structure they are

"""

import os

import ample_util
import pdb_edit
import residue_map


class Csymmatch(object):
    """Foo
    """
    
    def __init__( self ):
        return

    def run( self, nativePdb, placedPdb, refModelPdb, workdir=None ):
        """
        """

        self.workdir = workdir
        if not self.workdir:
            self.workdir = os.getcwd()
            
        
        self.refModelPdb = refModelPdb

        # Run a pass to find the # chains
        pdbedit = pdb_edit.PDBEdit()
        placedInfo = pdbedit.get_info( placedPdb )
        nativeInfo = pdbedit.get_info( nativePdb )
        native_chains = nativeInfo.models[ 0 ].chains
        placed_chains = placedInfo.models[ 0 ].chains # only ever one model in the refined pdb
        
        #print "got native chains ", native_chains
        #print "got refined chains ", refined_chains
            
        for nativeChainID in native_chains:
            #print "native_chain: {0}".format( nativeChainID )
                    
            if len( native_chains ) == 1:
                # Don't need to do owt as we are just using the native as is
                nativeChainPdb = nativePdb
            else:
                
                # Extract the chain from the pdb
                nativeChainPdb = ample_util.filename_append( filename=nativePdb, 
                                                             astr="chain{0}".format( nativeChainID ), 
                                                             directory=self.workdir
                                                             )
                pdbedit.extract_chain( nativePdb, nativeChainPdb, chainID=nativeChainID )
                
                assert os.path.isfile( nativeChainPdb )
            
            for placedChainID in placed_chains:
                
                #print "refined_chain: {0}".format( placedChainID )
                
                assert os.path.isfile( nativeChainPdb )
                
                # Extract the chain from the pdb
                placedChainPdb = ample_util.filename_append( filename=placedPdb,
                                                             astr="chain{0}".format( placedChainID ),
                                                             directory=self.workdir
                                                              )
                pdbedit.extract_chain( placedPdb, placedChainPdb, chainID=placedChainID, newChainID=nativeChainID )
                
                # Output filename
                astr = "chain{0}_csymmatch".format( placedChainID ) 
                csymmatchPdb = ample_util.filename_append( filename=nativeChainPdb,
                                                           astr=astr,
                                                           directory=self.workdir
                                                            )
                
                # Determine the best orientation/origin
                self.run_csymmatch( placedPdb=placedChainPdb,
                                    nativePdb=nativeChainPdb,
                                    csymmatchPdb=csymmatchPdb,
                                    nativeChainID=nativeChainID )
                
                # Check whether the residue numbers match
                resSeqMap = residue_map.residueSequenceMap()
                modelInfo = pdbedit.get_info( refModelPdb )
                resSeqMap.fromInfo( nativeInfo=nativeInfo,
                                    nativeChainID=nativeChainID,
                                    modelInfo=modelInfo,
                                    modelChainID=placedChainID
                                    )
                
                if not resSeqMap.resSeqMatch():
                    # We need to create a copy of the csymmatch pdb with numbering matching the native
                    csymmatchRenumberPdb = ample_util.filename_append( filename=csymmatchPdb, astr="_renumber" )
                    PE.match_resseq( nativeChainPdb, csymmatchPdb, csymmatchRenumberPdb, resMap=resMap )
                    csymmatchPdb = csymmatchRenumberPdb
                    
                # Now rename csymmatchPDB chain to X
                csymmatchXPdb = ample_util.filename_append( filename=csymmatchPdb, astr="X" )
                pdbedit.extract_chain( csymmatchPdb, csymmatchXPdb, chainID=placedChainID, newChainID='X', renumber=False )
                
                # Concatenate into one file
                joinedPdb = ample_util.filename_append( filename=csymmatchXPdb, astr="joined" )
                pdbedit.cat_pdbs( pdb1=nativeChainPdb, pdb2=csymmatchXPdb, pdbout=joinedPdb)

                
        return

    def run_csymmatch( self, placedPdb=None, nativePdb=None, csymmatchPdb=None, nativeChainID=None ):
        """FOO
        """
        
        workdir=os.getcwd()
        
#         resSeqMap = residue_map.residueSequenceMap( nativePdb, self.refModelPdb )
#         if not resSeqMap.resSeqMatch():
#             
#             # We need to create a copy of the native with numbering matching the model
#             n = os.path.splitext( os.path.basename( nativePdb ) )[0]
#             nativeRenumber = "{0}_renumber.pdb".format( n )
#             PE.match_resseq( nativePdb, model, nativeRenumber, resMap=resSeqMap )
#             
#             nativePdb = nativeRenumber

        # Now run csymmatch


        logfile = csymmatchPdb +".log"
        cmd= [ 'csymmatch',
              "-pdbin-ref",
              nativePdb,
              "-pdbin",
              placedPdb,
              "-pdbout",
              csymmatchPdb,
              "-connectivity-radius",
              "100",
              "-origin-hand" ]
        
        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=workdir, dolog=False)
        
        if retcode != 0:
            raise RuntimeError, "Error running command: {0}".format( " ".join(cmd) )
        return
    

nativePdb = "/Users/jmht/Documents/AMPLE/data/ncont/3PCV/3PCV.pdb"
placedPdb = "/Users/jmht/Documents/AMPLE/data/ncont/3PCV/refmac_phaser_loc0_ALL_poly_ala_trunc_2.822761_rad_1_UNMOD.pdb"
refModelPdb = "/Users/jmht/Documents/AMPLE/data/ncont/3PCVS_00000001.pdb"
workdir="/Users/jmht/Documents/AMPLE/data/ncont/3PCV"

c = Csymmatch()
c.run( nativePdb, placedPdb, refModelPdb, workdir=workdir )
