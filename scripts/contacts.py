#!/usr/bin/env python


"""



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


class Contacts(object):
    """Foo
    """
    
    def __init__( self ):
        self.ncontlog = None
        self.numContacts = 0
        self.inregister = 0
        self.ooregister = 0
        return
    
    def getContacts(self, nativePdb, placedPdb, refModelPdb, workdir=None ):
        
        try:
            self.run(nativePdb, placedPdb, refModelPdb, workdir=workdir )
            self.parse_ncontlog()
        except:
            print "ERROR RUNNING CONTACTS: ",placedPdb
            self.numContacts = -1
            self.inregister = -1
            self.ooregister = -1
        return
    
    def run( self, nativePdb, placedPdb, refModelPdb, workdir=None ):
        """
        """

        self.workdir = workdir
        if not self.workdir:
            self.workdir = os.getcwd()
            
        pdbedit = pdb_edit.PDBEdit()
        
        # Standardise pdb
        #nativePdbStd = ample_util.filename_append( filename=nativePdb, astr="std" )
        #pdbedit.standardise( inpdb=nativePdb, outpdb=nativePdbStd )
        #nativePdb = nativePdbStd
        
        # Run a pass to find the # chains
        nativeInfo = pdbedit.get_info( nativePdb )
        
        if len( nativeInfo.models ) > 1:
            raise RuntimeError,"> 1 model!"
        #native_chains = nativeInfo.models[ 0 ].chains
        #placed_chains = placedInfo.models[ 0 ].chains # only ever one model in the refined pdb
        
#         # Check numbering matches and match numbering if necessary
        resSeqMap = residue_map.residueSequenceMap()
        modelInfo = pdbedit.get_info( refModelPdb )
        # NEED TO FIX NAMING AS THIS IS WAY TOO CONFUSING!
        resSeqMap.fromInfo( nativeInfo=modelInfo,
                            nativeChainID=modelInfo.models[0].chains[0],
                            modelInfo=nativeInfo,
                            modelChainID=nativeInfo.models[0].chains[0]
                            )
         
        if not resSeqMap.resSeqMatch():
            #print "NUMBERING DOESN'T MATCH"
            #raise RuntimeError,"NUMBERING DOESN'T MATCH"
            # We need to create a copy of the placed pdb with numbering matching the native
            placedPdbRes = ample_util.filename_append( filename=placedPdb, astr="reseq", directory=self.workdir )
            pdbedit.match_resseq( targetPdb=placedPdb, sourcePdb=None, outPdb=placedPdbRes, resMap=resSeqMap )
            placedPdb = placedPdbRes
 
        # Output filename
        csymmatchPdb = ample_util.filename_append( filename=placedPdb,
                                                   astr="csymmatch",
                                                   directory=self.workdir
                                                    )
                
        # Determine the best orientation/origin
        self.run_csymmatch( placedPdb=placedPdb,
                            nativePdb=nativePdb,
                            csymmatchPdb=csymmatchPdb
                            )
        
            
        # Now rename csymmatchPdb chains
        csymmatchInfo = pdbedit.get_info( csymmatchPdb )
        fromChain = csymmatchInfo.models[0].chains
        toChain = [ c.lower() for c in fromChain ]
        
        csymmatchAaPdb = ample_util.filename_append( filename=csymmatchPdb, astr="ren" )
        pdbedit.rename_chains( inpdb=csymmatchPdb, outpdb=csymmatchAaPdb, fromChain=fromChain, toChain=toChain )
        
        # Concatenate into one file
        joinedPdb = ample_util.filename_append( filename=csymmatchAaPdb, astr="joined" )
        pdbedit.cat_pdbs( pdb1=nativePdb, pdb2=csymmatchAaPdb, pdbout=joinedPdb )
        
        # Run ncont
        self.run_ncont( pdbin=joinedPdb, sourceChains=fromChain, targetChains=toChain )
                
        return
    
    
    def run_over_chains( self, nativePdb, placedPdb, refModelPdb, workdir=None ):
        """
        """

        self.workdir = workdir
        if not self.workdir:
            self.workdir = os.getcwd()
            
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
                                    modelChainID=modelInfo.models[0].chains[0]
                                    )
                
                if not resSeqMap.resSeqMatch():
                    #print "NO MATCH"
                    # We need to create a copy of the csymmatch pdb with numbering matching the native
                    csymmatchRenumberPdb = ample_util.filename_append( filename=csymmatchPdb, astr="_renumber" )
                    pdbedit.match_resseq( nativeChainPdb, csymmatchPdb, csymmatchRenumberPdb, resMap=resSeqMap )
                    csymmatchPdb = csymmatchRenumberPdb
                    
                # Now rename csymmatchPDB chain to X
                csymmatchXPdb = ample_util.filename_append( filename=csymmatchPdb, astr="X" )
                pdbedit.extract_chain( csymmatchPdb, csymmatchXPdb, chainID=placedChainID, newChainID='X', renumber=False )
                
                # Concatenate into one file
                joinedPdb = ample_util.filename_append( filename=csymmatchXPdb, astr="joined" )
                pdbedit.cat_pdbs( pdb1=nativeChainPdb, pdb2=csymmatchXPdb, pdbout=joinedPdb)
                
        return

    def run_csymmatch( self, placedPdb=None, nativePdb=None, csymmatchPdb=None ):
        """FOO
        """
        
        logfile = csymmatchPdb +".log"
        cmd= [ 'csymmatch',
              "-pdbin-ref",
              nativePdb,
              "-pdbin",
              placedPdb,
              "-pdbout",
              csymmatchPdb,
              "-origin-hand" ]
        
#         "-connectivity-radius",
#         "100",

        retcode = ample_util.run_command(cmd=cmd, logfile=logfile, dolog=False)
        
        if retcode != 0:
            raise RuntimeError, "Error running command: {0}".format( " ".join(cmd) )
        return
    
    def run_ncont( self, pdbin=None, sourceChains=None, targetChains=None, maxdist=1.5 ):
        """FOO
        """
        
        self.ncontlog = pdbin +".ncont.log"
        cmd = [ "ncont", "xyzin", pdbin ]
        
        # Build up stdin
        stdin = ""
        stdin += "source {0}//CA\n".format( ",".join( sourceChains )  )  
        stdin += "target {0}//CA\n".format( ",".join( targetChains )  )  
        stdin += "maxdist {0}\n".format( maxdist )
        stdin += "cells 2\n"
        
        retcode = ample_util.run_command( cmd=cmd, logfile=self.ncontlog, directory=os.getcwd(), dolog=False, stdin=stdin )
        
        if retcode != 0:
            raise RuntimeError,"Error running ncont"
        
        return
    
    def parse_ncontlog(self, logfile=None ):
        
        if not logfile:
            logfile = self.ncontlog
            
        self.numContacts = 0
        clines = []
        
        capture=False
        with open( logfile, 'r' ) as f:#
            while True:
                line = f.readline().strip()
                
                if capture and not line:
                    break
                
                if "contacts found:" in line:
                    self.numContacts = int( line.split()[0] )
                
                if "NO CONTACTS FOUND." in line:
                    return False
                
                if "SOURCE ATOMS" in line:
                    capture=True
                    f.readline() # skip blank line
                    continue
                
                if capture:
                    clines.append( line )
            
        assert self.numContacts == len(clines)
        
        
        contacts = [] # Tuples of: chainID, resSeq, chainID, resSeq, dist
        # Only collect contacts for central cell
        for c in clines:
            fields = c.split()
            cell = int( fields[ 13 ] )
            if cell != 333:
                continue
            
            c1 = fields[ 0 ]
            r1 = fields[ 1 ]
            c2 = fields[ 6 ]
            r2 = fields[ 7 ]
            dist = float(  fields[ 12 ] )
            
            chainID1 = c1.split("/")[2]
            resSeq1 = int( r1.split("(")[0] )
            chainID2 = c2.split("/")[2].upper() # We converted to lower case
            resSeq2 = int( r2.split("(")[0] )
            
            contacts.append( (chainID1, resSeq1, chainID2, resSeq2, dist) )
    
        # Now count'em
        MINC = 3 # minimum contiguous to count
        self.inregister=0
        self.ooregister=0
        last1=None
        last2=None
        count=0
        register=True # true if in register, false if out
        for (chainID1, resSeq1, chainID2, resSeq2, dist) in contacts:
            
            # First
            if last1 == None:
                last1 = resSeq1
                last2 = resSeq2
                if resSeq1 != resSeq2:
                    register=False
                count+=1
                continue
            
            if resSeq1 == last1 + 1 and resSeq2 == last2 + 1:
                # Make sure we don't count continguous residues where the other residue doesn't match
                if ( resSeq1 != resSeq2 and register ) or ( resSeq1 == resSeq2 and not register ):
                    pass
                else:
                    #print "INCREMENTING"
                    count += 1
                    last1 = resSeq1
                    last2 = resSeq2
                    continue
                
            if count > MINC:
                #  end of a contiguous sequence
                #print "ADDING ",count
                if register:
                    self.inregister += count
                else:
                    self.ooregister += count
            
            # Either starting again or a random residue
            if resSeq1 == resSeq2:
                register=True
            else:
                register=False
            
            last1 = resSeq1
            last2 = resSeq2
            count = 1
        
        return

if __name__ == "__main__":
    root = "/Users/jmht/Documents/AMPLE/data/ncont"
    root = "/home/jmht/Documents/test/ncont/3PCV"
    nativePdb = root + "/3PCV_std.pdb"
    placedPdb = root + "/refmac_phaser_loc0_ALL_poly_ala_trunc_2.822761_rad_1_UNMOD.pdb"
    refModelPdb = root + "/S_00000001.pdb"
    
    root = "/home/jmht/Documents/test/ncont/1D7M"
    workdir = root + "/All_atom_trunc_5.131715_rad_2"
    nativePdb = root + "/1D7M.pdb"
    refModelPdb = root + "/S_00000001.pdb"
    placedPdb = workdir + "/phaser_loc0_ALL_All_atom_trunc_5.131715_rad_2_UNMOD.1.pdb"
    
    
    root = "/home/jmht/Documents/test/ncont/3PCV"
    nativePdb = root + "/3PCV.pdb"
    refModelPdb = root + "/S_00000001.pdb"
    
    workdir = root + "/All_atom_trunc_5.131715_rad_2"
    placedPdb = workdir + "/phaser_loc0_ALL_All_atom_trunc_5.131715_rad_2_UNMOD.1.pdb"
    
    workdir = root + "/poly_ala_trunc_2.822761_rad_1"
    placedPdb = workdir + "/phaser_loc0_ALL_poly_ala_trunc_2.822761_rad_1_UNMOD.1.pdb"
    
    workdir = root + "/All_atom_trunc_5.131715_rad_3"
    placedPdb = workdir + "/phaser_loc0_ALL_All_atom_trunc_5.131715_rad_3_UNMOD.1.pdb"
    
    os.chdir(workdir)
    
    # PE = pdb_edit.PDBEdit()
    # PE.cat_pdbs(pdb1=nativePdb, pdb2="csymmatchAll.pdb", pdbout="joined.pdb")
    # 
    # import sys
    # sys.exit()
    
    c = Contacts()
    #c.run( nativePdb, placedPdb, refModelPdb, workdir=workdir )
    logfile = "/home/jmht/Documents/test/ncont/3PCV/poly_ala_trunc_2.822761_rad_1/phaser_loc0_ALL_poly_ala_trunc_2.822761_rad_1_UNMOD.1_csymmatch_ren_joined.pdb.ncont.log"
    logfile = "/home/jmht/Documents/test/ncont/3PCV/All_atom_trunc_5.131715_rad_3/phaser_loc0_ALL_All_atom_trunc_5.131715_rad_3_UNMOD.1_csymmatch_ren_joined.pdb.ncont.log"
    c.parse_ncontlog( logfile=logfile )
    print c.ncontacts, c.inregister,c.ooregister
