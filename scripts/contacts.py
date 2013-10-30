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
import sys
import types

import ample_util
import csymmatch
import pdb_edit
import pdb_model
import residue_map


class ContactData(object):
    def __init__(self):
        self.numContacts = 0
        self.inregister = 0
        self.ooregister = 0
        self.origin = None
        self.allMatched = None
        self.ncontlog = None
        self.pdb = None
        self.csymmatchPdb = None
        return
    
    def __str__(self):
        """List the data attributes of this object"""
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
              isinstance(attr, types.FunctionType) ):
                me[slot] = attr
                
        s = "{0}\n".format( self.__repr__() )
        for k, v in me.iteritems():
            s += "{0} : {1}\n".format( k,v )
        return s

class Contacts(object):
    """Foo
    """
    
    def __init__( self ):
        self.ncontlog = None
        self.numContacts = 0
        self.inregister = 0
        self.ooregister = 0
        self.allMatched = []
        self.best = None
        return
    
    def getContacts(self, nativePdb=None, placedPdb=None, resSeqMap=None, nativeInfo=None, shelxePdb=None, workdir=None ):
        
#        try:
        self.run( nativePdb=nativePdb, placedPdb=placedPdb, resSeqMap=resSeqMap, nativeInfo=nativeInfo, shelxePdb=shelxePdb, workdir=workdir )
        self.parse_ncontlog()
#         except:
#             print "ERROR RUNNING CONTACTS: ",placedPdb
#             self.numContacts = -1
#             self.inregister = -1
#             self.ooregister = -1
#             self.allMatched = []
#             self.best = None
        return
    
    def run( self, nativePdb=None, placedPdb=None, resSeqMap=None, nativeInfo=None, shelxePdb=None, workdir=None ):
        """
        """

        self.workdir = workdir
        if not self.workdir:
            self.workdir = os.getcwd()
            
        pdbedit = pdb_edit.PDBEdit()
        
        if False:
            # Standardise pdb
            nativePdbStd = ample_util.filename_append( filename=nativePdb, astr="std" )
            pdbedit.standardise( inpdb=nativePdb, outpdb=nativePdbStd )
            nativePdb = nativePdbStd
            
            # Run a pass to find the # chains
            nativeInfo = pdbedit.get_info( nativePdb )
            
            if len( nativeInfo.models ) > 1:
                raise RuntimeError,"> 1 model!"
            
            # Check numbering matches and match numbering if necessary
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
 
        # Make a copy of placedPdb with chains renamed to lower case
        placedInfo = pdbedit.get_info( placedPdb )
        fromChain = placedInfo.models[0].chains
        toChain = [ c.lower() for c in fromChain ]
        placedAaPdb = ample_util.filename_append( filename=placedPdb, astr="ren", directory=self.workdir )
        pdbedit.rename_chains( inpdb=placedPdb, outpdb=placedAaPdb, fromChain=fromChain, toChain=toChain )

        # Get list of origins
        placedSpaceGroup = placedInfo.crystalInfo.spaceGroup
        if placedSpaceGroup != placedInfo.crystalInfo.spaceGroup:
            raise RuntimeError,"Mismatching space groups!"
        
        origins = pdb_model.alternateOrigins( placedSpaceGroup )
        #print "GOT ORIGINS ",origins
        
        # Add the shelxe origin to the list if it's not already in there
        csym = csymmatch.Csymmatch()
        csymmatchPdb = ample_util.filename_append( filename=shelxePdb, astr="csymmatch", directory=self.workdir )
        csym.run( refPdb=nativePdb, inPdb=shelxePdb, outPdb=csymmatchPdb )
        corig = csym.origin()
        
        if corig and corig not in origins:
            #print "{0} is not in {1}".format( corig, origins )
            origins.append( corig )
        
        # Loop over origins, move the placed pdb to the new origin and then run ncont
        self.best = ContactData()
        for i, origin in enumerate( origins  ):
            #print "GOT ORIGIN ",i,origin
            
            placedOriginPdb =  placedAaPdb
            if i != 0:
                # Move pdb to new origin
                placedOriginPdb = ample_util.filename_append( filename=placedAaPdb, astr="origin{0}".format(i), directory=self.workdir )
                pdbedit.translate( inpdb=placedAaPdb, outpdb=placedOriginPdb, ftranslate=origin )
            
            # Concatenate into one file
            joinedPdb = ample_util.filename_append( filename=placedOriginPdb, astr="joined", directory=self.workdir )
            pdbedit.merge( pdb1=nativePdb, pdb2=placedOriginPdb, pdbout=joinedPdb )
                
            # Run ncont
            # Need to get list of chains from Native as can't work out negate operator for ncont
            fromChain = nativeInfo.models[0].chains
            self.run_ncont( pdbin=joinedPdb, sourceChains=fromChain, targetChains=toChain )
            self.parse_ncontlog()
            
            if self.inregister + self.ooregister > self.best.inregister + self.best.ooregister:
                self.best.numContacts = self.numContacts
                self.best.inregister = self.inregister
                self.best.ooregister = self.ooregister
                self.best.allMatched = self.allMatched
                self.best.origin = origin
                self.best.ncontlog = self.ncontlog
                self.best.pdb = placedOriginPdb
                
            #print "GOT CONTACTS: {0} : {1} : {2}".format( self.numContacts, self.inregister, self.ooregister  )
            
        
        if self.best.pdb:
            # Just for info - run csymmatch so we can see the alignment
            csymmatchPdb = ample_util.filename_append( filename=self.best.pdb, astr="csymmatch_best", directory=self.workdir )
            csym.run( refPdb=nativePdb, inPdb=self.best.pdb, outPdb=csymmatchPdb, originHand=False )
            self.best.csymmatchPdb = csymmatchPdb
        else:
            self.best = None
#                 
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
            
        #print "LOG ",logfile
            
        self.numContacts = 0
        self.inregister = 0
        self.ooregister = 0
        self.allMatched = []
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
        
        #print "LINES ",clines
        
        contacts = [] # Tuples of: chainID, resSeq, chainID, resSeq, dist
        # Only collect contacts for central cell
        mycell = None
        for c in clines:
            fields = c.split()
            
            if len( fields ) != 15:
                # It seems we something get self contacts: 3GD8 - phaser_loc0_ALL_SCWRL_reliable_sidechains_trunc_12.483162_rad_3_UNMOD.1
                # don't really understand yet
                continue

            c1 = fields[ 0 ]
            r1 = fields[ 1 ]
            c2 = fields[ 6 ]
            r2 = fields[ 7 ]
            dist = float(  fields[ 12 ] )
            cell = int(  fields[ 13 ] )
            
            if not mycell:
                # initialise
                mycell = cell
            
            # central cell is 333
            if cell != mycell:
                #print "CHANGE OF CELL: {0} : {1}".format( cell, mycell )
                mycell = cell
                #raise RuntimeError,"Change of cell!"
            
            chainID1 = c1.split("/")[2]
            resSeq1 = int( r1.split("(")[0] )
            aa1 = r1.split("(")[1][:-2] # get amino acid and convert to single letter code
            aa1 = pdb_edit.three2one[ aa1 ] 
            chainID2 = c2.split("/")[2].upper() # We converted to lower case
            resSeq2 = int( r2.split("(")[0] )
            aa2 = r2.split("(")[1][:-2]
            aa2 = pdb_edit.three2one[ aa2 ]
            
            contacts.append( (chainID1, resSeq1, aa1, chainID2, resSeq2, aa2, dist, cell) )
    
        # Now count'em
        MINC = 3 # minimum contiguous to count
        last1=None
        last2=None
        count=0
        register=True # true if in register, false if out
        mycell=None
        thisMatched = []
        for i, (chainID1, resSeq1, aa1, chainID2, resSeq2, aa2, dist, cell) in enumerate( contacts ):
            
            #print "CONTACTS ",chainID1, resSeq1, aa1, chainID2, resSeq2, aa2, dist, cell
            # Initialise
            if i == 0:
                last1 = resSeq1
                last2 = resSeq2
                mycell = cell
                if resSeq1 != resSeq2:
                    register=False
                count = 1
                thisMatched = [ (chainID1, resSeq1, aa1, chainID2, resSeq2, aa2, dist, cell) ]
                continue
            
            # Is a contiguous residue and the register matches what we're reading
            # Make sure we don't count contiguous residues where the other residue doesn't match
            # Also where the cell changes
            if ( resSeq1 == last1 + 1 and resSeq2 == last2 + 1 ) \
                and  ( ( resSeq1 == resSeq2 and register ) or ( resSeq1 != resSeq2 and not register )  ) and mycell == cell:
                #print "INCREMENTING"
                count += 1
                thisMatched.append( (chainID1, resSeq1, aa1, chainID2, resSeq2, aa2, dist, cell) )
                last1 = resSeq1
                last2 = resSeq2
                
                # If this is the last one we want to drop through
                if i < len(contacts)-1:
                    continue
                
            if count > MINC:
                #  end of a contiguous sequence
                #print "ADDING ",count
                if register:
                    self.inregister += count
                else:
                    self.ooregister += count
                self.allMatched += thisMatched
            
            # Either starting again or a random residue
            if resSeq1 == resSeq2:
                register=True
            else:
                register=False
            
            last1 = resSeq1
            last2 = resSeq2
            mycell = cell
            thisMatched = [ (chainID1, resSeq1, aa1, chainID2, resSeq2, aa2, dist, cell) ]
            count = 1
        
        #print "GOT ALLMATCHED ",self.allMatched
            
        return

if __name__ == "__main__":
#     root = "/Users/jmht/Documents/AMPLE/data/ncont"
#     root = "/home/jmht/Documents/test/ncont/3PCV"
#     nativePdb = root + "/3PCV_std.pdb"
#     placedPdb = root + "/refmac_phaser_loc0_ALL_poly_ala_trunc_2.822761_rad_1_UNMOD.pdb"
#     refModelPdb = root + "/S_00000001.pdb"
#     
#     root = "/home/jmht/Documents/test/ncont/1D7M"
#     workdir = root + "/All_atom_trunc_5.131715_rad_2"
#     nativePdb = root + "/1D7M.pdb"
#     refModelPdb = root + "/S_00000001.pdb"
#     placedPdb = workdir + "/phaser_loc0_ALL_All_atom_trunc_5.131715_rad_2_UNMOD.1.pdb"
    
    
    root = "/home/jmht/Documents/test/ncont/3PCV"
    nativePdb = root + "/3PCV.pdb"
    refModelPdb = root + "/S_00000001.pdb"
    #workdir = root + "/All_atom_trunc_5.131715_rad_3"
    #placedPdb = workdir + "/phaser_loc0_ALL_All_atom_trunc_5.131715_rad_3_UNMOD.1.pdb"
    
    workdir = "/home/jmht/Documents/test/ncont/new/3GD8"
    nativePdb = workdir + "/3GD8_std.pdb"
    placedPdb = workdir + "/phaser_loc0_ALL_SCWRL_reliable_sidechains_trunc_12.483162_rad_3_UNMOD.1.pdb"
    refModelPdb = workdir + "/S_00000001.pdb"
    
    os.chdir(workdir)
    
    # PE = pdb_edit.PDBEdit()
    # PE.cat_pdbs(pdb1=nativePdb, pdb2="csymmatchAll.pdb", pdbout="joined.pdb")
    # 
    
    c = Contacts()
    c.run( nativePdb, placedPdb, refModelPdb, workdir=workdir )
    print c.best
    #logfile = "/home/jmht/Documents/test/ncont/3PCV/poly_ala_trunc_2.822761_rad_1/phaser_loc0_ALL_poly_ala_trunc_2.822761_rad_1_UNMOD.1_csymmatch_ren_joined.pdb.ncont.log"
    #logfile = "/home/jmht/Documents/test/ncont/3PCV/All_atom_trunc_5.131715_rad_3/phaser_loc0_ALL_All_atom_trunc_5.131715_rad_3_UNMOD.1_csymmatch_ren_joined.pdb.ncont.log"
    #logfile = "/Users/jmht/Documents/AMPLE/data/ncont/3GD8/phaser_loc0_ALL_All_atom_trunc_11.13199_rad_3_UNMOD.1_reseq_ren_origin2_joined.pdb.ncont.log"
    #c.parse_ncontlog( logfile=logfile )
    