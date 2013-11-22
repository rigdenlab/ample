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
import types
import unittest

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
        self.backwards = 0
        self.origin = None
        self.allMatched = None
        self.ncontLog = None
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
        self.ncontLog = None
        self.numContacts = 0
        self.contacts = None
        self.inregister = 0
        self.ooregister = 0
        self.backwards = 0
        self.allMatched = []
        self.best = None
        # testing
        self.originCompare = {}
        return
    
    def getContacts(self, nativePdb=None, placedPdb=None, resSeqMap=None, nativeInfo=None, shelxePdb=None, workdir=None ):
        
        self.run( nativePdb=nativePdb, placedPdb=placedPdb, resSeqMap=resSeqMap, nativeInfo=nativeInfo, shelxePdb=shelxePdb, workdir=workdir )
        
        return

    def helixFromContacts( self, dsspP=None ):
        """Return the sequence of the longest contiguous helix from the given contact data"""
        
        #print "GOT DATA ",self.best.allMatched
        #print "GOT DSSP ",dsspP.asDict()
        
        maxCounts = [] # list of maximum counts of contiguous helices in each group
        maxIndices = [] # Array of (start,stop) tuples for the max contiguous sequence in each group
        
        # Assign the secondary structure & work out the largest contiguous group
        for i, contactGroup in enumerate( self.best.allMatched ):
            count = 0
            thisCounts = []
            thisIndices = []
            start = 0 
            stop = 0
            #print "GROUP"
            for ic, c in enumerate( contactGroup ):
                
                ( chainId1, resSeq1, aa1, chainId2, resSeq2, aa2, dist, cell, symmetry ) = c
                
                #ss = getSS( c, dsspP )
                #print "DATA ",chainId1, resSeq1, aa1, chainId2, resSeq2, aa2, dist, cell
                ss = dsspP.getAssignment( resSeq1, aa1, chainId1 )
                
                #print "Looping through ", chainId1, resSeq1, aa1,ss
                if ss == 'H':
                    count += 1
                    stop = ic
                
                if ic == len( contactGroup ) - 1 or ss != 'H':
                    thisCounts.append( count )
                    thisIndices.append( ( start, stop ) )
                    count = 0
                    start = ic + 1
                    stop = ic
                    
            # Now add the maximum for that group
            #print "this Counts ",thisCounts
            #print "thisIndices ",thisIndices
            mc =  max( thisCounts )
            maxCounts.append( mc )
            
            maxIndices.append( thisIndices[ thisCounts.index( mc ) ] )
                    
            
        #print "GOT maxCounts ",maxCounts
        #print "GOT maxIndices ",maxIndices
        
        # Get the index of the group with the largest count
        gmax = maxCounts.index( max( maxCounts ) )
        start = maxIndices[ gmax ][0]
        stop = maxIndices[ gmax ][1]
        cg = self.best.allMatched[ gmax ]
        
        # Now get the sequence
        sequence = ""
        for i in range( start, stop + 1 ):
            c = cg[ i ]
            (chainId1, resSeq1, aa1, chainId2, resSeq2, aa2, dist, cell, symmetry ) = c
            sequence += aa1
            
        #print "GOT sequence ",sequence
        return sequence


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
            resSeqMap.fromInfo( refInfo=modelInfo,
                                refChainID=modelInfo.models[0].chains[0],
                                targetInfo=nativeInfo,
                                targetChainID=nativeInfo.models[0].chains[0]
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
        # Pythonic way of checking if any of the origins are floating
        floating = any(  map( lambda o: 'x' in o or 'y' in o or 'z' in o, origins  ) )
        
        # Add the shelxe origin to the list if it's not already in there
        csym = csymmatch.Csymmatch()
        corig = None
        if shelxePdb:
            csymmatchPdb = ample_util.filename_append( filename=shelxePdb, astr="csymmatch", directory=self.workdir )
            csym.run( refPdb=nativePdb, inPdb=shelxePdb, outPdb=csymmatchPdb )
            corig = csym.origin()
        
        #if not corig:
        #    print "NO CSYMMATCH ORIGIN"
        
        # For floating origins we use the csymmatch origin
        if floating:
            if not corig:
                # If csymmatch failed, we can't do owt
                self.best = None
                return False
            # Should check if the origin is acceptable, but that would require checking through all the 
            # alternate origins and seeing if tne non-floating axes had acceptable values  
            origins = [ corig ]
        else:
            if corig and corig not in origins:
                #print "csymmatch origin {0} is not in origins {1}".format( corig, origins )
                origins.append( corig )
        
        # Loop over origins, move the placed pdb to the new origin and then run ncont
        self.best = ContactData()
        self.originCompare = {}
        for i, origin in enumerate( origins  ):
            #print "GOT ORIGIN ",i,origin
            
            placedOriginPdb =  placedAaPdb
            if origin != [ 0.0, 0.0, 0.0 ]:
                # Move pdb to new origin
                #ostr="origin{0}".format(i)
                ostr="o{0}".format( origin ).replace(" ","" )
                placedOriginPdb = ample_util.filename_append( filename=placedAaPdb, astr=ostr, directory=self.workdir )
                pdbedit.translate( inpdb=placedAaPdb, outpdb=placedOriginPdb, ftranslate=origin )
            
            # Concatenate into one file
            joinedPdb = ample_util.filename_append( filename=placedOriginPdb, astr="joined", directory=self.workdir )
            pdbedit.merge( pdb1=nativePdb, pdb2=placedOriginPdb, pdbout=joinedPdb )
                
            # Run ncont
            # Need to get list of chains from Native as can't work out negate operator for ncont
            fromChain = nativeInfo.models[0].chains
            self.runNcont( pdbin=joinedPdb, sourceChains=fromChain, targetChains=toChain )
            self.parseNcontLog()
            self.countContacts()
            
            # Just for debugging
            self.originCompare[ "{0}".format( origin ).replace(" ","" ) ] = self.inregister + self.ooregister
            
            if self.inregister + self.ooregister > self.best.inregister + self.best.ooregister:
                self.best.numContacts = self.numContacts
                self.best.inregister = self.inregister
                self.best.ooregister = self.ooregister
                self.best.backwards = self.backwards
                self.best.contacts = self.contacts
                self.best.allMatched = self.allMatched
                self.best.origin = origin
                self.best.ncontLog = self.ncontLog
                self.best.pdb = placedOriginPdb
                
            #print "GOT CONTACTS: {0} : {1} : {2}".format( self.numContacts, self.inregister, self.ooregister  )
        
        if self.best.pdb:
            # Just for info - run csymmatch so we can see the alignment
            csymmatchPdb = ample_util.filename_append( filename=self.best.pdb, astr="csymmatch_best", directory=self.workdir )
            csym.run( refPdb=nativePdb, inPdb=self.best.pdb, outPdb=csymmatchPdb, originHand=False )
            self.best.csymmatchPdb = csymmatchPdb
            #if self.best.origin != corig:
            #    print "GOT DIFFERENT BEST ORIGIN {0} {1} FOR {2}".format( self.best.origin, corig, placedAaPdb )
            return True
        else:
            self.best = None
            return False

    def runNcont( self, pdbin=None, sourceChains=None, targetChains=None, maxdist=1.5 ):
        """FOO
        """
        
        self.ncontLog = pdbin +".ncont.log"
        cmd = [ "ncont", "xyzin", pdbin ]
        
        # Build up stdin
        stdin = ""
        stdin += "source {0}//CA\n".format( ",".join( sourceChains )  )  
        stdin += "target {0}//CA\n".format( ",".join( targetChains )  )  
        stdin += "maxdist {0}\n".format( maxdist )
        stdin += "cells 2\n"
        stdin += "sort target inc\n"
        
        retcode = ample_util.run_command( cmd=cmd, logfile=self.ncontLog, directory=os.getcwd(), dolog=False, stdin=stdin )
        
        if retcode != 0:
            raise RuntimeError,"Error running ncont"
        
        return
    
    def parseNcontLog( self, logfile=None ):
        """
        
        Lines are of format
        /1/B/1042(MET). / CA [ C]:  /1/b/ 988(GLU). / CA [ C]:   1.09 223 X-1/2,Y-1/2,Z
 
        """
        
        if not logfile:
            logfile = self.ncontLog
        
        #print "LOG ",logfile
            
        self.numContacts = 0
        self.contacts = None
        clines = []
        
        capture=False
        with open( logfile, 'r' ) as f:#
            while True:
                line = f.readline().rstrip()
                
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

        # Got data lines so now extract data
        # Could probably just do this in the reading loop now    
        contacts = [] 
        lastSource = None
        for c in clines:
            
            # Reconstruct lines with only the second contact using the data from the corresponding last complete line
            if not c[0:29].strip():
                #print "MATCH ACROSS TWO ATOMS"
                c = lastSource[0:29] + c[29:]
            else:
                lastSource = c
            
            # As the numbers overrun we can't split so we assume fixed format
            chainId1 = c[4]
            resSeq1 = int( c[6:10].strip() )
            aa1 = c[11:14]
            aa1 = pdb_edit.three2one[ aa1 ] # get amino acid and convert to single letter code
            chainId2 = c[32]
            resSeq2 = int( c[34:38].strip() )
            aa2 = c[39:42]
            aa2 = pdb_edit.three2one[ aa2 ]
            dist = float( c[56:62].strip() )
            cell = int( c[63:66])
            symmetry = c[67:]
            
            contacts.append( (chainId1, resSeq1, aa1, chainId2, resSeq2, aa2, dist, cell, symmetry ) )

#         # Put in dictionary by chains for easy access
#         contacts = {}
#         for c in contactsList:
#             
#             chainId1, resSeq1, aa1, chainId2, resSeq2, aa2, dist, cell, symmetry = c
#             
#             if not contacts.has_key( chainId1 ):
#                 # First entry in first source chain
#                 contacts[ chainId1 ] = { chainId2: [ c ] }
#                 continue
#             
#             if not contacts[ chainId1 ].has_key( chainId2 ):
#                 # Adding a new target chain
#                 contacts[ chainId1 ][ chainId2 ] = [ c ]
#                 continue
#             
#             # Adding to existing source & target chains
#             contacts[ chainId1 ][ chainId2 ].append( c )
#                 
# #         for sc in contacts.keys():
# #             for tc in contacts[ sc ]:
# #                 for c in contacts[ sc ][ tc ]:
# #                     print c
# 
        self.contacts = contacts
                    
        return contacts
    
    def countContacts( self ):
        
        if not self.contacts:
            return
        
        # Now count'em and put them into groups
        MINC = 3 # minimum contiguous to count
        self.inregister = 0
        self.ooregister = 0
        self.backwards = 0
        self.allMatched = []

        last1 = None
        last2 = None
        count = None
        register = True # true if in register, false if out
        backwards = False
        thisMatched = []
        
        for i, c in enumerate( self.contacts ):
            
            chainId1, resSeq1, aa1, chainId2, resSeq2, aa2, dist, cell, symmetry = c

            #print "CONTACT ",i, chainId1, resSeq1, aa1, chainId2, resSeq2, aa2, dist, cell, symmetry
            # Initialise
            if i == 0:
                last1 = resSeq1
                last2 = resSeq2
                if resSeq1 != resSeq2:
                    register = False
                count = 1
                thisMatched = [ ( chainId1, resSeq1, aa1, chainId2, resSeq2, aa2, dist, cell, symmetry ) ]
                continue
            
            # LOGIC HERE STILL NEEDS WORK
            # We are reading contiguous matches - forwards or backwards
            if ( resSeq1 == last1 + 1 and resSeq2 == last2 + 1 ) or \
               ( resSeq1 == last1 + 1 and resSeq2 == last2 - 1 ) or ( resSeq1 == last1 - 1 and resSeq2 == last2 + 1 ):
                
                # either in or oo register read
                if ( resSeq1 == resSeq2 and register ) or ( resSeq1 != resSeq2 and not register ):
                    
                    # Check if this is a change or part of a stretch
                    if ( ( resSeq1 == last1 + 1 and resSeq2 == last2 - 1 ) or ( resSeq1 == last1 - 1 and resSeq2 == last2 + 1 ) ) and not backwards:
                        backwards = True
                        #print "BACKWARDS ", self.ncontLog
                        
                    elif resSeq1 == last1 + 1 and resSeq2 == last2 + 1 and backwards:
                        backwards = False
                        
                    else:
                        count += 1
                        thisMatched.append( ( chainId1, resSeq1, aa1, chainId2, resSeq2, aa2, dist, cell, symmetry ) )
                        last1 = resSeq1
                        last2 = resSeq2
                        
                        # If this is the last one we want to drop through
                        if i < len( self.contacts ) - 1:
                            continue
                
            # Anything that doesn't continue didn't match
                
            if count >= MINC:
                #  end of a contiguous sequence
                if register:
                    assert not backwards
                    self.inregister += count
                else:
                    #print "adding {0} to ooregister {1}".format( count, self.ooregister )
                    self.ooregister += count
                    if backwards:
                        print "ADDING {0} to backwards log {1}".format( count, self.ncontLog )
                        self.backwards += count
                        
                self.allMatched.append( thisMatched )
            
            # Either starting again or a random residue
            if resSeq1 == resSeq2:
                register=True
            else:
                register=False
            
            last1 = resSeq1
            last2 = resSeq2
            thisMatched = [ (chainId1, resSeq1, aa1, chainId2, resSeq2, aa2, dist, cell, symmetry ) ]
            count = 1
        
        #print "GOT ALLMATCHED ",self.allMatched
            
        return
    
    def writeHelixFile(self, filename=None, dsspP=None ):
        if self.best:
            sequence = self.helixFromContacts( dsspP )
            with open( filename, 'w' ) as f:
                f.write( sequence+"\n" )
        
        return
    
    
class TestContacts( unittest.TestCase ):
    
    def setUp(self):
        
        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        self.ampleDir = os.sep.join( paths[ : -1 ] )
        self.testfilesDir = os.sep.join( paths[ : -1 ] + [ 'tests', 'testfiles' ] )
        
        return
    
    
    def testParse2(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont2.log" )
        
        c = Contacts()
        contacts = c.parseNcontLog( logfile=logfile )
        c.countContacts()
        
        self.assertEqual( c.numContacts, 10 )
        self.assertEqual( c.inregister, 0 )
        self.assertEqual( c.ooregister, 7 )
        self.assertEqual( c.backwards, 7 )
        
        return
    
    def testParse3(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont3.log" )
        
        c = Contacts()
        contacts = c.parseNcontLog( logfile=logfile )
        c.countContacts()
        
        self.assertEqual( c.numContacts, 14 )
        self.assertEqual( c.inregister, 0 )
        self.assertEqual( c.ooregister, 10 )
        self.assertEqual( c.backwards, 0 )
        
        return
    
    def testParse4(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont4.log" )
        
        c = Contacts()
        contacts = c.parseNcontLog( logfile=logfile )
        c.countContacts()
        
        self.assertEqual( c.numContacts, 56 )
        self.assertEqual( c.inregister, 0 )
        self.assertEqual( c.ooregister, 55 )
        self.assertEqual( c.backwards, 0 )
        
        return

if __name__ == "__main__":
    
    unittest.main()
    
    import sys
    sys.exit()


