#!/usr/bin/env python


"""



run_csymmatch
# Gives pdb with oriented model
create copy of csymmatch pdb with residue numbering matching native and rename chain to X
concatenate the native and reformatted csymmatch pdb
run ncont to generate contacts
parse ncont file to generate map & analyse whether placed bits match and what type of structure they are

"""

from operator import itemgetter
import os
import types
import unittest

import ample_util
import csymmatch
import dssp
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
        self.floatingOrigin = None
        self.csymmatchOrigin = None
        self.ncontLog = None
        self.pdb = None
        self.csymmatchPdb = None
        self.helix = None
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
        self.best = None
        # testing
        self.originCompare = {}
        return
    
    def getContacts(self, nativePdb=None, placedPdb=None, resSeqMap=None, nativeInfo=None, shelxePdb=None, workdir=None, dsspLog=None ):
        
        if not self.run( nativePdb=nativePdb, placedPdb=placedPdb, resSeqMap=resSeqMap, nativeInfo=nativeInfo, shelxePdb=shelxePdb, workdir=workdir ):
            return False
        
        
        #print "placedPdb ",placedPdb
        #print "nativePdb ",nativePdb
        #print "GOT MAP ",resSeqMap
        
        # We should have contact data to calculate a helix from
        assert dsspLog
        dsspP = dssp.DsspParser( dsspLog )
        sequence = self.helixFromContacts( dsspP=dsspP )
        if sequence:
            self.best.helix = sequence
        
        if self.best.pdb:
            
            csym = csymmatch.Csymmatch()
            # Just for info - run csymmatch so we can see the alignment
            csymmatchPdb = ample_util.filename_append( filename=self.best.pdb, astr="csymmatch_best", directory=self.workdir )
            csym.run( refPdb=nativePdb, inPdb=self.best.pdb, outPdb=csymmatchPdb, originHand=False )
            self.best.csymmatchPdb = csymmatchPdb
        
        return

    def helixFromContacts( self, dsspP=None, contacts=None, ):
        """Return the sequence of the longest contiguous helix from the given contact data
        
        
        source is the native, target the model
        
        Get start and stop indices of all contiguous chunks
        - we can match multiple chains in the model, but any match must be within a single native chain
        - the source can increment or decrement - the model only ever increments (ncont "sort target inc" argument)
        - matches can be in-register or out-of-register, but for finding the longest chunk we don't care
        
        startstop = [ (10,15), (17, 34), (38, 50) ]
        
        # Loop through all chunks, and if any two have a gap of < mingap, join the indices together.
        
        # Get the indices of the largest chunk
        
        # Get the corresponding AA sequence
        
        """
        
        #print "GOT DATA ",contacts
        #print "GOT DSSP ",dsspP.asDict()
        
        if contacts is None:
            contacts = self.best.contacts
            
        if not len( contacts ):
            return None
            
        #
        # Loop through the contacts finding the start, stop indices in the list of contacts of contiguous chunks
        #
        MINC         = 2
        startstop    = []
        lastChainId1 = None
        lastResSeq1  = None
        count        = None
        start        = None
        stop         = None
        backwards    = 0
        for i, c in enumerate( contacts ):
            
            # Assign the secondary structure
            ss = dsspP.getAssignment( c['resSeq1'], c['chainId1'], resName = c['aa1'] )
            
            #print "DATA: ",i, c['chainId1'], c['resSeq1'], c['aa1'], c['chainId2'], c['resSeq2'], c['aa2'], ss
            if i != 0:
                # For getting the longest segment we only care it's going up  by 1 and it's helix - we don't care what matches how
                if ( c['resSeq1'] == lastResSeq1 + 1 or c['resSeq1'] == lastResSeq1 - 1 ) and c['chainId1'] == lastChainId1 and ss =='H': # Native is incrementing
                    
                    if c['resSeq1'] == lastResSeq1 - 1:
                        backwards += 1
                    
                    if  not( c['resSeq1'] == lastResSeq1 + 1 and backwards >= 1 ):
                        count += 1
                        stop = c['resSeq1']
                        lastChainId1 = c['chainId1']
                        lastResSeq1 = c['resSeq1']
                                
                        # If this is the last one we want to drop through
                        if i < len( contacts ) - 1:
                            continue
                    
                # Anything that doesn't continue didn't match
                if count >= MINC:
                    startstop.append( ( lastChainId1, start, stop )  )
            
            # Either starting afresh or a random residue
            lastChainId1 = c['chainId1']
            lastResSeq1  = c['resSeq1']
            backwards = 0
            if ss == 'H':
                # Only count this one if its a helix
                count = 1
                start = c['resSeq1']
                stop  = c['resSeq1']
            else:
                count = 0
                start = c['resSeq1'] + 1
                stop  = c['resSeq1'] + 1
                
        # END LOOP
        
        if not len( startstop ):
            return None
        
        def join_chunks( chunk1, chunk2 ):
            """Take two chunks of contacts and see if the gap between them can be joined.
            
            If the chunks can't be joined, it returns the first chunk for adding to the list
            of chunks, and the second chunk for use in the subsequent join step.
            If the chunks can be joined, it returns None for the first chunk so that we know not
            to add anything.
            """
            
            #print "c1, c2",chunk1, chunk2
            MAXGAP = 3
            
            chainId1 = chunk1[0]
            start1   = chunk1[1]
            stop1    = chunk1[2]
            chainId2 = chunk2[0]
            start2   = chunk2[1]
            stop2    = chunk2[2]
            
            #print "GAP PARAM ",start1,stop1,chainId1,start2,stop2,chainId2
            #print "GAP WIDTH ",start2 - stop1 + 1
            
            # See if a suitable gap
            width = start2 - stop1 + 1
            if width < 1 or width > MAXGAP or chainId1 != chainId2:
                return ( chunk1, chunk2 )
            
            #print "POSSIBLE GAP"
                
            # Suitable gap so make sure it's all helix
            for resSeq in range( stop1 + 1, start2 - 1 ):
                ss = dsspP.getAssignment( resSeq, chainId1, resName = None )
                if ss != 'H':
                    return ( chunk1, chunk2 )
                
            
            #print "JOINED CHUNKS"
                
            return ( None, ( chainId1, start1, stop2 ) )
        
        #
        # Go through the start-stop chunks in pairs and see if they can be joied, creating
        # extended, which is the list of chunks with gaps filled-in
        #
        #print "GOT STARTSTOP ",startstop
        if len( startstop ) > 1:
            
            # Need to sort the chunks by chain and then start index.
            startstop.sort( key = itemgetter( 0, 1 ) ) # By chain
            #print "SORTED STARTSTOP ",startstop
            
            extended = []
            for i, newChunk in enumerate( startstop ):
                
                # initialise
                if i == 0:
                    toJoin = newChunk
                    continue
                
                chunk, toJoin = join_chunks( toJoin, newChunk  )
                
                if chunk is not None:
                    extended.append( chunk )
                
                # Last one needs to be handled specially
                if i == len( startstop ) - 1 and toJoin:
                    extended.append( toJoin )
                    
            # End Loop
            
            #print "GOT EXTENDED ",extended
            
            #
            # Find the biggest
            #
            biggest = sorted( extended, lambda x, y: abs( x[2] - x[1]) - abs(y[2] - y[1]), reverse = True )[0]
        
        else:
            biggest = startstop[ 0 ]
            
        #print "GOT BIGGEST ",biggest
        
        #
        # Get the sequence that the start, stop indices define
        #
        chainId     = biggest[0]
        startResSeq = min( biggest[1], biggest[2] ) # use min/max as could be running backwards
        stopResSeq  = max( biggest[1], biggest[2] )
        sequence = ""
        #s3 = []
        #print " s ",startResSeq
        #print " e ",stopResSeq
        for resSeq in range( startResSeq, stopResSeq + 1):
            resName  = dsspP.getResName( resSeq, chainId )
            sequence += resName
#             try:
#                 s3 .append( pdb_edit.one2three[ resName ] )
#             except KeyError:
#                 s3.append( "XXX" )
        #print "s3 "," ".join( s3 )
        
        #print "HELIX ",sequence
                
        return sequence
        
    def run( self, nativePdb=None, placedPdb=None, resSeqMap=None, nativeInfo=None, shelxePdb=None, workdir=None ):
        """
        """

        self.workdir = workdir
        if not self.workdir:
            self.workdir = os.getcwd()
            
        pdbedit = pdb_edit.PDBEdit()
        
        # Object to hold contact data
        self.best = ContactData()
        
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
        self.best.floatingOrigin = floating
        
        # Add the shelxe origin to the list if it's not already in there
        csym = csymmatch.Csymmatch()
        corig = None
        if shelxePdb:
            csymmatchPdb = ample_util.filename_append( filename=shelxePdb, astr="csymmatch", directory=self.workdir )
            csym.run( refPdb=nativePdb, inPdb=shelxePdb, outPdb=csymmatchPdb )
            corig = csym.origin()
            
        self.best.csymmatchOrigin = bool( corig )
        #if not corig:
        #    print "NO CSYMMATCH ORIGIN"
        
        # For floating origins we use the csymmatch origin
        if floating:
            if not corig:
                # If csymmatch failed, we can't do owt
                print "CSYMMATCH FAILED WITH FLOATING ORIGIN"
                return False
            # Should check if the origin is acceptable, but that would require checking through all the 
            # alternate origins and seeing if tne non-floating axes had acceptable values  
            origins = [ corig ]
        else:
            if corig and corig not in origins:
                #print "csymmatch origin {0} is not in origins {1}".format( corig, origins )
                origins.append( corig )
        
        # Loop over origins, move the placed pdb to the new origin and then run ncont
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
            
            # Save the result if there are more 'good' contacts than anything else
            if self.inregister + self.ooregister > self.best.inregister + self.best.ooregister:
                self.best.numContacts = self.numContacts
                self.best.inregister = self.inregister
                self.best.ooregister = self.ooregister
                self.best.backwards = self.backwards
                self.best.contacts = self.contacts
                self.best.origin = origin
                self.best.ncontLog = self.ncontLog
                self.best.pdb = placedOriginPdb
                
            #print "GOT CONTACTS: {0} : {1} : {2}".format( self.numContacts, self.inregister, self.ooregister  )
        
        if self.best.pdb:
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
            d = {}
            d['chainId1'] = c[4]
            d['resSeq1']  = int( c[6:10].strip() )
            aa1           = c[11:14]
            d['aa1']      = pdb_edit.three2one[ aa1 ] # get amino acid and convert to single letter code
            d['chainId2'] = c[32]
            d['resSeq2']  = int( c[34:38].strip() )
            aa2           = c[39:42]
            d['aa2']      = pdb_edit.three2one[ aa2 ]
            d['dist']     = float( c[56:62].strip() )
            d['cell']     = int( c[63:66])
            d['symmetry'] = c[67:]
            
            contacts.append( d )
            #contacts.append( (chainId1, resSeq1, aa1, chainId2, resSeq2, aa2, dist, cell, symmetry ) )

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

        #
        # Loop through the contacts finding the start, stop indices in the list of contacts of contiguous chunks
        #
        MINC         = 3
        startstop    = []
        lastChainId1 = None
        lastChainId2 = None
        lastResSeq1  = None
        lastResSeq2  = None
        count        = None
        start        = None
        stop         = None
        contacts = self.contacts
        backwards = 0
        for i, c in enumerate( contacts ):
            
            # Unpack contacts
            #chainId1, resSeq1, aa1, chainId2, resSeq2, aa2, dist, cell, symmetry = c
            
            #print "DATA: ",i, c['chainId1'], c['resSeq1'], c['aa1'], c['chainId2'], c['resSeq2'], c['aa2']
            if i != 0:
                # For getting the longest segment we only care it's changing by 1 - we test what matched how later
                if ( c['resSeq1'] == lastResSeq1 + 1 or c['resSeq1'] == lastResSeq1 - 1 ) and \
                   ( c['resSeq2'] == lastResSeq2 + 1 or c['resSeq2'] == lastResSeq2 - 1 ) and \
                     c['chainId1'] == lastChainId1 and c['chainId2'] == lastChainId2: 
                    
                    # Need to know when we are going backwards, but also need to know if we are in a stretch going backwards
                    # so we count how many we have been going backwards for
                    if c['resSeq1'] == lastResSeq1 - 1:
                        backwards += 1
                        
                    if  not( c['resSeq1'] == lastResSeq1 + 1 and backwards >= 1 ):
                        count += 1
                        stop = i
                        lastChainId1 = c['chainId1']
                        lastChainId2 = c['chainId2']
                        lastResSeq1 = c['resSeq1']
                        lastResSeq2 = c['resSeq2']
                            
                        # If this is the last one we want to drop through
                        if i < len( contacts ) - 1:
                            continue
                    
                # Anything that doesn't continue didn't match
                if count >= MINC:
                    startstop.append( ( start, stop )  )
            
            # Either starting afresh or a random residue
            lastChainId1 = c['chainId1']
            lastChainId2 = c['chainId2']
            lastResSeq1  = c['resSeq1']
            lastResSeq2  = c['resSeq2']
            count        = 1
            start        = i
            stop         = i
            backwards    = 0
            
        # END LOOP
        
        #print "GOT STARTSTOP ",startstop
        
        # Now categorise the chunks
        for ( start, stop ) in startstop:
            #print "NEW CHUNK"
            register    = True
            backwards   = False
            lastResSeq1 = None
            lastResSeq2 = None
            count       = 1
            #print "LEN ", stop - start + 1
            for i in range( start, stop + 1 ):
                
                c = contacts[ i ]
                #print "CONTACT: ",i, c['chainId1'], c['resSeq1'], c['aa1'], c['chainId2'], c['resSeq2'], c['aa2']
                
                if i != start:
                    if c['resSeq1'] == lastResSeq1 - 1 or c['resSeq2'] == lastResSeq2 - 1:
                        backwards = True
                        
                    # Sanity check
                    if backwards:
                        assert c['resSeq1'] == lastResSeq1 -1 or  c['resSeq2'] == lastResSeq2 - 1,"Going backwards but resSeq arent!"
                    if not register and not backwards:
                        assert c['resSeq1'] != c['resSeq2'],"not register and backwards and not matching"
                        
                    count += 1
                
                lastResSeq1  = c['resSeq1']
                lastResSeq2  = c['resSeq2']
                if c['resSeq1']  != c['resSeq2']:
                    register = False
            
            # Add total here
            if register:
                self.inregister += count
            else:
                self.ooregister += count
            
            if backwards:
                self.backwards += count
        
        # End categorising chunks
            
        return
    
    def writeHelixFile(self, filename=None ):
        if self.best.helix:
            with open( filename, 'w' ) as f:
                f.write( self.best.helix+"\n" )
            return True
        
        return False
    
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
    
    def testParse5(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont5.log" )
        
        c = Contacts()
        contacts = c.parseNcontLog( logfile=logfile )
        c.countContacts()
        
        self.assertEqual( c.numContacts, 77 )
        self.assertEqual( c.inregister, 19 )
        self.assertEqual( c.backwards, 16 )
        self.assertEqual( c.ooregister, 54 )
        
        return

    def testParse7(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont7.log" )
        
        c = Contacts()
        contacts = c.parseNcontLog( logfile=logfile )
        c.countContacts()
        
        self.assertEqual( c.numContacts, 18 )
        self.assertEqual( c.inregister, 0 )
        self.assertEqual( c.backwards, 0 )
        self.assertEqual( c.ooregister, 0 )
        
        return
    
    def testParse8(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont8.log" )
        
        c = Contacts()
        contacts = c.parseNcontLog( logfile=logfile )
        c.countContacts()
        
        self.assertEqual( c.numContacts, 9 )
        self.assertEqual( c.inregister, 0 )
        self.assertEqual( c.backwards, 0 )
        self.assertEqual( c.ooregister, 0 )
        
        return


    def testHelix5(self):

        logfile = os.path.join( self.testfilesDir, "ncont5.log" )
        
        import dssp
        dssplog = os.path.join( self.testfilesDir, "3RA3.dssp" )
        
        dsspP = dssp.DsspParser( dssplog )
        
        c = Contacts()
        contacts = c.parseNcontLog( logfile=logfile )
        sequence = c.helixFromContacts( contacts=contacts, dsspP=dsspP )
        
        self.assertEqual( "NARLKQEIAALEYEIAAL", sequence )
        
        return

if __name__ == "__main__":
    
    unittest.main()
    
    import sys
    sys.exit()


