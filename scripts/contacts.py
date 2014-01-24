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
    
    def getContacts( self, 
                     nativePdbInfo=None,
                     placedPdbInfo=None,
                     resSeqMap=None,
                     origins=None, 
                     workdir=None,
                     dsspLog=None ):
        
        if not self.run( nativePdbInfo=nativePdbInfo,
                         placedPdbInfo=placedPdbInfo,
                         resSeqMap=resSeqMap,
                         origins=origins, 
                         workdir=workdir ):
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
            csym.run( refPdb=nativePdbInfo.pdb, inPdb=self.best.pdb, outPdb=csymmatchPdb, originHand=False )
            self.best.csymmatchPdb = csymmatchPdb
        
        return

    def helixFromContacts( self, contacts=None, dsspP=None, minContig=2, maxGap=3  ):
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
            
        if contacts is None or not len( contacts ):
            return None
            
        #
        # Loop through the contacts finding the start, stop indices in the list of contacts of contiguous chunks
        #
        chunks = self.findChunks( contacts=contacts, dsspP=dsspP, ssTest=True, minContig=2 )
    
        if not chunks:
            return None
        
        def join_chunks( chunk1, chunk2, dsspP=None, maxGap=None ):
            """Take two chunks of contacts and see if the gap between them can be joined.
            
            If the chunks can't be joined, it returns the first chunk for adding to the list
            of chunks, and the second chunk for use in the subsequent join step.
            If the chunks can be joined, it returns None for the first chunk so that we know not
            to add anything.
            """
            
            assert dsspP and maxGap
            
            #print "CHECKING CHUNK ",chunk1, chunk2
            
            # See if a suitable gap
            width = abs( chunk2[ 'startResSeq' ] - chunk1[ 'stopResSeq' ] ) - 1
            if width < 1 or width > maxGap or chunk1[ 'chainId1' ] != chunk2[ 'chainId1' ]:
                return ( chunk1, chunk2 )
            
            # Suitable gap so make sure it's all helix
            for resSeq in range( chunk1[ 'stopResSeq' ] + 1, chunk2[ 'startResSeq' ] - 1 ):
                ss = dsspP.getAssignment( resSeq, chunk1[ 'chainId1' ], resName = None )
                if ss != 'H':
                    return ( chunk1, chunk2 )
            
            #print "JOINED CHUNKS", chunk1, chunk2
            
            joinedChunk =  { 'chainId1'    : chunk1[ 'chainId1' ],
                             'startIdx'    : chunk1[ 'startIdx'],
                             'startResSeq' : chunk1[ 'startResSeq' ],
                             'stopIdx'     : chunk2[ 'stopIdx'],
                             'stopResSeq'  : chunk2[ 'stopResSeq' ] }
                
            return ( None, joinedChunk )
        
        #
        # Go through the start-stop chunks in pairs and see if they can be joied, creating
        # extended, which is the list of chunks with gaps filled-in
        #
        #print "GOT CHUNKS ",chunks
        if len( chunks ) > 1:
            
            # Need to sort the chunks by chain and then startResSeq so that we can join anything on the same chain
            chunks.sort( key = itemgetter( 'chainId1', 'startResSeq' ) ) # By chain
            #print "SORTED CHUNKS ",chunks
            
            extended = []
            for i, newChunk in enumerate( chunks ):
                
                # initialise
                if i == 0:
                    toJoin = newChunk
                    continue
                
                chunk, toJoin = join_chunks( toJoin, newChunk, dsspP=dsspP, maxGap=maxGap  )
                
                if chunk is not None:
                    extended.append( chunk )
                
                # Last one needs to be handled specially
                if i == len( chunks ) - 1 and toJoin:
                    extended.append( toJoin )
                    
            # End Loop
            
            #print "GOT EXTENDED ",extended
            
            #
            # Find the biggest
            #
            biggest = sorted( extended, lambda x, y: abs( x['stopResSeq'] - x['startResSeq']) - abs(y['stopResSeq'] - y['startResSeq']), reverse = True )[0]
        
        else:
            biggest = chunks[ 0 ]
            
        #print "GOT BIGGEST ",biggest
        
        #
        # Get the sequence that the start, stop indices define
        #
        chainId     = biggest[ 'chainId1']
        startResSeq = min( biggest['startResSeq'], biggest['stopResSeq'] ) # use min/max as could be running backwards
        stopResSeq  = max( biggest['startResSeq'], biggest['stopResSeq'] )
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
        
    def run( self,
             nativePdbInfo=None,
             placedPdbInfo=None,
             resSeqMap=None,
             origins=None,
             workdir=None ):
        """
        """

        # Reset variables each time
        self.ncontLog = None
        self.numContacts = 0
        self.contacts = None
        self.inregister = 0
        self.ooregister = 0
        self.backwards = 0
        # Object to hold contact data
        self.best = ContactData()

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
            #modelInfo = pdbedit.get_info( refModelPdb )
            modelInfo = False
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
            placedPdbRes = ample_util.filename_append( filename=placedPdbInfo.pdb, astr="reseq", directory=self.workdir )
            pdbedit.match_resseq( targetPdb=placedPdbInfo.pdb, sourcePdb=None, outPdb=placedPdbRes, resMap=resSeqMap )
            placedPdb = placedPdbRes
        else:
            placedPdb = placedPdbInfo.pdb
 
        # Make a copy of placedPdb with chains renamed to lower case
        fromChain = placedPdbInfo.models[0].chains
        toChain = [ c.lower() for c in fromChain ]
        placedAaPdb = ample_util.filename_append( filename=placedPdb, astr="ren", directory=self.workdir )
        pdbedit.rename_chains( inpdb=placedPdb, outpdb=placedAaPdb, fromChain=fromChain, toChain=toChain )

        # Don't test as the labels could be differnt and it's not worth creating anther one just to test
        #placedSpaceGroup = placedInfo.crystalInfo.spaceGroup
        #if placedSpaceGroup != originInfo.currentSpaceGroup():
        #    raise RuntimeError,"Mismatching space groups!"
        
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
            pdbedit.merge( pdb1=nativePdbInfo.pdb, pdb2=placedOriginPdb, pdbout=joinedPdb )
                
            # Run ncont
            # Need to get list of chains from Native as can't work out negate operator for ncont
            fromChain = nativePdbInfo.models[0].chains
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
 
    def ssIsOK( self, contact, dsspP=None, ssTest=False ):
        if ssTest:
            ss = dsspP.getAssignment( contact['resSeq1'], contact['chainId1'], resName = contact['aa1'] )
            if ss == "H":
                return True
            return False
        else:
            return True

    def findChunks(self, contacts=None, dsspP=None, ssTest=False, minContig=3 ):
        
        if contacts is None:
            return False

        #
        # Loop through the contacts finding the start, stop indices in the list of contacts of contiguous chunks
        #
        chunks    = []
        lastChainId1 = None
        lastChainId2 = None
        lastResSeq1  = None
        lastResSeq2  = None
        count        = None
        startIdx     = None
        stopIdx      = None
        backwards = 0
        
        for i, c in enumerate( contacts ):
            
            ssOK = self.ssIsOK( c, dsspP=dsspP, ssTest=ssTest )
            
            #ss = dsspP.getAssignment( c['resSeq1'], c['chainId1'], resName = c['aa1'] )
            #print "DATA: ",i, c['chainId1'], c['resSeq1'], c['aa1'], c['chainId2'], c['resSeq2'], c['aa2'],ss, ssOK
            
            if i != 0:
                # For getting the longest segment we only care it's changing by 1 - we test what matched how later
                if ( c['resSeq1'] == lastResSeq1 + 1 or c['resSeq1'] == lastResSeq1 - 1 ) and \
                   ( c['resSeq2'] == lastResSeq2 + 1 or c['resSeq2'] == lastResSeq2 - 1 ) and \
                     c['chainId1'] == lastChainId1 and c['chainId2'] == lastChainId2 and ssOK: 
                    
                    # Need to know when we are going backwards, but also need to know if we are in a stretch going backwards
                    # So we can spot if there has been a change, so we count how many we have been going backwards for
                    if c['resSeq1'] == lastResSeq1 - 1:
                        backwards += 1
                        
                    if  not( c['resSeq1'] == lastResSeq1 + 1 and backwards >= 1 ):
                        count       += 1
                        stopIdx      = i
                        lastChainId1 = c['chainId1']
                        lastChainId2 = c['chainId2']
                        lastResSeq1  = c['resSeq1']
                        lastResSeq2  = c['resSeq2']
                            
                        # If this is the last one we want to drop through
                        if i < len( contacts ) - 1:
                            continue
                    
                # Anything that doesn't continue didn't match
                if count >= minContig:
                    d = { 'chainId1'    : lastChainId1,
                          'startIdx'    : startIdx,
                          'startResSeq' : contacts[ startIdx ][ 'resSeq1' ],
                          'stopIdx'     : stopIdx,
                          'stopResSeq'  : contacts[ stopIdx ][ 'resSeq1' ] }
                    chunks.append( d  )
            
            # Either starting afresh or a random residue
            lastChainId1 = c['chainId1']
            lastChainId2 = c['chainId2']
            lastResSeq1  = c['resSeq1']
            lastResSeq2  = c['resSeq2']
            backwards    = 0
            if ssOK:
                count        = 1
                startIdx     = i
                stopIdx      = i
            else:
                count    = 0
                startIdx = i + 1
                stopIdx  = i + 1
  
        return chunks
    
    def countContacts( self, contacts=None, dsspP=None ):
        
        
        if contacts is None:
            contacts = self.contacts
        if contacts is None:
            return
        
        chunks = self.findChunks( contacts=self.contacts, minContig=3 )
        
        if not chunks:
            return
        
        #print "GOT CHUNKS ",chunks
        
        # Now categorise the chunks
        self.inregister = 0
        self.ooregister = 0
        self.backwards = 0
        
        #for ( start, stop ) in startstop:
        for d in chunks:
            start = d['startIdx']
            stop  = d['stopIdx']
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


