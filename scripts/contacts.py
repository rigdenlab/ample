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


class ContactData(object):
    def __init__(self):

        
        # Hack to transfer data
        self.joinedPdb = None
        self.fromChains = None
        self.toChains = None
        
        self.origin = None
        self.originPdb = None
        
        self.numContacts = None # Used by both as parseNcontLog sets it
        self.aaNumContacts = 0
        self.rioNumContacts = 0
        self.rioInRegister  = 0
        self.rioOoRegister  = 0
        self.rioBackwards   = 0
        self.numGood     = 0
        
        self.contacts = None
        
        self.ncontLog = None
        self.helix = None
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
        self.best = None
        return
    
#     def getContacts( self, 
#                      nativePdbInfo=None,
#                      placedPdbInfo=None,
#                      resSeqMap=None,
#                      origins=None, 
#                      workdir=None,
#                      dsspLog=None ):
#         
#         if not self.run( nativePdbInfo=nativePdbInfo,
#                          placedPdbInfo=placedPdbInfo,
#                          resSeqMap=resSeqMap,
#                          origins=origins, 
#                          workdir=workdir ):
#             return False
#         
#         
#         #print "placedPdb ",placedPdb
#         #print "nativePdb ",nativePdb
#         #print "GOT MAP ",resSeqMap
#         
#         # We should have contact data to calculate a helix from
#         assert dsspLog
#         dsspP = dssp.DsspParser( dsspLog )
#         sequence = self.helixFromContacts( dsspP=dsspP )
#         if sequence:
#             self.best.helix = sequence
#         
#         if self.best.pdb:
#             csym = csymmatch.Csymmatch()
#             # Just for info - run csymmatch so we can see the alignment
#             csymmatchPdb = ample_util.filename_append( filename=self.best.pdb, astr="csymmatch_best", directory=self.workdir )
#             csym.run( refPdb=nativePdbInfo.pdb, inPdb=self.best.pdb, outPdb=csymmatchPdb, originHand=False )
#             self.best.csymmatchPdb = csymmatchPdb
#         
#         return

    def analyseRio( self, contactData, dsspP=None ):
        
        # Clear any data
        contactData.rioInRegister = 0
        contactData.rioOoRegister = 0
        contactData.rioBackwards = 0
        
        if contactData.numContacts == 0:
            return

        # Now categorise the chunks
        chunks = self.findChunks( contacts=contactData.contacts, minContig=3 )
        if not chunks:
            return 0 
        
        #print "GOT CHUNKS ",chunks
        
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
                
                c = contactData.contacts[ i ]
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
                contactData.rioInRegister += count
            else:
                contactData.rioOoRegister += count
            
            if backwards:
                contactData.rioBackwards += count
        
        # End categorising chunks
        return

    def calcAllAtom(self, contactData ):
        """Calculate the All Atom data using the origin and data saved previously"""

        # Clear out any set data
        contactData.numContacts   = 0
        contactData.aaNumContacts = 0
        contactData.contacts      = None
        
        # values from previous run to calc origin
        self.runNcont( pdbin=contactData.joinedPdb,
                       sourceChains=contactData.fromChains,
                       targetChains=contactData.toChains,
                       allAtom=True,
                       maxDist=0.5 )
        self.parseNcontLog( contactData )
        
        contactData.aaNumContacts = contactData.numContacts
        
        
        return contactData
    
    def calcRio(self, contactData ):
        """Calculate the RIO using the origin and data saved previously"""

        # Clear out any set data
        #contactData.numGood        = 0
        contactData.rioNumContacts = 0
        contactData.rioInRegister  = 0
        contactData.rioOoRegister  = 0
        contactData.rioBackwards   = 0
        contactData.contacts       = None
        
        # values from previous run to calc origin
        self.runNcont( pdbin=contactData.joinedPdb,
                       sourceChains=contactData.fromChains,
                       targetChains=contactData.toChains,
                        )
        self.parseNcontLog( contactData )
        self.analyseRio( contactData )
        
        contactData.rioNumContacts = contactData.numContacts
        
        return contactData
        
    def findOrigin(self,
                   nativePdbInfo=None,
                   mrPdbInfo=None,
                   resSeqMap=None,
                   origins=None,
                   allAtom=False,
                   workdir=os.getcwd() ):
        """Find the origin using the maximum number of contacts as metric"""
        
        self.workdir = workdir
        pdbedit = pdb_edit.PDBEdit()

        if not resSeqMap.resSeqMatch():
            #print "NUMBERING DOESN'T MATCH"
            #raise RuntimeError,"NUMBERING DOESN'T MATCH"
            # We need to create a copy of the placed pdb with numbering matching the native
            mrPdbRes = ample_util.filename_append( filename=mrPdbInfo.pdb, astr="reseq", directory=self.workdir )
            pdbedit.match_resseq( targetPdb=mrPdbInfo.pdb, sourcePdb=None, outPdb=mrPdbRes, resMap=resSeqMap )
            mrPdb = mrPdbRes
        else:
            mrPdb = mrPdbInfo.pdb
 
        # Make a copy of mrPdb with chains renamed to lower case
        ucChains = mrPdbInfo.models[0].chains
        toChains = [ c.lower() for c in ucChains ]
        placedAaPdb = ample_util.filename_append( filename=mrPdb, astr="ren", directory=self.workdir )
        pdbedit.rename_chains( inpdb=mrPdb, outpdb=placedAaPdb, fromChain=ucChains, toChain=toChains )

        # The list of chains in the native that we will be checking contacts from
        fromChains = nativePdbInfo.models[0].chains
        
        # Loop over origins, move the placed pdb to the new origin and then run ncont
        # Object to hold data on best origin
        self.data = None
        for origin in origins:
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
            
            # Set up object to hold data
            data            = ContactData()
            data.origin     = origin
            data.originPdb  = placedOriginPdb
            data.joinedPdb  = joinedPdb
            data.fromChains = fromChains
            data.toChains   = toChains
            data.numGood    = 0 # For holding the metric
            
            # Run ncont
            if allAtom:
                self.calcAllAtom( data )
                data.numGood = data.aaNumContacts
            else:
                self.calcRio( data )
                data.numGood = data.rioInRegister + data.rioOoRegister
            
            # Save the first origin and only update if we get a better score
            if not self.data or data.numGood > self.data.numGood:
                self.data = data

        # End loop over origins
        
        # Now need to calculate data for whichever one we didn't calculate
        if allAtom:
            self.calcRio( self.data )
        else:
            self.calcAllAtom( self.data )
        
        if self.data.numGood > 0:
            
            # If we got a match run csymmatch so we can see the result
            csym = csymmatch.Csymmatch()
            csymmatchPdb = ample_util.filename_append( filename=self.data.originPdb,
                                                       astr="csymmatch_best",
                                                       directory=self.workdir )
            csym.run( refPdb=nativePdbInfo.pdb,
                      inPdb=self.data.originPdb,
                      outPdb=csymmatchPdb,
                      originHand=False )
        
        return self.data

    def helixFromContacts( self, contacts, dsspLog, minContig=2, maxGap=3  ):
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
        if not len( contacts ):
            return None
        
        #print "GOT DSSP ",dsspP.asDict()
        # Parse the dssp Log
        dsspP = dssp.DsspParser( dsspLog )
        #
        # Loop through the contacts finding the start, stop indices in the list of contacts of contiguous chunks
        #
        chunks = self.findChunks( contacts=contacts, dsspP=dsspP, ssTest=True, minContig=2 )
        if not chunks:
            return None
        #
        # Go through the start-stop chunks in pairs and see if they can be joined, creating
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
                
                chunk, toJoin = self._join_chunks( toJoin, newChunk, dsspP=dsspP, maxGap=maxGap  )
                
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

    def _join_chunks( self, chunk1, chunk2, dsspP=None, maxGap=None ):
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

    def runNcont( self, pdbin=None, sourceChains=None, targetChains=None, maxDist=1.5, allAtom=False ):
        """FOO
        """
        
        if allAtom:
            self.ncontLog = pdbin + ".ncont_aa.log"
        else:
            self.ncontLog = pdbin + ".ncont_rio.log"
        
        cmd = [ "ncont", "xyzin", pdbin ]
        
        # Build up stdin
        stdin = ""
        # Need to use list of chains from Native as can't work out negate operator for ncont
        if allAtom:
                stdin += "source {0}//*\n".format( ",".join( sourceChains )  )  
                stdin += "target {0}//*\n".format( ",".join( targetChains )  ) 
        else:
            stdin += "source {0}//CA\n".format( ",".join( sourceChains )  )  
            stdin += "target {0}//CA\n".format( ",".join( targetChains )  )  
        stdin += "maxdist {0}\n".format( maxDist )
        stdin += "cells 2\n"
        stdin += "sort target inc\n"
        
        retcode = ample_util.run_command( cmd=cmd, logfile=self.ncontLog, directory=os.getcwd(), dolog=False, stdin=stdin )
        
        if retcode != 0:
            raise RuntimeError,"Error running ncont"
        
        return
    
    def parseNcontLog( self, contactData, logfile=None ):
        """
        
        Lines are of format
        /1/B/1042(MET). / CA [ C]:  /1/b/ 988(GLU). / CA [ C]:   1.09 223 X-1/2,Y-1/2,Z
 
        """
        
        if not logfile:
            logfile = self.ncontLog
        
        #print "LOG ",logfile
        
        
        contactData.contacts = None
        contactData.numContacts = 0
        clines = []
        
        capture=False
        with open( logfile, 'r' ) as f:#
            while True:
                line = f.readline().rstrip()
                
                if capture and not line:
                    break
                
                if "contacts found:" in line:
                    contactData.numContacts = int( line.split()[0] )
                
                if "NO CONTACTS FOUND." in line:
                    return False
                
                if "SOURCE ATOMS" in line:
                    capture=True
                    f.readline() # skip blank line
                    continue
                
                if capture:
                    clines.append( line )
            
        assert contactData.numContacts == len(clines)
        #print "LINES ",clines

        contacts = [] 
        # Got data lines so now extract data
        # Could probably just do this in the reading loop now    
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
            # Will get key error for all-atom if there is solvent etc.
            try:
                d['aa1']      = pdb_edit.three2one[ aa1 ] # get amino acid and convert to single letter code
            except KeyError:
                d['aa1']      = 'X'
            d['chainId2'] = c[32]
            d['resSeq2']  = int( c[34:38].strip() )
            aa2           = c[39:42]
            try:
                d['aa2']      = pdb_edit.three2one[ aa2 ]
            except KeyError:
                d['aa2']      = 'X'
            d['dist']     = float( c[56:62].strip() )
            d['cell']     = int( c[63:66])
            d['symmetry'] = c[67:]
            
            contacts.append( d )
    
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
    
            #print "GOT CONTACTS"
            #for c in contacts:
            #    print "chainId1 {0} resSeq1 {1} chainId2 {2} resSeq2 {3}\n".format(  c['chainId1'], c['resSeq1'], c['chainId2'], c['resSeq2']  )
    
        contactData.contacts = contacts
                    
        return
    
    def scoreOrigin(self,
                    origin=None,
                    mrPdbInfo=None,
                    nativePdbInfo=None,
                    resSeqMap=None,
                    workdir=os.getcwd()
                     ):
        
        self.workdir = workdir
        pdbedit = pdb_edit.PDBEdit()

        if not resSeqMap.resSeqMatch():
            #print "NUMBERING DOESN'T MATCH"
            #raise RuntimeError,"NUMBERING DOESN'T MATCH"
            # We need to create a copy of the placed pdb with numbering matching the native
            mrPdbRes = ample_util.filename_append( filename=mrPdbInfo.pdb, astr="reseq", directory=self.workdir )
            pdbedit.match_resseq( targetPdb=mrPdbInfo.pdb, sourcePdb=None, outPdb=mrPdbRes, resMap=resSeqMap )
            mrPdb = mrPdbRes
        else:
            mrPdb = mrPdbInfo.pdb
 
        # Make a copy of mrPdb with chains renamed to lower case
        ucChains = mrPdbInfo.models[0].chains
        toChains = [ c.lower() for c in ucChains ]
        mrAaPdb = ample_util.filename_append( filename=mrPdb, astr="ren", directory=self.workdir )
        pdbedit.rename_chains( inpdb=mrPdb, outpdb=mrAaPdb, fromChain=ucChains, toChain=toChains )

        # The list of chains in the native that we will be checking contacts from
        fromChains = nativePdbInfo.models[0].chains
        
        mrOriginPdb =  mrAaPdb
        if origin != [ 0.0, 0.0, 0.0 ]:
            # Move pdb to new origin
            #ostr="origin{0}".format(i)
            ostr="o{0}".format( origin ).replace(" ","" )
            mrOriginPdb = ample_util.filename_append( filename=mrAaPdb, astr=ostr, directory=self.workdir )
            pdbedit.translate( inpdb=mrAaPdb, outpdb=mrOriginPdb, ftranslate=origin )
        
        # Concatenate into one file
        joinedPdb = ample_util.filename_append( filename=mrOriginPdb, astr="joined", directory=self.workdir )
        pdbedit.merge( pdb1=nativePdbInfo.pdb, pdb2=mrOriginPdb, pdbout=joinedPdb )
            
        # Run ncont
        data = ContactData()
        data.origin = origin
        data.originPdb = mrOriginPdb
        data.joinedPdb = joinedPdb
        data.fromChains = fromChains
        data.toChains = toChains

        # First get AllAtom score        
        self.calcAllAtom( data )
        
        # Then score RIO
        self.calcRio( data )
        #data.numGood = data.inregister + data.ooregister
        
        return data
 
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
    
class TestContacts( unittest.TestCase ):
    
    def setUp(self):
        
        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        self.ampleDir = os.sep.join( paths[ : -1 ] )
        self.testfilesDir = os.sep.join( paths[ : -1 ] + [ 'tests', 'testfiles' ] )
        
        return

    def testParse1(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont1.log" )
        
        c = Contacts()
        c.parseNcontLog( logfile=logfile )
        c.analyseRio()
        
        self.assertEqual( c.numContacts, 26 )
        self.assertEqual( c.inregister, 0 )
        self.assertEqual( c.ooregister, 0 )
        self.assertEqual( c.backwards, 0 )
        
        return
    
    def testParse2(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont2.log" )
        
        c = Contacts()
        c.parseNcontLog( logfile=logfile )
        c.analyseRio()
        
        self.assertEqual( c.numContacts, 10 )
        self.assertEqual( c.inregister, 0 )
        self.assertEqual( c.ooregister, 7 )
        self.assertEqual( c.backwards, 7 )
        
        return
    
    def testParse3(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont3.log" )
        
        c = Contacts()
        c.parseNcontLog( logfile=logfile )
        c.analyseRio()
        
        self.assertEqual( c.numContacts, 14 )
        self.assertEqual( c.inregister, 0 )
        self.assertEqual( c.ooregister, 10 )
        self.assertEqual( c.backwards, 0 )
        
        return
    
    def testParse4(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont4.log" )
        
        c = Contacts()
        c.parseNcontLog( logfile=logfile )
        c.analyseRio()
        
        self.assertEqual( c.numContacts, 56 )
        self.assertEqual( c.inregister, 0 )
        self.assertEqual( c.ooregister, 55 )
        self.assertEqual( c.backwards, 0 )
        
        return
    
    def testParse5(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont5.log" )
        
        c = Contacts()
        c.parseNcontLog( logfile=logfile )
        c.analyseRio()
        
        self.assertEqual( c.numContacts, 77 )
        self.assertEqual( c.inregister, 19 )
        self.assertEqual( c.backwards, 16 )
        self.assertEqual( c.ooregister, 54 )
        
        return

    def testParse7(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont7.log" )
        
        c = Contacts()
        c.parseNcontLog( logfile=logfile )
        c.analyseRio()
        
        self.assertEqual( c.numContacts, 18 )
        self.assertEqual( c.inregister, 0 )
        self.assertEqual( c.backwards, 0 )
        self.assertEqual( c.ooregister, 0 )
        
        return
    
    def testParse8(self):
        
        logfile = os.path.join( self.testfilesDir, "ncont8.log" )
        
        c = Contacts()
        c.parseNcontLog( logfile=logfile )
        c.analyseRio()
        
        self.assertEqual( c.numContacts, 9 )
        self.assertEqual( c.inregister, 0 )
        self.assertEqual( c.backwards, 0 )
        self.assertEqual( c.ooregister, 0 )
        
        return


    def testHelix5(self):

        logfile = os.path.join( self.testfilesDir, "ncont5.log" )
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


