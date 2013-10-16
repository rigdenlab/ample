'''
Useful manipulations on PDB files
'''

# Python imports
import os
import sys
import types
import unittest

# our imports
import ample_util
import pdb_edit
import pdb_model

class residueSequenceMap( object ):
    """Class for handling mapping between model and native residue indices.
    
     """
    
    def __init__( self, nativePdb=None, modelPdb=None ):
        
        self.modelResSeq = []
        self.modelSequence = None
        self.modelCAlphaMask = []
        self.modelBbMask = []
        self.modelOffset = None # Where the matched part of the sequences starts in the model
        self._modelIncomparable = None # List of atoms in the model that cannot be compared to the native
        
        self.nativeResSeq = []
        self.nativeSequence = None
        self.nativeCAlphaMask = []
        self.nativeBbMask = []
        self.nativeOffset = None
        self._nativeIncomparable = None # List of atoms in the model that cannot be compared to the model
        
        self.lenMatch = None # The number of residues that match
        
        # Like this just for testing
        if nativePdb and modelPdb:
            self.calc_map(nativePdb, modelPdb)
        
        return
    
    def native2model( self, nativeResSeq, extrapolate=True ):
        """Return the model resSeq for the given native resSeq
        If extrapolate is True, we calculate a residue number for the missing residue.
        """
        
        # Find where this residue is in the native
        native_idx = self.nativeResSeq.index( nativeResSeq )
        
        # See how far into the matching region we are
        match_idx = native_idx - self.nativeOffset
        
        if match_idx < 0:
            # We are before the matching region - see if there are valid model residues
            i = self.modelOffset + match_idx
            if i >= 0:
                return self.modelResSeq[ i ]
            
            # Need to calcualte a new index
            first = self.modelResSeq[ 0 ]
            new = first + i
            return new
        
        elif match_idx > self.lenMatch - 1:
            beyond = match_idx + 1 - self.lenMatch
            extra = len(self.modelResSeq) - ( self.modelOffset + self.lenMatch )
            if beyond <= extra:
                return self.modelResSeq[ self.modelOffset + self.lenMatch + beyond ]
            else:
                # Need to calculate a new index
                last = self.modelResSeq[ -1 ]
                toadd = beyond - extra
                return last + toadd
        else:
            # Residue is within matched region
            model_idx = self.modelOffset + match_idx
            return self.modelResSeq[ model_idx ]
        
        #return self._modelResSeqMap[ self._nativeResSeqMap.index( nativeResSeq ) ]
    
    def model2native( self, modelResSeq):
        """Return the native resSeq for the given native resSeq"""
        
        # Find where this residue is in the native
        native_idx = self.modelResSeq.index( modelResSeq )
        
        # See how far into the matching region we are
        match_idx = native_idx - self.modelOffset
        
        if match_idx < 0:
            # We are before the matching region - see if there are valid model residues
            i = self.nativeOffset + match_idx
            if i >= 0:
                return self.nativeResSeq[ i ]
            
            # Need to calcualte a new index
            first = self.nativeResSeq[ 0 ]
            new = first + i
            return new
        
        elif match_idx > self.lenMatch - 1:
            beyond = match_idx + 1 - self.lenMatch
            # the number of extra atoms present in the native beyond the end of the match
            extra = len(self.nativeResSeq) - ( self.nativeOffset + self.lenMatch )
            if beyond <= extra:
                return self.nativeResSeq[ self.nativeOffset + self.lenMatch + beyond ]
            else:
                # Need to calculate a new index
                last = self.nativeResSeq[ -1 ]
                toadd = beyond - extra
                return last + toadd
        else:
            # Residue is within matched region
            model_idx = self.nativeOffset + match_idx
            return self.nativeResSeq[ model_idx ]
        
        #return self._nativeResSeqMap[ self._modelResSeqMap.index( modelResSeq ) ]

    def modelIncomparable(self):
        """Return a list of the indices of any residues that are in the model but not in the native or
        residues for which there are no c-Alpha atoms in one or other of the models"""
        
        if self._modelIncomparable == None:
            self._modelIncomparable = []
            for i, resSeq in enumerate( self.modelResSeq ):
                
                # Before the start of the matching region
                if i < self.modelOffset:
                    self._modelIncomparable.append( resSeq )
                    continue
                
                # In matching region but no C-alpha
                if self.modelCAlphaMask[ i ]:
                    self._modelIncomparable.append( resSeq )
                    continue
                
                # Matching residues in native
                nativeResSeq = self.model2native( resSeq )
                try:
                    j = self.nativeResSeq.index( nativeResSeq )
                except ValueError:
                    # A residue that isn't actually in the native
                    self._modelIncomparable.append( resSeq )
                    continue
                
                # No C-Alpha
                if self.nativeCAlphaMask[ j ]:
                    self._modelIncomparable.append( resSeq )
                    continue
            
        return self._modelIncomparable
        #return [ resId for i, resId in enumerate( self._modelResSeqMap ) if ( self._nativeResSeqMap[ i ] == None or self.cAlphaMask[ i ] ) ]

    def nativeIncomparable(self):
        """Return a list of the indices of any residues that are in the native but not in the model or
        residues for which there are no c-Alpha atoms in one or other of the pdbs"""
        
        if self._nativeIncomparable == None:
            self._nativeIncomparable = []
            for i, resSeq in enumerate( self.nativeResSeq ):
                
                # Before the start of the matching region
                if i < self.nativeOffset:
                    self._nativeIncomparable.append( resSeq )
                    continue
                
                # In matching region but no C-alpha
                if self.nativeCAlphaMask[ i ]:
                    self._nativeIncomparable.append( resSeq )
                    continue
                
                # Matching residues in model
                modelResSeq = self.native2model( resSeq )
                try:
                    j = self.modelResSeq.index( modelResSeq )
                except ValueError:
                    # A residue that isn't actually in the native
                    self._nativeIncomparable.append( resSeq )
                    continue
                
                # No C-Alpha
                if self.modelCAlphaMask[ j ]:
                    self._nativeIncomparable.append( resSeq )
                    continue
            
        return self._nativeIncomparable
        #return [ resId for i, resId in enumerate( self._modelResSeqMap ) if ( self._nativeResSeqMap[ i ] == None or self.cAlphaMask[ i ] ) ]


    def __str__(self):
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
              isinstance(attr, types.FunctionType) ):
                me[slot] = attr
        
        s = self.__repr__() + "\n"
        for k in sorted( me.keys() ):
            s += "{0}: {1}\n".format( k, me[k]  )
        return s
    
    def calc_map( self, nativePdb, modelPdb ):
        
        self.nativeSequence, self.nativeResSeq, self.nativeCAlphaMask = self.read_pdb( nativePdb )
        self.modelSequence, self.modelResSeq, self.modelCAlphaMask = self.read_pdb( modelPdb )
        
        self._calc_map()
        
        return
    
    def _calc_map( self ):
        """Return a ResSeqMap mapping the index of a residue in the model to the corresponding residue in the native.
        Only works if 1 chain in either file and with standard residues
        """
        
        if len(self.nativeSequence) < 10 or len(self.modelSequence) < 10:
            raise RuntimeError,"Very short sequences - this will not work!"
        
        # The window of AA we used to check for a match    
        PROBE_LEN = 10
        
        # MAXINSET is the max number of AA into the sequence that we will go searching for a match - i.e. if more
        # then MAXINSET AA at the start are non-matching, we won't find the match 
        l = len( self.modelSequence ) if len( self.modelSequence ) < len( self.nativeSequence ) else len( self.nativeSequence )
        MAXINSET=30 if l > 30 else ( l - PROBE_LEN )
        
        got=False
        for modelOffset in range( MAXINSET + 1 ):
            probe = self.modelSequence[ modelOffset : modelOffset+PROBE_LEN ]
            #print "PROBE ",probe
            for nativeOffset in range( MAXINSET + 1 ):
                #print "TEST ",self.nativeSequence[ nativeOffset:nativeOffset+PROBE_LEN ]
                if self.nativeSequence[ nativeOffset:nativeOffset+PROBE_LEN ] == probe:
                    got=True
                    break
            
            if got:
#                 print "GOT MODEL MATCH AT i,j ",modelOffset,nativeOffset
                break
            
        if not got:
            raise RuntimeError,"Could not calculate map!"
        
        self.modelOffset = modelOffset
        self.nativeOffset = nativeOffset
        
        # Work out how long the match is
        count=0
        for i in range( self.modelOffset, len( self.modelSequence ) ):
            
            modelAA = self.modelSequence[ i ]
            
            y = self.nativeOffset + count
            if y >= len( self.nativeSequence ):
                break
            
            nativeAA = self.nativeSequence[ y ]
            #print "matching {0} to {1}".format( modelAA, nativeAA )
            if modelAA != nativeAA:
                break 
            
            count +=1
        
        self.lenMatch = count
            
        return
    
    def fromInfo(self, nativeInfo=None, nativeChainID=None, modelInfo=None, modelChainID=None, model=0):
        """Create a map from 2 info objects"""
        
        # Determine index of chain so we know where to get the data from
        nativeIdx = nativeInfo.models[ model ].chains.index( nativeChainID )

        self.nativeResSeq = nativeInfo.models[ model ].resSeqs[ nativeIdx ]
        self.nativeSequence = nativeInfo.models[ model ].sequences[ nativeIdx ]
        self.nativeCAlphaMask = nativeInfo.models[ model ].caMask[ nativeIdx ]
        self.nativeBbMask = nativeInfo.models[ model ].bbMask[ nativeIdx ]
        self.nativeOffset = None
        self._nativeIncomparable = None

        modelIdx = modelInfo.models[ model ].chains.index( modelChainID )
        self.modelResSeq = modelInfo.models[ model ].resSeqs[ modelIdx ]
        self.modelSequence = modelInfo.models[ model ].sequences[ modelIdx ]
        self.modelCAlphaMask = modelInfo.models[ model ].caMask[ modelIdx ]
        self.modelBbMask =  modelInfo.models[ model ].bbMask[ modelIdx ]
        self.modelOffset = None
        self._modelIncomparable = None
        
        self.lenMatch = None
        
        self._calc_map()
        
        return
    
    
    def read_pdb( self, pdb ):
        """Get sequence as string of 1AA
        get list of matching resSeq
        """
        
        atomTypes = [] # For checking we have all required atom types
    
        resSeq = []
        resName = []
        _atomTypes = []
        atomTypesList = []
        
        chain=None
        readingResSeq=None
        readingResName=None
        for line in open( pdb ):
            
            if line.startswith("MODEL"):
                raise RuntimeError,"FOUND MULTI_MODEL FILE!"
            
            if line.startswith("TER"):
                break
            
            if line.startswith("ATOM"):
                
                atom = pdb_model.PdbAtom( line )
                
                if not chain:
                    chain = atom.chainID
                
                if atom.chainID != chain:
                    raise RuntimeError," FOUND ADDITIONAL CHAIN"
                    break
                    
                # First atom in first residue
                if readingResSeq == None:
                    readingResSeq = atom.resSeq
                    readingResName = atom.resName
                    _atomTypes.append( atom.name.strip() )
                    continue
                
                if readingResSeq != atom.resSeq:
                    # Adding a new residue
                    
                    # Add the atom we've just finished reading
                    resName.append( readingResName )
                    resSeq.append( readingResSeq )
                    atomTypesList.append( _atomTypes )
                    
                    # Reset
                    readingResSeq = atom.resSeq
                    readingResName = atom.resName
                    _atomTypes = [ atom.name.strip() ]
                else:
                    if atom.name not in _atomTypes:
                        _atomTypes.append( atom.name.strip() )
                        
        # End reading loop

        # Add the atom we've just finished reading
        resName.append( readingResName )
        resSeq.append( readingResSeq )
        atomTypesList.append( _atomTypes )
        
        sequence = ""
        # Build up the sequence
        for n in resName:
            sequence += pdb_edit.three2one[ n ]
        
        # Build up the mask
        cAlphaMask = []
        for atomTypes in atomTypesList:
            if 'CA' not in atomTypes:
                cAlphaMask.append( True )
            else:
                cAlphaMask.append( False )
        
        return ( sequence, resSeq, cAlphaMask )
    
    
    def resSeqMatch(self):
        """Return true if the residue numbering between the model and native over the aligned region is the same"""
        
        #print self.modelResSeq[ self.modelOffset : self.modelOffset + self.lenMatch ]
        #print self.nativeResSeq[ self.nativeOffset : self.nativeOffset + self.lenMatch ]
        return self.modelResSeq[ self.modelOffset : self.modelOffset + self.lenMatch ] == self.nativeResSeq[ self.nativeOffset : self.nativeOffset + self.lenMatch ]


class Test(unittest.TestCase):



    def testResSeqMap1(self):
        """See if we can sort out the indexing between the native and model"""

        resSeqMap = residueSequenceMap()
        
        resSeqMap.modelSequence = ['G', 'G', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'F', 'F', 'F', 'F', 'F', 'F']
        resSeqMap.modelResSeq = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        resSeqMap.modelCAlphaMask = [False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False]
        
        resSeqMap.nativeSequence = [ 'H', 'H', 'H', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A' ]
        resSeqMap.nativeResSeq = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        resSeqMap.nativeCAlphaMask = [False, False, False, False, False, False, True, False, False, False, False, False, False]
        
        resSeqMap._calc_map()
        
        self.assertEqual( resSeqMap.modelOffset, 2)
        self.assertEqual( resSeqMap.nativeOffset, 3)
        self.assertEqual( resSeqMap.lenMatch, 10)
        
        self.assertEqual( resSeqMap.native2model( 0 ), -6 )
        self.assertEqual( resSeqMap.native2model( 3 ), -3 )
        self.assertRaises( ValueError, resSeqMap.native2model, 13 )
        
        self.assertEqual( resSeqMap.model2native( 1 ), 7 )
        self.assertEqual( resSeqMap.model2native( 12 ), 18 )
        self.assertEqual( resSeqMap.model2native( 6 ), 12 )
        
        self.assertEqual( resSeqMap.modelIncomparable(), [-5,-4,-1, 0, 7, 8, 9, 10, 11, 12] )
        self.assertEqual( resSeqMap.nativeIncomparable(), [ 0, 1, 2, 5, 6 ] )
        
        # Check ends match up
        m1 = resSeqMap.modelResSeq[ resSeqMap.modelOffset ]
        n1 = resSeqMap.model2native( m1 )
        self.assertEqual( m1, resSeqMap.native2model(n1) )
        re = resSeqMap.nativeResSeq[ resSeqMap.nativeOffset + resSeqMap.lenMatch - 1  ]
        self.assertEqual( resSeqMap.native2model( re ), resSeqMap.modelResSeq[ resSeqMap.modelOffset + resSeqMap.lenMatch - 1  ] )
        
        
        return
    
    def testRefSeqMap2(self):
        """See if we can sort out the indexing between the native and model"""
        
        
        nativePdb = "../tests/testfiles//2XOV.pdb"
        modelPdb = "../tests/testfiles/2XOV_S_00000001.pdb" 
        
        resSeqMap = residueSequenceMap( nativePdb, modelPdb )
        
        # Check ends match up
        m1 = resSeqMap.modelResSeq[ resSeqMap.modelOffset ]
        n1 = resSeqMap.model2native( m1 )
        self.assertEqual( m1, resSeqMap.native2model(n1) )
        re = resSeqMap.nativeResSeq[ resSeqMap.nativeOffset + resSeqMap.lenMatch - 1  ]
        self.assertEqual( resSeqMap.native2model( re ), resSeqMap.modelResSeq[ resSeqMap.modelOffset + resSeqMap.lenMatch - 1  ] )
        
        return
    
    def testResSeqMap3(self):
        """See if we can sort out the indexing between the native and model"""
        
        
        nativePdb = "../tests/testfiles/2UUI.pdb"
        modelPdb = "../tests/testfiles/2UUI_S_00000001.pdb"
        
        PE = pdb_edit.PDBEdit()
        chainA = "2UUI_A.pdb"
        PE.extract_chain( nativePdb, chainA, chainID='A' )
        chainAstd = "2UUI_A_std.pdb"
        PE.standardise(chainA, chainAstd)
        
        resSeqMap = residueSequenceMap( chainA, modelPdb )
        
        nativeMask = [ False ] * 155 + [ True ]
        self.assertEqual( resSeqMap.nativeCAlphaMask, nativeMask)
        
        self.assertEqual( resSeqMap.native2model(10), 16  )
        self.assertEqual( resSeqMap.model2native(155), 149 )
        
        # Check ends match up
        m1 = resSeqMap.modelResSeq[ resSeqMap.modelOffset ]
        n1 = resSeqMap.model2native( m1 )
        self.assertEqual( m1, resSeqMap.native2model(n1) )
        re = resSeqMap.nativeResSeq[ resSeqMap.nativeOffset + resSeqMap.lenMatch - 1  ]
        self.assertEqual( resSeqMap.native2model( re ), resSeqMap.modelResSeq[ resSeqMap.modelOffset + resSeqMap.lenMatch - 1  ] )
        
        os.unlink( chainA )
        os.unlink( chainAstd )
        
        return
    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
