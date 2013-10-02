'''
Useful manipulations on PDB files
'''

# Python imports
import os
import sys
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
        self.modelOffset = None # Where the matched part of the sequences starts in the model
        self._modelIncomparable = None # List of atoms in the model that cannot be compared to the native
        
        self.nativeResSeq = []
        self.nativeSequence = None
        self.nativeCAlphaMask = []
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
        
        s = "residueSequenceMap: {0}\n".format( self.__repr__() )
        s += "Model: {0}\n".format( self._modelResSeqMap )
        s+= "Native: {0}\n".format( self._nativeResSeqMap )
        s+= "cAlphaMask: {0}\n".format( self.cAlphaMask )
        
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
        for i in range( self.modelOffset, len( self.modelSequence ) ):
            modelAA = self.modelSequence[ i ]
            
            y = self.nativeOffset + i
            if y >= len( self.nativeSequence ):
                break
            
            nativeAA = self.nativeSequence[ y ]
            #print "matching {0} to {1}".format( modelAA, nativeAA )
            if modelAA != nativeAA:
                break 
        
        self.lenMatch = i
        #print "GOT LEN MATCH ",self.lenMatch
            
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

    def XtestRefSeqMap(self):
        """See if we can sort out the indexing between the native and model"""
        
        
        nativePdb = "/media/data/shared/TM/3RLB/3RLB.pdb"
        modelPdb = "/media/data/shared/TM/3RLB/models/S_00000001.pdb" 

        #modelSeq = "MHHHHHHHHAMSNSKFNVRLLTEIAFMAALAFIISLIPNTVYGWIIVEIACIPILLLSLRRGLTAGLVGGLIWGILSMITGHAYILSLSQAFLEYLVAPVSLGIAGLFRQKTAPLKLAPVLLGTFVAVLLKYFFHFIAGIIFWSQYAWKGWGAVAYSLAVNGISGILTAIAAFVILIIFVKKFPKLFIHSNY"
        #modelIdx = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192 ]
        #nativeSeq = "NVRLLTEIAFMAALAFIISLIPNTVYGWIIVEIACIPILLLSLRRGLTAGLVGGLIWGILSMITGHAYILSLSQAFLEYLVAPVSLGIAGLFRQKTAPLKLAPVLLGTFVAVLLKYFFHFIAGIIFWSQYAWKGWGAVAYSLAVNGISGILTAIAAFVILIIFVKKFPKLFIHSNY"
        nativeIdx = [ 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182 ]
        
        nativeResSeq = [ None ] * 16 + nativeIdx
        
        PE = pdb_edit.PDBEdit()
        
        tmpf = ample_util.tmpFileName()+".pdb"
        PE.extract_chain( nativePdb, tmpf, chainID='A' )
        
        refSeqMap = PE.get_resseq_map( tmpf, modelPdb )
        
        self.assertEqual( refSeqMap._nativeResSeqMap, nativeResSeq, "map doesn't match")
        
        os.unlink( tmpf )
        
        return
    
    def XtestRefSeqMap2(self):
        """See if we can sort out the indexing between the native and model"""
        
        
        nativePdb = "/media/data/shared/TM/3U2F/3U2F.pdb"
        modelPdb = "/media/data/shared/TM/3U2F/models/S_00000001.pdb" 

        nativeResSeq = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, None, None ]
        
        PE = pdb_edit.PDBEdit()
        
        tmp1 = ample_util.tmpFileName()+".pdb"
        PE.standardise( nativePdb, tmp1 )
        
        tmp2 = ample_util.tmpFileName()+".pdb"
        PE.extract_chain( tmp1, tmp2, chainID='K' )
        
        refSeqMap = PE.get_resseq_map( tmp2, modelPdb )
        
        self.assertEqual( refSeqMap._nativeResSeqMap, nativeResSeq, "map doesn't match")
        
        os.unlink( tmp1 )
        os.unlink( tmp2 )
        
        return

    def testResSeqMap3(self):
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
        
        
        return
    
    def testResSeqMap4(self):
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
        
        os.unlink( chainA )
        os.unlink( chainAstd )
        
        return

    def XtestKeepMatching(self):
        """XX"""
        
        nativePdb = "/media/data/shared/TM/3U2F/3U2F.pdb"
        modelPdb = "/media/data/shared/TM/3U2F/models/S_00000001.pdb"
        refinedPdb = "/media/data/shared/TM/3U2F/ROSETTA_MR_0/MRBUMP/cluster_1/search_poly_ala_trunc_0.21093_rad_2_molrep_mrbump/data/loc0_ALL_poly_ala_trunc_0.21093_rad_2/unmod/mr/molrep/refine/refmac_molrep_loc0_ALL_poly_ala_trunc_0.21093_rad_2_UNMOD.pdb"
# 
        PE = pdb_edit.PDBEdit()
        tmp1 = ample_util.tmpFileName()+".pdb"
        PE.standardise( nativePdb, tmp1 )
         
        nativec1 = ample_util.tmpFileName()+".pdb"
        PE.extract_chain( tmp1, nativec1, chainID='K' )
        resSeqMap = PE.get_resseq_map( nativec1, modelPdb )
#         
        refinedc1 = "refinedc1.pdb"
        PE.extract_chain( refinedPdb, refinedc1, chainID='A' )
        print "nativec1 ",nativec1
        matching = "matching.pdb"
        PE.keep_matching( refpdb=refinedc1, targetpdb=nativec1, outpdb=matching, resSeqMap=resSeqMap )
       
#         refpdb = "ref.pdb"
#         targetpdb = "target.pdb"
#         PE._keep_matching( refpdb=refpdb, targetpdb=targetpdb, outpdb=matching, resSeqMap=resSeqMap )
         
        for f in [ tmp1, nativec1, refinedc1, matching ]:
            os.unlink( f )
         
        return
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
