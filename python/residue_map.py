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
    
    The mapping is based on the model, so it starts at the first model residue and continues to the end.
    
    If the native is longer than the model, we ignore any residues that are outside the model - so you
    will get an index error if you query for residues beyond the length of the model.
    
     """
    
    def __init__( self, nativePdb, modelPdb ):
        
        # Matching and ordered lists of the resSeq in the model and native
        # The model list is "complete" so the first residue in the model will be at pos 0 and
        # the last item will be the last residue
        self.modelResSeq = []
        self.nativeResSeq = []
        
        # An array mapped onto the resSeq arrays that has True where there are no cALpha atoms
        self.cAlphaMask = []
        
        self.calc_map(nativePdb, modelPdb)
        
        return
    
    def native2model( self, nativeResSeq ):
        """Return the model resSeq for the given native resSeq"""
        return self.modelResSeq[ self.nativeResSeq.index( nativeResSeq ) ]
    
    def model2native( self, modelResSeq):
        """Return the native resSeq for the given native resSeq"""
        return self.nativeResSeq[ self.modelResSeq.index( modelResSeq ) ]

    def incomparable(self):
        """Return a list of the indices of any residues that are in the model but not in the native or
        residues for which there are no c-Alpha atoms in one or other of the models"""
        return [ resId for i, resId in enumerate( self.modelResSeq ) if ( self.nativeResSeq[ i ] == None or self.cAlphaMask[ i ] ) ]
    
    def __str__(self):
        
        s = "residueSequenceMap: {0}\n".format( self.__repr__() )
        s += "Model: {0}\n".format( self.modelResSeq )
        s+= "Native: {0}\n".format( self.nativeResSeq )
        s+= "cAlphaMask: {0}\n".format( self.cAlphaMask )
        
        return s
    
    def calc_map( self, nativePdb, modelPdb ):
        
        self.native_seq, self.native_idx, self.native_camask = self.read_pdb( nativePdb )
        self.model_seq, self.model_idx, self.model_camask = self.read_pdb( modelPdb )
        
        self._calc_map()
        
        return
    
    def _calc_map( self ):
        """Return a ResSeqMap mapping the index of a residue in the model to the corresponding residue in the native.
        Only works if 1 chain in either file and with standard residues
        """
        
        if len(self.native_seq) < 10 or len(self.model_seq) < 10:
            raise RuntimeError,"Very short sequences - this will not work!"
        
        # The window of AA we used to check for a match    
        PROBE_LEN = 10
        
        # MAXINSET is the max number of AA into the sequence that we will go searching for a match - i.e. if more
        # then MAXINSET AA are non-matching, we won't find the match 
        l = len( self.model_seq ) if len( self.model_seq ) < len( self.native_seq ) else len( self.native_seq )
        MAXINSET=30 if l > 30 else ( l - PROBE_LEN )
        
        got=False
        for model_i in range( MAXINSET + 1 ):
            probe = self.model_seq[ model_i : model_i+PROBE_LEN ]
            #print "PROBE ",probe
            for native_i in range( MAXINSET + 1 ):
                #print "TEST ",self.native_seq[ native_i:native_i+PROBE_LEN ]
                if self.native_seq[ native_i:native_i+PROBE_LEN ] == probe:
                    got=True
                    break
            
            if got:
#                 print "GOT MODEL MATCH AT i,j ",model_i,native_i
                break
            
        if not got:
            raise RuntimeError,"Could not calculate map!"
        
        # Now we know where they start we can sort out the indicies
        # map goes from the model -> native. For any in model that are not in native we set them to None
        self.modelResSeq = self.model_idx
        
        # reset masks
        self.cAlphaMask = [ False ] * len( self.model_seq )
        
        for i in range( len( self.model_seq ) ):
            
            if i < model_i:
                # These are residues that are present in the model but not in the native
                self.nativeResSeq.append( None )
                
                if self.model_camask[ i ]:
                    self.cAlphaMask[ i ] = True
                    
                continue
            
            pos = i - model_i + native_i
            if pos >= len( self.native_idx ):
                self.nativeResSeq.append( None )
            else:
                self.nativeResSeq.append( self.native_idx[ pos ] )
                # Add any masks
                if self.model_camask[ i ]:
                    self.cAlphaMask[ i ] = True
                if self.native_camask[ pos ]:
                    self.cAlphaMask[ i ] = True
        
        if len(self.nativeResSeq) != len(self.modelResSeq):
            raise RuntimeError, "Mismatching maps: {0}".format( self ) 
            
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
        
        self.assertEqual( refSeqMap.nativeResSeq, nativeResSeq, "map doesn't match")
        
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
        
        self.assertEqual( refSeqMap.nativeResSeq, nativeResSeq, "map doesn't match")
        
        os.unlink( tmp1 )
        os.unlink( tmp2 )
        
        return

    def XtestResSeqMap3(self):
        """See if we can sort out the indexing between the native and model"""

        resSeqMap = residueSequenceMap()
        
        resSeqMap.model_seq = ['G', 'G', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'F', 'F', 'F', 'F', 'F', 'F']
        resSeqMap.model_idx = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        resSeqMap.model_camask = [False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False]
        resSeqMap.model_bbmask = [False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False]
        
        resSeqMap.native_seq = [ 'H', 'H', 'H', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A' ]
        resSeqMap.native_idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        resSeqMap.native_camask = [False, False, False, False, False, False, True, False, False, False, False, False, False]
        resSeqMap.native_bbmask = [True, False, False, False, False, False, False, False, False, False, False, False, False]
        
        resSeqMap._calc_map()
        
        native = [None, None, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, None, None, None, None, None, None]
        cAlphaMask = [False, False, False, False, True, True, False, False, False, False, False, False, False, False, False, False, False, False]
        
        self.assertEqual( resSeqMap.modelResSeq, resSeqMap.model_idx)
        self.assertEqual( resSeqMap.nativeResSeq, native)
        self.assertEqual( resSeqMap.cAlphaMask, cAlphaMask)
        
        
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
        model = [ i for i in range(1,157)]
        native = [ i - 5 for i in range(156)]
        cAlphaMask = [ False ] * 155 + [ True ]
        
        self.assertEqual( resSeqMap.modelResSeq, model)
        self.assertEqual( resSeqMap.nativeResSeq, native)
        self.assertEqual( resSeqMap.cAlphaMask, cAlphaMask)
        
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
