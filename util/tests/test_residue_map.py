"""Test functions for util.residue_map"""

import os
import unittest
from ample.util import pdb_edit
from ample.util import residue_map
from ample.testing import constants

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        cls.ample_dir = constants.AMPLE_DIR
        cls.tests_dir=os.path.join(cls.ample_dir,"testing")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')

    def test_resSeqMap1(self):
        # See if we can sort out the indexing between the native and model

        resSeqMap = residue_map.residueSequenceMap()
        
        resSeqMap.targetSequence = ['G', 'G', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'F', 'F', 'F', 'F', 'F', 'F']
        resSeqMap.targetResSeq = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        resSeqMap.targetCAlphaMask = [False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False]
        
        resSeqMap.refSequence = [ 'H', 'H', 'H', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'G', 'G' ]
        resSeqMap.refResSeq = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 ]
        resSeqMap.refCAlphaMask = [False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, False]
        
        resSeqMap._calc_map()
        
        self.assertEqual( resSeqMap.targetOffset, 2)
        self.assertEqual( resSeqMap.refOffset, 3)
        self.assertEqual( resSeqMap._lenMatch(), 10)
        
        self.assertEqual( resSeqMap.ref2target( 0 ), -6 )
        self.assertEqual( resSeqMap.ref2target( 3 ), -3 )
        
        self.assertEqual( resSeqMap.target2ref( 1 ), 7 )
        self.assertEqual( resSeqMap.target2ref( 12 ), 18 )
        self.assertEqual( resSeqMap.target2ref( 6 ), 12 )
        
        self.assertEqual( resSeqMap.targetIncomparable(), [-5,-4,-1, 0, 7, 8, 9, 10, 11, 12] )
        self.assertEqual( resSeqMap.refIncomparable(), [ 0, 1, 2, 5, 6, 12, 13, 14, 15, 16, 17, 18 ] )
        
        # Check ends match up
        m1 = resSeqMap.targetResSeq[ resSeqMap.targetOffset ]
        n1 = resSeqMap.target2ref( m1 )
        self.assertEqual( m1, resSeqMap.ref2target(n1) )
        re = resSeqMap.refResSeq[ resSeqMap.refOffset + resSeqMap.lenMatch - 1  ]
        self.assertEqual( resSeqMap.ref2target( re ), resSeqMap.targetResSeq[ resSeqMap.targetOffset + resSeqMap.lenMatch - 1  ] )
    
    def test_resSeqMap2(self):
        # See if we can sort out the indexing between the native and model
        
        
        nativePdb = os.path.join(self.testfiles_dir,"2XOV.pdb")
        modelPdb = os.path.join(self.testfiles_dir,"2XOV_S_00000001.pdb")
        
        resSeqMap = residue_map.residueSequenceMap( nativePdb, modelPdb )
        
        self.assertEqual( 181, resSeqMap._lenMatch() )
        # Check ends match up
        m1 = resSeqMap.targetResSeq[ resSeqMap.targetOffset ]
        n1 = resSeqMap.target2ref( m1 )
        self.assertEqual( m1, resSeqMap.ref2target(n1) )
        re = resSeqMap.refResSeq[ resSeqMap.refOffset + resSeqMap.lenMatch - 1  ]
        self.assertEqual( resSeqMap.ref2target( re ), resSeqMap.targetResSeq[ resSeqMap.targetOffset + resSeqMap.lenMatch - 1  ] )
    
    def test_resSeqMap3(self):
        # See if we can sort out the indexing between the native and model
        
        nativePdb = os.path.join(self.testfiles_dir,"2UUI.pdb")
        modelPdb = os.path.join(self.testfiles_dir,"2UUI_S_00000001.pdb")
        
        chainA = "2UUI_A.pdb"
        pdb_edit.extract_chain( nativePdb, chainA, chainID='A' )
        chainAstd = "2UUI_A_std.pdb"
        pdb_edit.standardise(chainA, chainAstd)
        
        resSeqMap = residue_map.residueSequenceMap( chainA, modelPdb )
        
        self.assertEqual( 156, resSeqMap._lenMatch() )

        
        nativeMask = [ False ] * 155 + [ True ]
        self.assertEqual( resSeqMap.refCAlphaMask, nativeMask)
        
        self.assertEqual( resSeqMap.ref2target(10), 16  )
        self.assertEqual( resSeqMap.target2ref(155), 149 )
        
        # Check ends match up
        m1 = resSeqMap.targetResSeq[ resSeqMap.targetOffset ]
        n1 = resSeqMap.target2ref( m1 )
        self.assertEqual( m1, resSeqMap.ref2target(n1) )
        re = resSeqMap.refResSeq[ resSeqMap.refOffset + resSeqMap.lenMatch - 1  ]
        self.assertEqual( resSeqMap.ref2target( re ), resSeqMap.targetResSeq[ resSeqMap.targetOffset + resSeqMap.lenMatch - 1  ] )
        
        os.unlink( chainA )
        os.unlink( chainAstd )
    
    def test_resSeqMap4(self):
        # See if we can sort out the indexing between the native and model
        
        
        nativePdb = os.path.join(self.testfiles_dir,"1K33.pdb")
        modelPdb = os.path.join(self.testfiles_dir,"1K33_S_00000001.pdb")
        
        nativePdbStd = "1K33_std.pdb"
        pdb_edit.standardise( nativePdb, nativePdbStd )
        
        nativeInfo = pdb_edit.get_info( nativePdbStd )
        modelInfo = pdb_edit.get_info( modelPdb )
        
        resSeqMap = residue_map.residueSequenceMap( )
        resSeqMap.fromInfo( nativeInfo, 'A', modelInfo, 'A' )
        
        os.unlink( nativePdbStd )

if __name__ == "__main__":
    unittest.main()
