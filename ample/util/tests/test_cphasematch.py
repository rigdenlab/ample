"""Test functions for cphasematch"""

import os
import unittest
from ample import constants
from ample.util import cphasematch

class TestContacts( unittest.TestCase ):

    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share,'testfiles')
        
    def test_cphasematch(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()        
        cp = cphasematch.Cphasematch()
        native_mtz =  os.path.join(self.testfiles_dir, "toxd_59.1.mtz")
        placed_mtz =  os.path.join(self.testfiles_dir, "phaser_loc0_ALL_c1_tl49_r1_allatom_UNMOD.1.mtz")
        
        
        cp.run(native_mtz, placed_mtz)
        self.assertEqual('88.8476', cp.before_origin)
        self.assertEqual('62.4833', cp.after_origin)
        
        os.unlink('all_phases.mtz')
        os.unlink('cphasematch.mtz')
        os.unlink('c1_tl49_r1_allato_cphasematch.log')
 
        return

if __name__ == "__main__":
    unittest.main()
