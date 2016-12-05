"""Test functions for cphasematch"""

import os
import unittest
from ample import constants
from ample.util import cphasematch
import iotbx.mtz

class Test( unittest.TestCase ):

    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share,'testfiles')
        
    def test_merge_mtz(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()        
        mtz1 =  os.path.join(self.ample_share, "examples", 'toxd-example', 'input', '1dtx.mtz')
        mtz1_labels = ['FP', 'SIGFP']
        mtz2 =  os.path.join(self.testfiles_dir, "phaser_loc0_ALL_c1_tl49_r1_allatom_UNMOD.1.mtz")
        mtz2_labels = ['FC', 'PHIC']
        merged_mtz, mtz_labels = cphasematch.merge_mtz(mtz1, mtz1_labels, mtz2, mtz2_labels)
        
        input_labels = set(mtz1_labels + mtz2_labels)
        
        # Get all the labels
        miller_dict = iotbx.mtz.object(file_name=merged_mtz).as_miller_arrays_dict()
        ref_labels = set([d[2] for d in miller_dict.keys()])
        
        self.assertEqual(len(input_labels - ref_labels), 0, "Labels were not identical")
        
        return
    
    def test_merge_mtz2(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()        
        mtz1 =  os.path.join(self.ample_share, "examples", 'toxd-example', 'input', '1dtx.mtz')
        mtz1_labels = ['FP', 'SIGFP']
        mtz2 =  os.path.join(self.ample_share, "examples", 'toxd-example', 'input', '1dtx.mtz')
        mtz2_labels = ['FP', 'SIGFP']

        merged_mtz, mtz_labels = cphasematch.merge_mtz(mtz1, mtz1_labels, mtz2, mtz2_labels)
        
        # Get all the labels
        miller_dict = iotbx.mtz.object(file_name=merged_mtz).as_miller_arrays_dict()
        merged_labels = set([d[2] for d in miller_dict.keys()])
        
        # An integer equal to the file number should have been added
        ref_labels = set(['FP','SIGFP','FP2', 'SIGFP2'])
        
        self.assertEqual(len(merged_labels - ref_labels), 0, "Labels were not identical")
        
        return
    
    def test_cphasematch_pdb(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()        
        native_pdb =  os.path.join(self.ample_share, "examples", 'toxd-example', 'input', '1DTX.pdb')
        native_mtz =  os.path.join(self.ample_share, "examples", 'toxd-example', 'input', '1dtx.mtz')
        mr_mtz =  os.path.join(self.testfiles_dir, "phaser_loc0_ALL_c1_tl49_r1_allatom_UNMOD.1.mtz")
        before_origin, after_origin, change_of_hand, origin_shift = cphasematch.calc_phase_error_pdb(native_pdb,
                                                                                                     native_mtz,
                                                                                                     mr_mtz,
                                                                                                     f_label='FP',
                                                                                                     sigf_label='SIGFP')
        self.assertEqual(89.6763, before_origin)
        self.assertEqual(62.4365, after_origin)
        self.assertEqual([0.0, 0.0, 0.5], origin_shift)
        return
    
    def test_cphasematch_mtz(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()        
        native_mtz_phased =  os.path.join(self.testfiles_dir, "toxd_59.1.mtz")
        mr_mtz =  os.path.join(self.testfiles_dir, "phaser_loc0_ALL_c1_tl49_r1_allatom_UNMOD.1.mtz")
        before_origin, after_origin, change_of_hand, origin_shift = cphasematch.calc_phase_error_mtz(native_mtz_phased,
                                                                                                     mr_mtz,
                                                                                                     f_label='FP',
                                                                                                     sigf_label='SIGFP')
        self.assertEqual(88.8476, before_origin)
        self.assertEqual(62.4833, after_origin)
        self.assertEqual([0.0, 0.5, 0.0], origin_shift)
        return
    
    @unittest.skip("Work in progress")
    def test_cphasematch_mtz_origin(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()        
        native_mtz_phased =  os.path.join(self.testfiles_dir, "toxd_59.1.mtz")
        mr_mtz =  os.path.join(self.testfiles_dir, "phaser_loc0_ALL_c1_tl49_r1_allatom_UNMOD.1.mtz")
        origin = [0.0, 0.5, 0.0]
        before_origin, after_origin, change_of_hand, origin_shift = cphasematch.calc_phase_error_mtz(native_mtz_phased,
                                                                                                     mr_mtz,
                                                                                                     origin=origin)
        # Can't test exact equality as cphasematch and cctbx return slightly different errors
        self.assertAlmostEqual(88.8476, before_origin, 0)
        self.assertAlmostEqual(62.4833, after_origin, 0)
        self.assertEqual([0.0, 0.5, 0.0], origin_shift)
        return
    
    @unittest.skip("Work in progress")
    def test_cphasematch_shelxe(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()        
        native_mtz_phased =  os.path.join(self.testfiles_dir, "toxd_59.1.mtz")
        shelxe_mtz =  os.path.join(self.testfiles_dir, "shelxe_phaser_loc0_ALL_c1_t49_r1_polyAla_UNMOD.mtz")
        origin = [0.0, 0.0, 0.0]
        before_origin, after_origin, change_of_hand, origin_shift = cphasematch.calc_phase_error_mtz(native_mtz_phased,
                                                                                                     shelxe_mtz,
                                                                                                     origin=origin)
        # Can't test exact equality as cphasematch and cctbx return slightly different errors
        self.assertAlmostEqual(88.8476, before_origin, 0)
        self.assertAlmostEqual(62.4833, after_origin, 0)
        self.assertEqual([0.0, 0.5, 0.0], origin_shift)
        return
    

if __name__ == "__main__":
    unittest.main()
