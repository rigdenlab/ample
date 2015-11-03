'''
Created on 2 Feb 2015

@author: jmht
'''
import os
import shutil
import sys
import unittest

# our imports
import ample_util
import mtz_util

if not "CCP4" in os.environ.keys(): raise RuntimeError('CCP4 not found')
mrbumpd=os.path.join(os.environ['CCP4'],"share","mrbump","include","parsers")
sys.path.insert(0,mrbumpd)
import parse_shelxe

def shelxe_origin(shelxe_exe, native_pdb, native_mtz, mr_pdb):
    if not ample_util.is_exe(shelxe_exe): raise RuntimeError,"Cannot find shelxe executable: {0}".format(shelxe_exe)
    if not os.path.isfile(native_pdb): raise RuntimeError,"Cannot find native_pdb: {0}".format(native_pdb)
    if not os.path.isfile(native_mtz): raise RuntimeError,"Cannot find native_mtz: {0}".format(native_mtz)
    if not os.path.isfile(mr_pdb): raise RuntimeError,"Cannot find mr_pdb: {0}".format(mr_pdb)
    
    stem="shelxe-input" # stem name for all shelxe files
    hkl_file=stem+".hkl"
    hkl_file=mtz_util.to_hkl(native_mtz, hkl_file=hkl_file)
        
    # Rename nativePdb and mrPdb
    shutil.copyfile(mr_pdb, stem+".pda")
    shutil.copyfile(native_pdb, stem+".ent")
    trace_cycles=0
    frac_solvent=0.5
    cmd=[shelxe_exe,'shelxe-input.pda','-a{0}'.format(trace_cycles),'-q', '-s{0}'.format(frac_solvent),'-o','-n','-t0','-m0','-x']
    logfile='shelxe.log'
    ret = ample_util.run_command(cmd=cmd, logfile=logfile, directory=None, dolog=True, stdin=None)
    if not ret==0:
        raise RuntimeError,"Error running shelxe - see log: {0}".format(logfile)
    else:
        for ext in ['.pda','.hkl','.ent','.pdo','.phs','.lst','_trace.ps']:
            os.unlink(stem+ext)
    
    sp=parse_shelxe.ShelxeLogParser(logfile)
    os.unlink(logfile)
    originShift=[ o*-1 for o in sp.originShift ]
    return originShift

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = cls.thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -1 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        
        cls.shelxe_exe=ample_util.find_exe('shelxe')

        return

    def test1BYZ(self):
        pdb=os.path.join(self.testfiles_dir,"1BYZ.pdb")
        mtz=os.path.join(self.testfiles_dir,"1BYZ-cad.mtz")
        mrPdb=os.path.join(self.testfiles_dir,"1BYZ_phaser_loc0_ALL_poly_ala_trunc_0.486615_rad_1_UNMOD.1.pdb")
        origin=shelxe_origin(self.shelxe_exe,pdb,mtz,mrPdb)
        self.assertEqual(origin,[0.326, 0.19, 0.275])
        return
    
    def test1D7M(self):
        pdb=os.path.join(self.testfiles_dir,"1D7M.pdb")
        mtz=os.path.join(self.testfiles_dir,"1D7M-cad.mtz")
        mrPdb=os.path.join(self.testfiles_dir,"1D7M_phaser_loc0_ALL_SCWRL_reliable_sidechains_trunc_5.241154_rad_1_UNMOD.1.pdb")
        origin=shelxe_origin(self.shelxe_exe,pdb,mtz,mrPdb)
        self.assertEqual(origin,[-0.0, -0.0, 0.5])
        return

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()