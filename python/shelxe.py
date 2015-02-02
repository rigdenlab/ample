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

# mrbump import
if not "CCP4" in os.environ.keys():
    raise RuntimeError('CCP4 not found')
mrbumpd=os.path.join(os.environ['CCP4'],"include","mrbump","include","parsers")
mrbumpd="/opt/mrbump-trunk/include/parsers"
sys.path.insert(0,mrbumpd)
import parse_shelxe

def shelxeOrigin(shelxeExe,nativePdb,nativeMtz,mrPdb):
    
    stem="shelxe-input" # stem name for all shelxe files
    
    print "Parsing MTZ file {0} to determine column labels".format(nativeMtz)
    F,SIGF,FREE=mtz_util.getLabels(nativeMtz)
    
    nativeHkl=stem+".hkl"
    print "Creating HKL format file".format(nativeHkl)
    
    cmd=['mtz2various','HKLIN',nativeMtz,'HKLOUT', nativeHkl]
    logfile="mtz2various.log"
    stdin  = """LABIN FP={0} SIGFP={1} FREE={2}
OUTPUT SHELX
FSQUARED
END""".format(F,SIGF,FREE)
    
    ret = ample_util.run_command(cmd=cmd, logfile=logfile, directory=None, dolog=False, stdin=stdin)
    if not ret==0:
        raise RuntimeError,"Error converting {0} to HKL format - see log: {1}".format(nativeMtz,logfile)
    else:
        os.unlink(logfile)
        
    # Rename nativePdb and mrPdb
    shutil.copyfile(mrPdb, stem+".pda")
    shutil.copyfile(nativePdb, stem+".ent")
    traceCycles=0
    fracSolvent=0.5
    cmd=[shelxeExe,'shelxe-input.pda','-a{0}'.format(traceCycles),'-q', '-s{0}'.format(fracSolvent),'-o','-n','-t0','-m0','-x']
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
        origin=shelxeOrigin(self.shelxe_exe,pdb,mtz,mrPdb)
        self.assertEqual(origin,[0.326, 0.19, 0.275])
        return
    
    def test1D7M(self):
        pdb=os.path.join(self.testfiles_dir,"1D7M.pdb")
        mtz=os.path.join(self.testfiles_dir,"1D7M-cad.mtz")
        mrPdb=os.path.join(self.testfiles_dir,"1D7M_phaser_loc0_ALL_SCWRL_reliable_sidechains_trunc_5.241154_rad_1_UNMOD.1.pdb")
        origin=shelxeOrigin(self.shelxe_exe,pdb,mtz,mrPdb)
        self.assertEqual(origin,[-0.0, -0.0, 0.5])
        return


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()