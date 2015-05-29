'''
Created on 26 May 2015

@author: jmht
'''
import ample_sequence
import ample_util
import glob
import logging
import os
import pdb_edit
import shutil
import sys
import unittest

class Theseus(object):
    
    def __init__(self, work_dir=None, theseus_exe=None):
        
        self.theseus_exe = theseus_exe
        if not os.path.exists(self.theseus_exe) and os.access(self.theseus_exe, os.X_OK):
            raise RuntimeError,"Cannot find theseus_exe: {0}".format(self.theseus_exe)
        
        self.logger = logging.getLogger()
        
        self.work_dir = None
        self.variance_log = None
        self.superposed_models = None
        self.aligned_models = None
        
        self._set_work_dir(work_dir)
        return
    
    def _set_work_dir(self,work_dir):
        if work_dir: self.work_dir = work_dir
        if not self.work_dir: self.work_dir = os.getcwd()
        if not os.path.isdir(self.work_dir): os.mkdir(self.work_dir)
        return self.work_dir
    
    def alignment_file(self, models, alignment_file=None):
        if not alignment_file: alignment_file = os.path.join(self.work_dir,'homologs.fasta')
        all_seq = ample_sequence.Sequence(pdb=models[0])
        for model in models[1:]: all_seq += ample_sequence.Sequence(pdb=model)
        all_seq.write_fasta(alignment_file,pdbname=True)
        return alignment_file

    def align_models(self, models, work_dir=None, basename=None, homologs=False):
        self._set_work_dir(work_dir)
        if homologs:
            # Theseus expects all the models to be in the directory that it is run in as the string
            # given in the fasta header is used to construct the file names of the aligned pdb files
            # If a full or relative path is given (e.g. /foo/bar.pdb), it tries to create files called "basename_/foo/bar.pdb"
            # We therefore copy the models in and then delete them afterwards
            alignment_file = self.alignment_file(models)
            copy_models = [ os.path.join(self.work_dir,os.path.basename(m)) for m in models ]
            for orig, copy in zip(models, copy_models): shutil.copy(orig, copy)
        
        if not basename: basename = 'theseus'

        cmd = [ self.theseus_exe, '-a0', '-r', basename ]
        if homologs:
            cmd += [ '-A', alignment_file ]
            cmd += [ os.path.basename(m) for m in copy_models ]
        else:
            cmd += models
        
        self.theseus_log=os.path.join(self.work_dir,"theseus.log")
        retcode = ample_util.run_command(cmd,
                                         logfile=self.theseus_log,
                                         directory=self.work_dir)
        if retcode != 0:
            msg = "non-zero return code for theseus in align_models!\n See log: {0}".format(self.theseus_log)
            self.logger.critical(msg)
            raise RuntimeError, msg
        
        self.variance_file = os.path.join(self.work_dir,'{0}_variances.txt'.format(basename))
        self.superposed_models = os.path.join(self.work_dir,'{0}_sup.pdb'.format(basename))
        if homologs:
            self.aligned_models = [ os.path.join(self.work_dir,"theseus_{0}".format(os.path.basename(m))) for m in copy_models ]
            for m in copy_models: os.unlink(m)
        
        return self.superposed_models

    def var_by_res(self):
        """Return a list of tuples: (resSeq,variance)"""
        
        #--------------------------------
        # get variations between pdbs
        #--------------------------------
        if not os.path.isfile(self.variance_file):
            raise RuntimeError,"Cannot find theseus variance file: {0} Please check the log: {1}".format(self.variance_file,
                                                                                                         self.theseus_log)
        variances=[]
        with open(self.variance_file) as f:
            for i, line in enumerate(f):
                # Skip header
                if i==0: continue

                line=line.strip()
                if not line: continue # Skip blank lines

                #print line
                tokens=line.split()
                # Different versions of theseus may have a RES card first, so need to check
                if tokens[0]=="RES":
                    idxidx=1
                    idxResSeq=3
                    idxVariance=4
                else:
                    idxidx=0
                    idxResSeq=2
                    idxVariance=3
                idx = int(tokens[idxidx])
                assert idx == i,"Index and atom lines don't match! {0} : {1}".format(idx,i) # paranoid check
                # Theseus counts from 1, we count from 0
                idx -= 1
                resSeq = int(tokens[idxResSeq])
                variance = float(tokens[idxVariance])
                variances.append((idx,resSeq,variance))
                
        return variances


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
        
        cls.theseus_exe=ample_util.find_exe('theseus')

        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        formatter = logging.Formatter('%(message)s')
        ch.setFormatter(formatter)
        root.addHandler(ch)

        return

    def testAlignModels(self):
        """Test we can reproduce the original thresholds"""
        os.chdir(self.thisd)
        
        models = glob.glob(os.path.join(self.testfiles_dir,'models','*.pdb'))
        work_dir = os.path.join(self.tests_dir,'theseus_align')
        homologs = False
        rtheseus = Theseus(work_dir=work_dir,theseus_exe=self.theseus_exe)
        rtheseus.align_models(models,homologs=homologs)
        var_by_res = rtheseus.var_by_res()
        ref = [(0, 1, 58.093855), (1, 2, 49.037612), (2, 3, 49.9941), (3, 4, 41.759792), (4, 5, 37.227847), 
               (5, 6, 27.3795), (6, 7, 25.348492), (7, 8, 25.799824), (8, 9, 22.432552), (9, 10, 23.265923), 
               (10, 11, 23.050341), (11, 12, 20.235297), (12, 13, 18.29234), (13, 14, 16.800248), (14, 15, 16.07131), 
               (15, 16, 10.678152), (16, 17, 10.77214), (17, 18, 6.206533), (18, 19, 5.665422), (19, 20, 3.152776), 
               (20, 21, 1.860673), (21, 22, 0.705796), (22, 23, 0.185265), (23, 24, 0.116864), (24, 25, 0.10304), 
               (25, 26, 0.041897), (26, 27, 0.039639), (27, 28, 0.084557), (28, 29, 0.090918), (29, 30, 0.044023), 
               (30, 31, 0.022351), (31, 32, 0.018216), (32, 33, 0.081551), (33, 34, 0.085867), (34, 35, 0.65924), 
               (35, 36, 0.874891), (36, 37, 1.772875), (37, 38, 3.509799), (38, 39, 4.900248), (39, 40, 7.269827), 
               (40, 41, 10.250037), (41, 42, 15.287341), (42, 43, 23.696658), (43, 44, 30.767263), (44, 45, 35.877388), 
               (45, 46, 36.17383), (46, 47, 41.678763), (47, 48, 53.733142), (48, 49, 56.091737), (49, 50, 50.992341), 
               (50, 51, 70.052003), (51, 52, 60.586176), (52, 53, 43.26875), (53, 54, 58.925292), (54, 55, 74.272249), 
               (55, 56, 59.917994), (56, 57, 56.26098), (57, 58, 80.235351), (58, 59, 86.028986)]
        self.assertEqual(var_by_res,ref)
        shutil.rmtree(work_dir)
        return
    
    def testAlignModelsHomo(self):
        """Test we can reproduce the original thresholds"""
        os.chdir(self.thisd)

        work_dir = os.path.join(self.tests_dir,'theseus_align_homo')
        os.mkdir(work_dir)
        pdb_list = [ '1D7M.pdb', '1GU8.pdb', '2UUI.pdb', '1K33.pdb' ,'1BYZ.pdb' ]
        models = []
        tokeep_idx = [ i for i in range(12) ]
        for pdb in pdb_list:
            pdbin = os.path.join(self.testfiles_dir,pdb)
            name = os.path.splitext(pdb)[0]
            pdbout = os.path.join(self.testfiles_dir,"{0}_cut.pdb".format(name))
            pdb_edit.select_residues(pdbin, pdbout, tokeep_idx=tokeep_idx)
            models.append(pdbout)

        homologs = True
        rtheseus = Theseus(work_dir=work_dir,theseus_exe=self.theseus_exe)
        rtheseus.align_models(models,homologs=homologs)
        var_by_res = rtheseus.var_by_res()
        ref = [(0, 243, 9.918397), (1, 244, 3.897504), (2, 245, 1.877927), (3, 246, 2.004033), (4, 247, 1.24683), 
               (5, 248, 0.753177), (6, 249, 0.005146), (7, 250, 0.02917), (8, 251, 0.04054), (9, 252, 0.027774), 
               (10, 253, 0.093861), (11, 254, 0.0)]
        self.assertEqual(var_by_res,ref)
        
        self.assertTrue(all([os.path.isfile(m) for m in rtheseus.aligned_models]))
        # clean up
        for m in models: os.unlink(m)
        shutil.rmtree(work_dir)
        return

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()