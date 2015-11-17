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
        """Create an alignment file for the models - this is based on the assumption they are all the same length
        but may have different residues"""
        if not alignment_file: alignment_file = os.path.join(self.work_dir,'homologs.fasta')
        all_seq = ample_sequence.Sequence(pdb=models[0])
        for model in models[1:]: all_seq += ample_sequence.Sequence(pdb=model)
        if not all(map(lambda x: x == len(all_seq.sequences[0]), [ len(s) for s in all_seq.sequences ])):
            raise RuntimeError,'PDB files are not all of the same length!\n{0}'.format(models)
        all_seq.write_fasta(alignment_file,pdbname=True)
        return alignment_file

    def align_models(self, models, work_dir=None, basename=None, homologs=False, alignment_file=None):
        self._set_work_dir(work_dir)
        if not basename: basename = 'theseus'
        if homologs:
            # Theseus expects all the models to be in the directory that it is run in as the string given in 
            # the fasta header is used to construct the file names of the aligned pdb files. If a full or 
            # relative path is given (e.g. /foo/bar.pdb), it tries to create files called "basename_/foo/bar.pdb"
            # We therefore copy the models in and then delete them afterwards
            if not alignment_file: alignment_file = self.alignment_file(models)
            copy_models = [ os.path.join(self.work_dir,os.path.basename(m)) for m in models ]
            for orig, copy in zip(models, copy_models): shutil.copy(orig, copy)
        
        # -Z included so we don't line the models up to the principle axis and can compare the ensembles
        #cmd = [ self.theseus_exe, '-a0', '-r', basename, '-Z', '-o', os.path.basename(copy_models[0]) ]
        cmd = [ self.theseus_exe, '-a0', '-r', basename ]
        if homologs:
            cmd += [ '-A', alignment_file ]
            cmd += [ os.path.basename(m) for m in copy_models ]
        else:
            cmd += [ os.path.relpath(m,self.work_dir) for m in models ]
        
        self.theseus_log=os.path.join(self.work_dir,"theseus.log")
        retcode = ample_util.run_command(cmd,
                                         logfile = self.theseus_log,
                                         directory = self.work_dir)
        if retcode != 0:
            msg = "non-zero return code for theseus in align_models!\n See log: {0}".format(self.theseus_log)
            self.logger.critical(msg)
            raise RuntimeError, msg
        
        self.variance_file = os.path.join(self.work_dir,'{0}_variances.txt'.format(basename))
        self.superposed_models = os.path.join(self.work_dir,'{0}_sup.pdb'.format(basename))
        if homologs:
            # Horrible - need to rename the models so that they match the names in the alignment file
            self.aligned_models = []
            for m in copy_models:
                mb = os.path.basename(m)
                aligned_model = os.path.join(self.work_dir,"{0}_{1}".format(basename,mb))
                os.unlink(m)
                os.rename(aligned_model, os.path.join(self.work_dir,mb))
                self.aligned_models.append(mb)
        
        return self.superposed_models

    def var_by_res(self, homologs=False):
        """Return a list of tuples: (resSeq,variance)"""
        
        #--------------------------------
        # get variations between pdbs
        #--------------------------------
        if not os.path.isfile(self.variance_file):
            raise RuntimeError,"Cannot find theseus variance file: {0} Please check the log: {1}".format(self.variance_file,
                                                                                                         self.theseus_log)
        variances=[]
        core_count = 0
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
                    idxCore = 7
                else:
                    idxidx=0
                    idxResSeq=2
                    idxVariance=3
                    idxCore = 6
                    
                if homologs and (len(tokens) < idxCore + 1 or tokens[idxCore] != 'CORE'): continue
                if homologs:
                    idx = core_count
                    core_count += 1
                else:
                    idx = int(tokens[idxidx]) - 1 # Theseus counts from 1, we count from 0
                
                #assert idx == i,"Index and atom lines don't match! {0} : {1}".format(idx,i) # paranoid check
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
        # Below with theseus 3.1.1 on osx 10.9.5
        ref = [(0, 1, 55.757593), (1, 2, 46.981238), (2, 3, 47.734236), (3, 4, 39.857326), (4, 5, 35.477433),
               (5, 6, 26.066719), (6, 7, 24.114493), (7, 8, 24.610988), (8, 9, 21.187142), (9, 10, 21.882375),
               (10, 11, 21.622263), (11, 12, 18.680601), (12, 13, 16.568074), (13, 14, 14.889583), (14, 15, 13.889769),
               (15, 16, 8.722903), (16, 17, 8.719501), (17, 18, 4.648107), (18, 19, 4.263961), (19, 20, 2.338545),
               (20, 21, 1.412784), (21, 22, 0.57754), (22, 23, 0.204917), (23, 24, 0.226518), (24, 25, 0.162323),
               (25, 26, 0.068066), (26, 27, 0.057023), (27, 28, 0.135811), (28, 29, 0.145613), (29, 30, 0.081845),
               (30, 31, 0.051059), (31, 32, 0.045182), (32, 33, 0.112322), (33, 34, 0.102072), (34, 35, 0.446003),
               (35, 36, 0.504418), (36, 37, 1.276947), (37, 38, 2.641781), (38, 39, 4.336794), (39, 40, 6.484846),
               (40, 41, 9.559536), (41, 42, 14.467942), (42, 43, 22.818975), (43, 44, 29.55385), (44, 45, 34.692256),
               (45, 46, 35.141769), (46, 47, 40.41399), (47, 48, 52.268871), (48, 49, 54.535848), (49, 50, 49.527155),
               (50, 51, 67.9861), (51, 52, 58.661069), (52, 53, 41.802971), (53, 54, 57.085415), (54, 55, 71.944127),
               (55, 56, 57.893953), (56, 57, 54.34137), (57, 58, 77.736775), (58, 59, 83.279371)]
        
        self.assertEqual([x[0] for x in var_by_res],[x[0] for x in ref])
        self.assertEqual([x[1] for x in var_by_res],[x[1] for x in ref])
        for i,(t,r) in enumerate(zip([x[2] for x in var_by_res], [x[2] for x in ref])):
            self.assertTrue(abs(t-r) < 0.0001,"Mismatch for: {0} {1} {2}".format(i,t,r))
            
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
        
        #models = [ os.path.join(self.testfiles_dir,pdb) for pdb in pdb_list]

        homologs = True
        rtheseus = Theseus(work_dir=work_dir, theseus_exe=self.theseus_exe)
        rtheseus.align_models(models, homologs=homologs)
        var_by_res = rtheseus.var_by_res()
        # Below with theseus 3.1.1 on osx 10.9.5
        ref  = [(0, 243, 8.049061), (1, 244, 2.614031), (2, 245, 1.343609), (3, 246, 2.261761), (4, 247, 1.112115),
                (5, 248, 0.574936), (6, 249, 0.03114), (7, 250, 0.002894), (8, 251, 0.002314), (9, 252, 0.002174),
                (10, 253, 0.016252), (11, 254, 0.109965)]

        self.assertEqual([x[0] for x in var_by_res],[x[0] for x in ref])
        self.assertEqual([x[1] for x in var_by_res],[x[1] for x in ref])
        for i,(t,r) in enumerate(zip([x[2] for x in var_by_res], [x[2] for x in ref])):
            self.assertTrue(abs(t-r) < 0.0001,"Mismatch for: {0} {1} {2}".format(i,t,r))

        self.assertTrue(all([os.path.isfile(os.path.join(work_dir,m)) for m in rtheseus.aligned_models]))
        # clean up
        for m in models: os.unlink(m)
        shutil.rmtree(work_dir)
        return

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()