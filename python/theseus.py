'''
Created on 26 May 2015

@author: jmht
'''
import collections
import glob
import logging
import os
import shutil
import sys

# local imports
from ample.util import ample_util
from ample.util import sequence_util
from ample.util import pdb_edit

# We create this here otherwise it causes problems with pickling
TheseusVariances = collections.namedtuple('TheseusVariances', ['idx', 'resName', 'resSeq', 'variance', 'stdDev', 'rmsd', 'core'])

class Theseus(object):
    
    def __init__(self, work_dir=None, theseus_exe=None):
        
        self.theseus_exe = theseus_exe
        if theseus_exe is None or not os.path.exists(self.theseus_exe) and os.access(self.theseus_exe, os.X_OK):
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
        all_seq = sequence_util.Sequence(pdb=models[0])
        for model in models[1:]: all_seq += sequence_util.Sequence(pdb=model)
        if not all(map(lambda x: x == len(all_seq.sequences[0]), [ len(s) for s in all_seq.sequences ])):
            raise RuntimeError,'PDB files are not all of the same length!\n{0}'.format(models)
        all_seq.write_fasta(alignment_file,pdbname=True)
        return alignment_file

    def superpose_models(self, models, work_dir=None, basename=None, homologs=False, alignment_file=None):
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
        
        self.theseus_log = os.path.join(self.work_dir,"tlog_{0}.log".format(basename))
        retcode = ample_util.run_command(cmd,
                                         logfile = self.theseus_log,
                                         directory = self.work_dir)
        if retcode != 0:
            msg = "non-zero return code for theseus in superpose_models!\n See log: {0}".format(self.theseus_log)
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
    
    def parse_variances(self, variance_file):
        if not os.path.isfile(variance_file): raise RuntimeError,"Cannot find theseus variance file: {0}".format(variance_file)
        data = []
        
        with open(variance_file) as f:
            for i, line in enumerate(f):
                # Skip header
                if i==0: continue

                line=line.strip()
                if not line: continue # Skip blank lines

                #print line
                tokens = line.split()
                # Different versions of theseus may have a RES card first, so need to check
                if tokens[0]=="RES":
                    idxidx = 1
                    idxResName = 2
                    idxResSeq = 3
                    idxVariance = 4
                    idxStdDev = 5
                    idxRmsd = 6
                    idxCore = 7
                else:
                    idxidx = 0
                    idxResName = 1
                    idxResSeq = 2
                    idxVariance = 3
                    idxStdDev = 4
                    idxRmsd = 5
                    idxCore = 6
                
                # Core may or may not be there
                core = False
                if len(tokens) > idxCore and tokens[idxCore] == 'CORE': core = True
                data.append( TheseusVariances( idx = int(tokens[idxidx]) - 1, # Theseus counts from 1, we count from 0,
                                               resName = tokens[idxResName],
                                               resSeq = int(tokens[idxResSeq]),
                                               variance = float(tokens[idxVariance]),
                                               stdDev = float(tokens[idxStdDev]),
                                               rmsd = float(tokens[idxRmsd]),
                                               core = core ) )
        return data

    def var_by_res(self, core=False):
        """Return a namedtuple with variance data"""
        if not os.path.isfile(self.variance_file):
            raise RuntimeError,"Cannot find theseus variance file: {0} Please check the log: {1}".format(self.variance_file, self.theseus_log)
        var_by_res = self.parse_variances(self.variance_file)
        if core:
            return [v for v in var_by_res if v.core]
        else:
            return var_by_res

