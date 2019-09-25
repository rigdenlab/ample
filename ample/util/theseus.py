'''
Created on 26 May 2015

@author: jmht
'''
import collections
import logging
import os
import shutil

from ample.util import ample_util
from ample.util import sequence_util

# We create this here otherwise it causes problems with pickling
TheseusVariances = collections.namedtuple(
    'TheseusVariances', ['idx', 'resName', 'resSeq', 'variance', 'stdDev', 'rmsd', 'core']
)

_logger = logging.getLogger(__name__)


class Theseus(object):
    """Class to run THESEUS to superpose pdb files and determine the per-residue variances.
    
    .. _THESEUS website:
        http://www.theseus3d.org/
    """

    def __init__(self, work_dir=None, theseus_exe=None):

        if theseus_exe is None:
            if 'CCP4' in os.environ:
                theseus_exe = os.path.join(os.environ['CCP4'], 'bin', 'theseus')
        if theseus_exe is None or not os.path.exists(theseus_exe) and os.access(theseus_exe, os.X_OK):
            raise RuntimeError("Cannot find theseus_exe: {0}".format(theseus_exe))
        self.theseus_exe = theseus_exe
        self.work_dir = None
        self.var_by_res = None
        self.variance_log = None
        self.variance_log_test = None  # For mocking up a test
        self.superposed_models = None
        self.aligned_models = None
        self._set_work_dir(work_dir)
        return

    def _set_work_dir(self, work_dir):
        if work_dir:
            self.work_dir = work_dir
        if not self.work_dir:
            self.work_dir = os.getcwd()
        if not os.path.isdir(self.work_dir):
            os.mkdir(self.work_dir)
        return self.work_dir

    def alignment_file(self, models, alignment_file=None):
        """Create an alignment file for the models - this is based on the assumption they are all the same length
        but may have different residues"""
        if not alignment_file:
            alignment_file = os.path.join(self.work_dir, 'homologs.fasta')
        all_seq = sequence_util.Sequence(pdb=models[0])
        for model in models[1:]:
            all_seq += sequence_util.Sequence(pdb=model)
        if not all(map(lambda x: x == len(all_seq.sequences[0]), [len(s) for s in all_seq.sequences])):
            raise RuntimeError('PDB files are not all of the same length!\n{0}'.format(models))
        all_seq.write_fasta(alignment_file, pdbname=True)
        return alignment_file

    def superpose_models(self, models, work_dir=None, basename='theseus', homologs=False, alignment_file=None):
        """Superpose models and return the ensemble. Also set superposed_models and var_by_res variables.
        
        This also sets the `superposed_models` and `var_by_res` parameters.

        Parameters
        ----------
        models : :obj:`list`
            List of pdb files to be superposed.
        work_dir: str
            The directory to run theseus in and generate all the output files
        basename : str
            The stem that will be used to name all files
        homologs : bool
            True if the pdbs are homologous models as opposed to ab initio ones
        alignment_file : str
            An externally generated alignment file for homolgous models in FASTA format
            
        Returns
        -------
        superposed_models : a pdb file containing an ensemble of the superposed models 
        
        """
        self._set_work_dir(work_dir)
        if homologs:
            # Theseus expects all the models to be in the directory that it is run in as the string given in
            # the fasta header is used to construct the file names of the aligned pdb files. If a full or
            # relative path is given (e.g. /foo/bar.pdb), it tries to create files called "basename_/foo/bar.pdb"
            # We therefore copy the models in and then delete them afterwards
            if not alignment_file:
                alignment_file = self.alignment_file(models)
            copy_models = [os.path.join(self.work_dir, os.path.basename(m)) for m in models]
            for orig, copy in zip(models, copy_models):
                shutil.copy(orig, copy)
            models = copy_models

        # -Z included so we don't line the models up to the principle axis and -o so that they all line
        # up with the first model
        # cmd = [ self.theseus_exe, '-a0', '-r', basename ]
        cmd = [self.theseus_exe, '-a0', '-r', basename, '-Z', '-o', os.path.basename(models[0])]
        if homologs:
            cmd += ['-A', alignment_file]
            cmd += [os.path.basename(m) for m in models]
        else:
            # Not sure why we had relpath - fails some of the tests so changing
            # cmd += [ os.path.relpath(m,self.work_dir) for m in models ]
            cmd += models

        self.theseus_log = os.path.join(self.work_dir, "tlog_{0}.log".format(basename))
        retcode = ample_util.run_command(cmd, logfile=self.theseus_log, directory=self.work_dir)
        if retcode != 0:
            raise RuntimeError(
                "non-zero return code for theseus in superpose_models!\n See log: {0}".format(self.theseus_log)
            )

        self.variance_log = os.path.join(self.work_dir, '{0}_variances.txt'.format(basename))
        self.superposed_models = os.path.join(self.work_dir, '{0}_sup.pdb'.format(basename))
        if homologs:
            # Horrible - need to rename the models so that they match the names in the alignment file
            self.aligned_models = []
            for m in copy_models:
                mb = os.path.basename(m)
                aligned_model = os.path.join(self.work_dir, "{0}_{1}".format(basename, mb))
                os.unlink(m)
                os.rename(aligned_model, os.path.join(self.work_dir, mb))
                self.aligned_models.append(mb)

        # Set the variances
        self.var_by_res = self.parse_variances()
        return self.superposed_models

    def parse_variances(self):
        # The variance_log_test variable may be set if we are mocking out the variances for testing
        if self.variance_log_test:
            variance_log = self.variance_log_test
        else:
            variance_log = self.variance_log

        if not os.path.isfile(variance_log):
            raise RuntimeError("Cannot find theseus variance log: {0}".format(variance_log))
        var_by_res = []
        with open(variance_log) as f:
            for i, line in enumerate(f):

                if i == 0:
                    continue

                line = line.strip()
                if not line:
                    continue

                tokens = line.split()
                # Different versions of theseus may have a RES card first, so need to check
                if tokens[0] == "RES":
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
                isCore = False
                if len(tokens) > idxCore and tokens[idxCore] == 'CORE':
                    isCore = True
                var_by_res.append(
                    TheseusVariances(
                        idx=int(tokens[idxidx]) - 1,  # Theseus counts from 1, we count from 0,
                        resName=tokens[idxResName],
                        resSeq=int(tokens[idxResSeq]),
                        variance=float(tokens[idxVariance]),
                        stdDev=float(tokens[idxStdDev]),
                        rmsd=float(tokens[idxRmsd]),
                        core=isCore,
                    )
                )
        return var_by_res
