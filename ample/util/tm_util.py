#!/usr/bin/env ccp4-python

from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "28 Jul 2016"
__version__ = 1.0

import itertools
import logging
import operator
import os
import random
import string
import sys
import tempfile
import warnings

from ample.parsers import alignment_parser
from ample.parsers import tm_parser
from ample.util import ample_util
from ample.util import pdb_edit

from pyjob import Job
from pyjob.misc import make_script

try:
    from Bio import PDB
    from Bio import SeqIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

logger = logging.getLogger(__name__)


class ModelData(object):
    """Class to store model data"""
    __slots__ = ('model_name', 'structure_name', 'model_fname', 'structure_fname', 
                 'log_fname', 'tmscore', 'rmsd', 'nr_residues_common', 'gdtts', 
                 'gdtha', 'maxsub', 'seq_id')
    
    def __init__(self, model_name, structure_name, model_fname, structure_fname, log_fname, tmscore, rmsd, nr_residues_common):
        self.model_name = model_name
        self.structure_name = structure_name
        self.model_fname = model_fname
        self.structure_fname = structure_fname
        self.log_fname = log_fname
        self.tmscore = tmscore
        self.rmsd = rmsd
        self.nr_residues_common = nr_residues_common
        self.gdtts = 0.0
        self.gdtha = 0.0
        self.maxsub = 0.0
        self.seq_id = 0

    def _asdict(self):
        """Convert the object to a dictionary"""
        return {k: getattr(self, k) for k in self.__slots__}


class TMapps(object):
    """
    Superclass of TMscore and TMalign

    Attributes
    ----------
    entries : list
       List containing the TMscore entries on a per-model basis
    structure : str
       Path to the reference structure
    method : str
       One of [TMscore|TMalign]
    executable : str
       Path to the TMscore executable
    work_dir : str
       Path to the working directory

    Todo
    ----
    * Function to return the entries as numpy matrix
    * Function to return the entries in a pandas dataframe

    """

    def __init__(self, executable, method, wdir=".", **kwargs):
        """
        Parameters
        ----------
        executable : str
           Path to the executable
        method : str
           One of [TMscore|TMalign]
        work_dir : str
           Path to the working directory

        """
        self.entries = []
        self.executable = executable
        self.method = method.lower()
        self.tmp_dir = None
        self.work_dir = os.path.abspath(wdir)

        if 'submit_cluster' in kwargs and kwargs['submit_cluster']:
            self._qtype = kwargs['submit_qtype']
        else:
            self._qtype = "local"
        self._queue = kwargs['submit_queue'] if 'submit_queue' in kwargs else None
        self._nproc = kwargs['nproc'] if 'nproc' in kwargs else 1
        self._max_array_jobs = kwargs['submit_max_array'] if 'submit_max_array' in kwargs else None

    def comparison(self, models, structures):
        """
        Compare a list of model structures to a second list of reference structures

        Parameters
        ----------
        models : list
           List containing the paths to the model structure files
        structures : list
           List containing the paths to the reference structure files

        Returns
        -------
        entries : list
           List of TMscore data entries on a per-model basis

        """

        if len(models) < 1 or len(structures) < 1:
            msg = 'No model structures provided' if len(models) < 1 else \
                'No reference structures provided'
            logger.critical(msg)
            raise RuntimeError(msg)

        elif len(structures) == 1:
            logger.info('Using single structure provided for all model comparisons')
            structures = [structures[0] for _ in xrange(len(models))]

        elif len(models) != len(structures):
            msg = "Unequal number of models and structures!"
            logger.critical(msg)
            raise RuntimeError(msg)

        if self.method == "tmalign":
            pt = tm_parser.TMalignLogParser()
        elif self.method == "tmscore":
            pt = tm_parser.TMscoreLogParser()
        else:
            msg = "Invalid method selected: %s", self.method
            logger.critical(msg)
            raise RuntimeError(msg)

        logger.info('Using algorithm: {0}'.format(self.method))
        logger.info('------- Evaluating decoys -------')
        data_entries, job_scripts, log_files = [], [], []
        for model_pdb, structure_pdb in zip(models, structures):
            model_name = os.path.splitext(os.path.basename(model_pdb))[0]
            structure_name = os.path.splitext(os.path.basename(structure_pdb))[0]
            stem = "_".join([model_name, structure_name, self.method])

            if os.path.isfile(model_pdb) and os.path.isfile(structure_pdb):
                data_entries.append([model_name, structure_name, model_pdb, structure_pdb])
                script = make_script([self.executable, model_pdb, structure_pdb], prefix="tmscore_", stem=stem)
                job_scripts.append(script)
                log_files.append(os.path.splitext(script)[0] + ".log")
            else:
                if not os.path.isfile(model_pdb):
                    logger.warning("Cannot find: %s", model_pdb)
                if not os.path.isfile(structure_pdb):
                    logger.warning("Cannot find: %s", structure_pdb))
                continue
            
        logger.info('Executing TManalysis scripts')
        j = Job(self._qtype)
        j.submit(job_scripts, nproc=self._nproc, max_array_jobs=self._max_array_jobs, 
                 queue=self._queue, name="tmscore")
        j.wait(interval=1)

        self.entries = []
        for entry, log, script in zip(data_entries, log_files, job_scripts):
            try:
                pt.reset()
                pt.parse(log)
            except Exception:
                logger.critical("Error processing the %s log file: %s", self.method, log)
                log = "None"
            model_name, structure_name, model_pdb, structure_pdb = entry
            _entry = self._store(model_name, structure_name, model_pdb, structure_pdb, log, pt)
            self.entries.append(_entry)
            os.unlink(script)

        return self.entries

    def _get_iterator(self, all_vs_all):
        """

        Parameters
        ----------
        all_vs_all: bool

        Returns
        -------
        function

        """
        if all_vs_all:
            logger.info("All-vs-all comparison of models and structures")
            return itertools.product
        else:
            logger.info("Direct comparison of models and structures")
            return itertools.izip

    def _store(self, model_name, structure_name, model_pdb, structure_pdb, logfile, pt):
        model = ModelData(
            model_name, structure_name, model_pdb, structure_pdb, 
            logfile, pt.tm, pt.rmsd, pt.nr_residues_common
        )
        if hasattr(pt, 'gdtts'):
            model.gdtts = pt.gdtts
        if hasattr(pt, 'gdtha'):
            model.gdtha = pt.gdtha
        if hasattr(pt, 'maxsub'):
            model.maxsub = pt.maxsub
        if hasattr(pt, 'seq_id'):
            model.seq_id = pt.seq_id
        return model._asdict()

    @staticmethod
    def binary_avail(binary):
        """Check if TM binary is available

        Paramaters
        ----------
        binary : str
           The binary name of `TMscore` or `TMalign`

        Returns
        -------
        bool
        
        Raises
        ------
        ValueError
           The binary is not `TMalign` or `TMscore`

        """
        if binary.lower() == 'tmalign':
            exe_name = "TMalign" + ample_util.EXE_EXT
        elif binary.lower() == 'tmscore':
            exe_name = "TMscore" + ample_util.EXE_EXT
        else:
            raise ValueError('Provide one of TMalign or TMscore')
        try:
            ample_util.find_exe(exe_name)
        except:
            return False
        return True


class TMalign(TMapps):
    """
    Wrapper to handle TMalign scoring for one or more structures

    Examples
    --------
    >>> models = ["<MODEL_1>", "<MODEL_2>", "<MODEL_3>"]
    >>> references = ["<REFERENCE_1>", "<REFERENCE>", "<REFERENCE>"]
    >>> tm = TMalign("<PATH_TO_EXE>")
    >>> entries = tm.compare_structures(models, references)

    """

    def __init__(self, executable, wdir=None, **kwargs):
        super(TMalign, self).__init__(executable, "TMalign", wdir=wdir, **kwargs)

    def compare_structures(self, models, structures, all_vs_all=False):
        """
        Compare a list of model structures to a second list of reference structures

        Parameters
        ----------
        models : list
           List containing the paths to the model structure files
        structures : list
           List containing the paths to the reference structure files
        all_vs_all : bool
           Flag to compare all models against all structures

        Returns
        -------
        entries : list

        """
        if len(structures) == 1:
            logger.info('Using single structure provided for all model comparisons')
            structures = [structures[0] for _ in xrange(len(models))]

        models_to_compare, structures_to_compare = [], []
        combination_iterator = self._get_iterator(all_vs_all)

        for (model, structure) in combination_iterator(models, structures):
            models_to_compare.append(model)
            structures_to_compare.append(structure)

        return self.comparison(models_to_compare, structures_to_compare)
        

class TMscore(TMapps):
    """
    Wrapper to handle TMscoring for one or more structures

    Examples
    --------
    >>> models = ["<MODEL_1>", "<MODEL_2>", "<MODEL_3>"]
    >>> references = ["<REFERENCE_1>", "<REFERENCE>", "<REFERENCE>"]
    >>> tm = TMscore("<PATH_TO_EXE>")
    >>> entries = tm.compare_structures(models, references)

    """

    def __init__(self, executable, wdir=None, **kwargs):
        super(TMscore, self).__init__(executable, "TMscore", wdir=wdir, **kwargs)
        self.tmp_dir = os.path.join(self.work_dir, "tm_util_pdbs")
        if not os.path.isdir(self.tmp_dir):
            os.mkdir(self.tmp_dir)

    def compare_structures(self, models, structures, fastas=None, all_vs_all=False):
        """
        Compare a list of model structures to a second list of reference structures

        Parameters
        ----------
        models : list
           List containing the paths to the model structure files
        structures : list
           List containing the paths to the reference structure files
        fastas : list
           List containing the paths to the FASTA files
        all_vs_all : bool
           Flag to compare all models against all structures

        Returns
        -------
        entries : list
           List of TMscore data entries on a per-model basis

        Notes
        -----
        If a FASTA sequence is provided, a much more accurate comparison can be carried out. However, to by-pass this
        there is also an option to run the comparison without it. This might work just fine for larger models.

        """
        if not BIOPYTHON_AVAILABLE:
            raise RuntimeError("Biopython is not available")
        if structures and len(structures) == 1:
            logger.info('Using single structure provided for all model comparisons')
            structures = [structures[0] for _ in xrange(len(models))]

        if fastas and len(fastas) == 1:
            logger.info('Using single FASTA provided for all model comparisons')
            fastas = [fastas[0] for _ in xrange(len(models))]

        models_to_compare, structures_to_compare = [], []
        combination_iterator = self._get_iterator(all_vs_all)

        if fastas:

            for (model, structure, fasta) in combination_iterator(models, structures, fastas):
                fasta_record = list(SeqIO.parse(open(fasta, 'r'), 'fasta'))[0]
                fasta_data = [(i + 1, j) for i, j in enumerate(str(fasta_record.seq))]

                model_data = list(self._pdb_info(model))
                structure_data = list(self._pdb_info(structure))

                for fasta_pos in fasta_data:
                    if fasta_pos not in model_data:
                        model_data.append((fasta_pos[0], "-"))
                model_data.sort(key=operator.itemgetter(0))

                aln_parser = alignment_parser.AlignmentParser()
                fasta_structure_aln = aln_parser.align_sequences("".join(zip(*fasta_data)[1]),
                                                                 "".join(zip(*structure_data)[1]))

                to_remove = []
                _alignment = zip("".join(zip(*model_data)[1]), fasta_structure_aln[1])
                for i, (model_res, structure_res) in enumerate(_alignment):
                    if model_res == "-" and structure_res == "-":
                        to_remove.append(i)
                for i in reversed(to_remove):
                    _alignment.pop(i)

                model_aln = "".join(zip(*_alignment)[0])
                structure_aln = "".join(zip(*_alignment)[1])

                if len(model_aln) != len(structure_aln):
                    msg = "Unequal lengths of your model and structure sequences"
                    logger.critical(msg)
                    raise RuntimeError(msg)

                pdb_combo = self._mod_structures(model_aln, structure_aln, model, structure)

                models_to_compare.append(pdb_combo[0])
                structures_to_compare.append(pdb_combo[1])

        else:
            for (model, structure) in combination_iterator(models, structures):

                model_data = list(self._pdb_info(model))
                structure_data = list(self._pdb_info(structure))

                alignment = alignment_parser.AlignmentParser().align_sequences("".join(zip(*model_data)[1]),
                                                                               "".join(zip(*structure_data)[1]))
                alignment = zip(alignment[0], alignment[1])

                model_aln = "".join(zip(*alignment)[0])
                structure_aln = "".join(zip(*alignment)[1])

                if len(model_aln) != len(structure_aln):
                    msg = "Unequal lengths of your model and structure sequences"
                    logger.critical(msg)
                    raise RuntimeError(msg)

                pdb_combo = self._mod_structures(model_aln, structure_aln, model, structure)

                models_to_compare.append(pdb_combo[0])
                structures_to_compare.append(pdb_combo[1])

        return self.comparison(models_to_compare, structures_to_compare)

    def _mod_structures(self, model_aln, structure_aln, model_pdb, structure_pdb):
        """

        Parameters
        ----------
        model_aln : str
           A string containing the aligned sequence of the model
        structure_aln : str
           A string containing the alignment sequence of the structure
        model_pdb : str
           The path to the model pdb file
        structure_pdb : str
           The path to the structure pdb file

        Returns
        -------
        model_pdb_ret : str
           The path to the modified model pdb file
        structure_pdb_ret : str
           The path to the modified structure pdb file

        """
        random_suffix = ''.join(random.SystemRandom().choice(string.ascii_lowercase + string.digits) for _ in range(10))

        model_name = os.path.basename(model_pdb).rsplit(".", 1)[0]
        model_pdb_ret = os.path.join(self.tmp_dir, "_".join([model_name, random_suffix, "mod.pdb"]))

        structure_name = os.path.basename(structure_pdb).rsplit(".", 1)[0]
        structure_pdb_ret = os.path.join(self.tmp_dir, "_".join([structure_name, random_suffix, "mod.pdb"]))

        if os.path.isfile(model_pdb_ret) or os.path.isfile(structure_pdb_ret):
            msg = "Comparison structures exist. Move, delete or rename before continuing"
            logger.critical(msg)
            raise RuntimeError(msg)

        _model_pdb_tmp_stage1 = ample_util.tmp_file_name(delete=False, directory=self.tmp_dir, suffix=".pdb")
        _model_pdb_tmp_stage2 = ample_util.tmp_file_name(delete=False, directory=self.tmp_dir, suffix=".pdb")

        _structure_pdb_tmp_stage1 = ample_util.tmp_file_name(delete=False, directory=self.tmp_dir, suffix=".pdb")
        _structure_pdb_tmp_stage2 = ample_util.tmp_file_name(delete=False, directory=self.tmp_dir, suffix=".pdb")

        model_gaps = self._find_gaps(model_aln)
        structure_gaps = self._find_gaps(structure_aln)

        pdb_edit.renumber_residues_gaps(model_pdb, _model_pdb_tmp_stage1, model_gaps)
        pdb_edit.renumber_residues_gaps(structure_pdb, _structure_pdb_tmp_stage1, structure_gaps)

        model_gaps_indeces = [i+1 for i, is_gap in enumerate(model_gaps) if is_gap]
        structure_gaps_indeces = [i + 1 for i, is_gap in enumerate(structure_gaps) if is_gap]

        pdb_edit.select_residues(_model_pdb_tmp_stage1, _model_pdb_tmp_stage2, delete=structure_gaps_indeces)
        pdb_edit.select_residues(_structure_pdb_tmp_stage1, _structure_pdb_tmp_stage2, delete=model_gaps_indeces)

        pdb_edit.renumber_residues(_model_pdb_tmp_stage2, model_pdb_ret)
        pdb_edit.renumber_residues(_structure_pdb_tmp_stage2, structure_pdb_ret)

        _model_data = list(self._pdb_info(model_pdb_ret))
        _structure_data = list(self._pdb_info(structure_pdb_ret))

        if set(_model_data) != set(_structure_data):
            msg = "Residues in model and structure non-identical. Affected PDBs {0} - {1}".format(model_name, structure_name)
            logger.critical(msg)
            raise RuntimeError(msg)

        map(os.unlink, [
            _model_pdb_tmp_stage1, _model_pdb_tmp_stage2, _structure_pdb_tmp_stage1, _structure_pdb_tmp_stage2
        ])

        return model_pdb_ret, structure_pdb_ret

    def _find_gaps(self, seq):
        """
        Identify gaps in the protein chain

        Parameters
        ----------
        seq : str
           String of amino acids

        Returns
        -------
        indexes : list
           List of booleans that contain gaps

        """
        return [True if char == "-" else False for char in seq]

    def _pdb_info(self, pdb):
        """
        Obtain the pdb indeces and residue names

        Parameters
        ----------
        pdb : str
           The path to a PDB file

        Yields
        ------
        list
            A list containing per residue information

        """
        if not BIOPYTHON_AVAILABLE: raise RuntimeError("Biopython is not available")
        with warnings.catch_warnings():
            logger.debug("Suppressing BIOPYTHON warnings for PDBParser")
            warnings.simplefilter("ignore")
            structure = PDB.PDBParser().get_structure("pdb", pdb)

        # Weird problem with get_...() methods if MODEL is not explicitly stated in
        # PDB file. Thus, just iterate over everything return when top chain completed
        for model in structure:
            for i, chain in enumerate(model):
                if i > 0:
                    return
                for residue in chain:
                    hetero, res_seq, _ = residue.get_id()

                    if hetero.strip():
                        logger.debug("Hetero atom detected in {0}: {1}".format(pdb, res_seq))
                        continue

                    resname_three = residue.resname
                    if resname_three == "MSE":
                        resname_three = "MET"
                    resname_one = PDB.Polypeptide.three_to_one(resname_three)

                    yield (res_seq, resname_one)

    def _residue_one(self, pdb):
        """
        Find the first residue index in a pdb structure

        Parameters
        ----------
        pdb : str
           Path to a structure file in PDB format

        Returns
        -------
        int
           Residue sequence index of first residue in structure file

        """
        for line in open(pdb, 'r'):
            if line.startswith("ATOM"):
                return int(line.split()[5])


def main():
    import argparse
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument('--allvall', action="store_true", help="All vs all comparison")
    parser.add_argument('--log', default='info', choices=['debug', 'info', 'warning', 'error'],
                        help="logging level (defaults to 'warning')")
    parser.add_argument('--purge', action="store_true", help="Remove all temporary files")
    parser.add_argument('--rundir', default='.', help="Run directory")
    parser.add_argument('-f', '--fasta', dest="fastas", nargs="+", default=None)
    parser.add_argument('-m', '--models', dest="models", nargs="+", required=True)
    parser.add_argument('-s', '--structures', dest="structures", nargs="+", required=True)
    parser.add_argument('-t', '--threads', default=1, type=int)
    tmapp = parser.add_mutually_exclusive_group(required=True)
    tmapp.add_argument('-tm', '--tmscore', dest="tmscore", type=str)
    tmapp.add_argument('-ta', '--tmalign', dest="tmalign", type=str)
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log.upper(), None), format='%(levelname)s: %(message)s')
    
    if not os.path.isdir(args.rundir):
        os.mkdir(args.rundir)

    if args.fastas:
        args.fastas = [os.path.abspath(m) for m in args.fastas]
    args.models = [os.path.abspath(m) for m in args.models]
    args.structures = [os.path.abspath(m) for m in args.structures]
    
    kwargs = dict(wdir=args.rundir, nproc=args.threads)
    if args.tmalign:
        tmapp = TMalign(args.tmalign, **kwargs)
        tmapp.compare_structures(args.models, args.structures, all_vs_all=args.allvall)
    elif args.tmscore:
        tmapp = TMscore(args.tmscore, **kwargs)
        tmapp.compare_structures(args.models, args.structures, fastas=args.fastas, all_vs_all=args.allvall)
    else:
        raise RuntimeError("You're doomed if you get here!")

    df = pd.DataFrame(tmapp.entries) 
    df.to_csv(os.path.join(args.rundir, 'tm_results.csv'), index=False)
    logger.info("Final TMscore table:\n\n" + df.to_string(index=False) + "\n")

    if args.purge and tmapp.tmp_dir:
        from shutil import rmtree
        rmtree(tmapp.tmp_dir)

    return tmapp.entries

if __name__ == "__main__":
    main()
    sys.exit(0)
