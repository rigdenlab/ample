#!/usr/bin/env ccp4-python

from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "28 Jul 2016"
__version__ = 2.0

import itertools
import logging
import operator
import os
import random
import string
import sys
import warnings

from ample.parsers import alignment_parser
from ample.parsers import tm_parser
from ample.util import ample_util
from ample.util import pdb_edit

try:
    from Bio import PDB
    from Bio import SeqIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


TMSCORE_EXE = None

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
        
    def _to_dict(self):
        """Convert the object to a dictionary"""
        dictionary = {'model_name': self.model_name, 'structure_name': self.structure_name,
                      'model_fname': self.model_fname, 'structure_fname': self.structure_fname,
                      'log_fname': self.log_fname, 'tmscore': self.tmscore, 'rmsd': self.rmsd,
                      'nr_residues_common': self.nr_residues_common, 'gdtts': self.gdtts,
                      'gdtha': self.gdtha, 'maxsub': self.maxsub, 'seq_id': self.seq_id}
        return dictionary


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

    def __init__(self, executable, method, wdir=None):
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
        self.method = method.lower()
        self.work_dir = wdir if wdir else os.getcwd()

        # Quick check if the executable exists
        if ample_util.is_exe(executable):
            self.executable = executable
        else:
            try:
                self.executable = ample_util.find_exe(executable)
            except:
                msg = "Cannot find provided executable: {0}".format(executable)
                logger.critical(msg)
                raise RuntimeError(msg)

        return

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
            msg = "Unequal number of models and structures"
            logger.critical(msg)
            raise RuntimeError(msg)

        # Create a logfile parser
        if self.method == "tmalign":
            pt = tm_parser.TMalignLogParser()
        elif self.method == "tmscore":
            pt = tm_parser.TMscoreLogParser()
        else:
            msg = "Invalid method selected: ", self.method
            logger.critical(msg)
            raise RuntimeError(msg)

        # =======================================================================
        # Iterate through the structure files and execute the TMscore comparisons
        # =======================================================================

        logger.info('Using algorithm: {0}'.format(self.method))
        logger.info('-------Evaluating decoys/models-------')
        entries = []

        for model_pdb, structure_pdb in zip(models, structures):

            model_name = os.path.splitext(os.path.basename(model_pdb))[0]
            structure_name = os.path.splitext(os.path.basename(structure_pdb))[0]

            logger.debug("Working on: {0} - {1}".format(model_name, structure_name))

            if not os.path.isfile(model_pdb):
                logger.warning("Cannot find: {0}".format(model_pdb))
                continue
            elif not os.path.isfile(structure_pdb):
                logger.warning("Cannot find: {0}".format(structure_pdb))
                continue

            # TODO: Spawn the jobs across a number of CPUs. ample_util.workers_util.run_scripts() maybe?
            log_name = "{0}_{1}_{2}.log".format(model_name, structure_name, self.method)
            log = os.path.join(self.work_dir, log_name)
            cmd = [self.executable, model_pdb, structure_pdb]
            ample_util.run_command(cmd, logfile=log, directory=self.work_dir)

            try:
                # Reset the TM log parser to default values
                pt.reset()
                # Parse the TM method logfile to extract the data
                pt.parse(log)
            except Exception:
                msg = "Issues processing the {0} log file: {1}".format(self.method, log)
                logger.critical(msg)
                log = "None"

            _entry = self._store(model_name, structure_name, model_pdb, structure_pdb, log, pt)
            entries.append(_entry)

        self.entries = entries
        return entries

    def _get_iterator(self, all_vs_all):
        """

        Parameters
        ----------
        all_vs_all: bool

        Returns
        -------
        iterator : function

        """
        # Use different itertools functions depending on the comparison type
        if all_vs_all:
            logger.info("All-vs-all comparison of models and structures")
            iterator = itertools.product     # yields an iterator of all unique combinations
        else:
            logger.info("Direct comparison of models and structures")
            iterator = itertools.izip        # yields a zipped iterator
        return iterator

    def _store(self, model_name, structure_name, model_pdb, structure_pdb, logfile, pt):
        # Generic data that either both parsers have or are defined independently of the parser
        model = ModelData(model_name, structure_name, model_pdb, structure_pdb, logfile, pt.tm, pt.rmsd, pt.nr_residues_common)
        # Specific attributes that either but not both parsers have
        if hasattr(pt, 'gdtts'):
            model.gdtts = pt.gdtts
        if hasattr(pt, 'gdtha'):
            model.gdtha = pt.gdtha
        if hasattr(pt, 'maxsub'):
            model.maxsub = pt.maxsub
        if hasattr(pt, 'seq_id'):
            model.seq_id = pt.seq_id

        return model._to_dict()

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
            TMSCORE_EXE = ample_util.find_exe(exe_name)
        except:
            TMSCORE_EXE = None
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
    >>> entries = tm.compare_to_structure(models, references)
    """

    def __init__(self, executable, wdir=None):
        super(TMalign, self).__init__(executable, "TMalign", wdir=wdir)

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
        # Check what we are comparing
        if len(structures) == 1:
            logger.info('Using single structure provided for all model comparisons')
            structures = [structures[0] for _ in xrange(len(models))]

        # The models parsed forward to the comparison
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
    >>> entries = tm.compare_to_structure(models, references)
    """

    def __init__(self, executable, wdir=None):
        super(TMscore, self).__init__(executable, "TMscore", wdir=wdir)

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
        if not BIOPYTHON_AVAILABLE: raise RuntimeError("Biopython is not available")
        # Check what we are comparing
        if structures and len(structures) == 1:
            logger.info('Using single structure provided for all model comparisons')
            structures = [structures[0] for _ in xrange(len(models))]

        if fastas and len(fastas) == 1:
            logger.info('Using single FASTA provided for all model comparisons')
            fastas = [fastas[0] for _ in xrange(len(models))]

        # The models parsed forward to the comparison
        models_to_compare, structures_to_compare = [], []
        combination_iterator = self._get_iterator(all_vs_all)

        if fastas:

            # Determine the iterator and create the combinations to be compared
            for index, (model, structure, fasta) in enumerate(combination_iterator(models, structures, fastas)):

                # Extract the FASTA sequence
                fasta_record = list(SeqIO.parse(open(fasta, 'r'), 'fasta'))[0]
                fasta_data = [(i + 1, j) for i, j in enumerate(str(fasta_record.seq))]

                # Extract some information from each PDB structure file
                model_data = list(self._pdb_info(model))
                structure_data = list(self._pdb_info(structure))

                # Align the model sequence to the FASTA sequence
                for fasta_pos in fasta_data:
                    if not fasta_pos in model_data:
                        model_data.append((fasta_pos[0], "-"))
                model_data.sort(key=operator.itemgetter(0))

                # Align the structure_sequence to the model sequence
                aln_parser = alignment_parser.AlignmentParser()
                fasta_structure_aln = aln_parser.align_sequences("".join(zip(*fasta_data)[1]),
                                                                 "".join(zip(*structure_data)[1]))

                # Remove parts of the alignment that are gaps in both sequences
                to_remove = []
                _alignment = zip("".join(zip(*model_data)[1]), fasta_structure_aln[1])
                for index, (model_res, structure_res) in enumerate(_alignment):
                    if model_res == "-" and structure_res == "-":
                        to_remove.append(index)
                for i in reversed(to_remove):
                    _alignment.pop(i)

                model_aln = "".join(zip(*_alignment)[0])
                structure_aln = "".join(zip(*_alignment)[1])

                if len(model_aln) != len(structure_aln):
                    msg = "Unequal lengths of your model and structure sequences"
                    logger.critical(msg)
                    raise RuntimeError(msg)

                # Modify the structures based on the aligned sequences
                pdb_combo = self._mod_structures(model_aln, structure_aln, model, structure)

                # Save the models
                models_to_compare.append(pdb_combo[0])
                structures_to_compare.append(pdb_combo[1])

        else:
            # No FASTA processing etc, pure comparisons of sequences
            for index, (model, structure) in enumerate(combination_iterator(models, structures)):

                # Extract some information from each PDB structure file
                model_data = list(self._pdb_info(model))
                structure_data = list(self._pdb_info(structure))

                # Align the sequences to see how much of the predicted decoys are in the xtal
                alignment = alignment_parser.AlignmentParser().align_sequences("".join(zip(*model_data)[1]),
                                                                               "".join(zip(*structure_data)[1]))
                # Redundant but identical to if
                alignment = zip(alignment[0], alignment[1])

                model_aln = "".join(zip(*alignment)[0])
                structure_aln = "".join(zip(*alignment)[1])

                if len(model_aln) != len(structure_aln):
                    msg = "Unequal lengths of your model and structure sequences"
                    logger.critical(msg)
                    raise RuntimeError(msg)

                # Modify the structures based on the aligned sequences
                pdb_combo = self._mod_structures(model_aln, structure_aln, model, structure)

                # Save the models9
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

        # ================================
        # File definitions
        # ================================

        # Create a storage for the files
        work_dir_mod = os.path.join(self.work_dir, "tm_util_pdbs")
        if not os.path.isdir(work_dir_mod):
            os.mkdir(work_dir_mod)

        # Create a random file suffix to avoid overwriting file names if duplicate
        # Taken from http://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python
        random_suffix = ''.join(random.SystemRandom().choice(string.ascii_lowercase + string.digits) for _ in range(10))

        # File names and output files
        model_name = os.path.basename(model_pdb).rsplit(".", 1)[0]
        model_pdb_ret = os.path.join(work_dir_mod, "_".join([model_name, random_suffix, "mod.pdb"]))

        structure_name = os.path.basename(structure_pdb).rsplit(".", 1)[0]
        structure_pdb_ret = os.path.join(work_dir_mod, "_".join([structure_name, random_suffix, "mod.pdb"]))

        # Check if the files we are to create for comparison do not exist
        if os.path.isfile(model_pdb_ret) or os.path.isfile(structure_pdb_ret):
            msg = "Comparison structures exist. Move, delete or rename before continuing"
            logger.critical(msg)
            raise RuntimeError(msg)

        # Create temporary files
        _model_pdb_tmp_stage1 = ample_util.tmp_file_name(delete=False, directory=work_dir_mod, suffix=".pdb")
        _model_pdb_tmp_stage2 = ample_util.tmp_file_name(delete=False, directory=work_dir_mod, suffix=".pdb")

        _structure_pdb_tmp_stage1 = ample_util.tmp_file_name(delete=False, directory=work_dir_mod, suffix=".pdb")
        _structure_pdb_tmp_stage2 = ample_util.tmp_file_name(delete=False, directory=work_dir_mod, suffix=".pdb")

        # ==================================
        # File manipulation and modification
        # ==================================

        # Get the gap positions in both sequences
        model_gaps = self._find_gaps(model_aln)
        structure_gaps = self._find_gaps(structure_aln)

        # Renumber the pdb files - required in case there are any gaps
        pdb_edit.renumber_residues_gaps(model_pdb, _model_pdb_tmp_stage1, model_gaps)
        pdb_edit.renumber_residues_gaps(structure_pdb, _structure_pdb_tmp_stage1, structure_gaps)

        # Determine the gap indeces
        model_gaps_indeces = [i+1 for i, is_gap in enumerate(model_gaps) if is_gap]
        structure_gaps_indeces = [i + 1 for i, is_gap in enumerate(structure_gaps) if is_gap]

        # Use gaps of other sequence to even out
        pdb_edit.select_residues(_model_pdb_tmp_stage1, _model_pdb_tmp_stage2, delete=structure_gaps_indeces)
        pdb_edit.select_residues(_structure_pdb_tmp_stage1, _structure_pdb_tmp_stage2, delete=model_gaps_indeces)

        # Renumber the pdb files - required by TMscore binary
        pdb_edit.renumber_residues(_model_pdb_tmp_stage2, model_pdb_ret)
        pdb_edit.renumber_residues(_structure_pdb_tmp_stage2, structure_pdb_ret)

        # ==================================
        # Checks and validations
        # ==================================

        # Extract some information from each PDB structure file
        _model_data = list(self._pdb_info(model_pdb_ret))
        _structure_data = list(self._pdb_info(structure_pdb_ret))

        # Make sure our structures contain the same residues with correct indeces
        if set(_model_data) != set(_structure_data):
            msg = "Structure file modification did not work. Affected PDBs {0} - {1}".format(model_name, structure_name)
            logger.critical(msg)
            raise RuntimeError(msg)

        # Remove the temporary files
        for f in [_model_pdb_tmp_stage1, _model_pdb_tmp_stage2, _structure_pdb_tmp_stage1, _structure_pdb_tmp_stage2]:
            os.unlink(f)

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
        indeces : list
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
        structure = PDB.PDBParser().get_structure("pdb", pdb)

        # Weird problem with get_...() methods if MODEL is not explicitly stated in
        # PDB file. Thus, just iterate over everything return when top chain completed
        for model in structure:
            for index, chain in enumerate(model):
                if index > 0: return
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
        index : int
           Residue sequence index of first residue in structure file
        """
        for line in open(pdb, 'r'):
            if line.startswith("ATOM"):
                line = line.split()
                return int(line[5])


def tm_available():
    """
    Check if TM binary is available and set TMSCORE_EXE module variable.

    Returns
    -------
    bool
    """
    warnings.warn("This method will be deprecated in future releases. Use `TMapps.binary_avail()` instead.", DeprecationWarning)

    exe_name = "TMscore" + ample_util.EXE_EXT
    try:
        TMSCORE_EXE = ample_util.find_exe(exe_name)
    except:
        TMSCORE_EXE = None
        return False
    return True


def main():
    import argparse
    import shutil
    parser = argparse.ArgumentParser()
    parser.add_argument('--allvall', action="store_true", help="All vs all comparison")
    parser.add_argument('--log', default='info', choices=['debug', 'info', 'warning', 'error'],
                        help="logging level (defaults to 'warning')")
    # parser.add_argument('--purge', action="store_true", help="Remove the run directory")
    parser.add_argument('--rundir', default=os.getcwd(), help="Run directory")
    parser.add_argument('-f', '-fasta', dest="fastas", nargs="+", default=None)
    parser.add_argument('-m', '-models', dest="models", nargs="+", required=True)
    parser.add_argument('-s', '-structures', dest="structures", nargs="+", required=True)
    tmapp = parser.add_mutually_exclusive_group()
    tmapp.add_argument('-tm', '-tmscore', dest="tmscore", type=str)
    tmapp.add_argument('-ta', '-tmalign', dest="tmalign", type=str)
    args = parser.parse_args()

    # Set up some very basic logging - Logging taken from
    # http://stackoverflow.com/questions/30824981/do-i-need-to-explicitly-check-for-name-main-before-calling-getlogge
    logging.basicConfig(level=getattr(logging, args.log.upper(), None), format='%(levelname)s: %(message)s')
    
    # Make the run directory if it doesn't exist
    if not os.path.isdir(args.rundir):
        os.mkdir(args.rundir)
    
    # Perform the comparison
    if args.tmalign:
        entries = TMalign(args.tmalign, wdir=args.rundir).compare_structures(args.models, args.structures, all_vs_all=args.allvall)
    elif args.tmscore:
        entries = TMscore(args.tmscore, wdir=args.rundir).compare_structures(args.models, args.structures, fastas=args.fastas, all_vs_all=args.allvall)
    else:
        entries = None
    
    # Remove run directory if required
    # if args.purge:
    #     shutil.rmtree(os.path.join(args.rundir))

    # Do a much fancier table print statement if pandas is installed
    # TODO: Print statement from entries dictionary if we don't have pandas
    try:
        import pandas
        PANDAS_AVAILABLE = True
    except ImportError:
        import csv
        PANDAS_AVAILABLE = False
    
    csv_file = os.path.join(args.rundir, 'tm_results.csv')
    if PANDAS_AVAILABLE and entries:
        df = pandas.pandas.DataFrame(data=entries)
        df.sort_values("tmscore", inplace=True, ascending=False)
        df.reset_index(drop=True)
        logger.info("Results table:\n{0}".format(df[["model_name", "structure_name", "tmscore"]].to_string()))
        df.to_csv(csv_file)
    elif entries:
        with open(csv_file, 'wb') as f_out:
            w = csv.DictWriter(f_out, entries[0].keys())
            w.writeheader()
            for entry in entries:
                w.writerow(entry)

    return entries

if __name__ == "__main__":
    main()
    sys.exit(0)
