#!/usr/bin/env ccp4-python

import filecmp
import itertools
import logging
import os

from ample.parsers import alignment_parser
from ample.parsers import tm_parser
from ample.util import ample_util
from ample.util import pdb_edit

__author__ = "Felix Simkovic"
__date__ = "09.11.2015"
__version__ = "2.0"

LOGGER = logging.getLogger(__name__)

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
        self.executable = executable
        self.method = method.lower()
        self.work_dir = wdir if wdir else os.getcwd()

        self.entries = []
        return

    def comparison(self, models, structures, all_vs_all=False, keep_modified_structures=False, identical_sequences=False):
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
        keep_modified_structures : bool
           Flag to delete intermediate, modified structure files
        identical_sequences : bool
           Flag to avoid any modification of files due to sequence identity

        Returns
        -------
        entries : list
           List of TMscore data entries on a per-model basis
        """

        if len(models) < 1 or len(structures) < 1:
            msg = 'No model structures provided' if len(models) < 1 else \
                'No reference structures provided'
            LOGGER.critical(msg)
            raise RuntimeError(msg)

        elif len(structures) == 1:
            LOGGER.info('Using single structure provided for all model comparisons')
            structures = [structures[0] for _ in xrange(len(models))]

        elif len(models) != len(structures):
            msg = "Unequal number of models and structures"
            LOGGER.critical(msg)
            raise RuntimeError(msg)

        # Use different itertools functions depending on the comparison type
        if all_vs_all:
            LOGGER.info("All-vs-all comparison of models and structures")
            iterator = itertools.product     # yields an iterator of all unique combinations
        else:
            LOGGER.info("Direct comparison of models and structures")
            iterator = itertools.izip        # yields a zipped iterator

        # Create a logfile parser
        if self.method == "tmalign":
            pt = tm_parser.TMalignLogParser()
        elif self.method == "tmscore":
            pt = tm_parser.TMscoreLogParser()
        else:
            msg = "Invalid method selected: ", self.method
            LOGGER.critical(msg)
            raise RuntimeError(msg)

        # =======================================================================
        # Iterate through the structure files and execute the TMscore comparisons
        # =======================================================================

        LOGGER.info('-------Evaluating decoys/models-------')
        LOGGER.info('Using algorithm: ', self.method)
        entries = []

        for model, structure in iterator(models, structures):

            # Quick file-check to check for similarity [skip identical files]
            identical_structures = self._identical_structures(model, structure)

            model_name = os.path.splitext(os.path.basename(model))[0]
            structure_name = os.path.splitext(os.path.basename(structure))[0]

            LOGGER.debug("Working on: {0} - {1}".format(model_name, structure_name))

            if not os.path.isfile(model):
                LOGGER.warning("Cannot find: {0}".format(model))
                continue
            elif not os.path.isfile(structure):
                LOGGER.warning("Cannot find: {0}".format(structure))
                continue

            # Modify structures to be identical as required by TMscore binary
            if not (identical_sequences and identical_structures and self.method != "tmalign"):
                model_mod = os.path.join(self.work_dir, model_name + "_mod.pdb")
                structure_mod = os.path.join(self.work_dir, model_name + "_" + structure_name + "_mod.pdb")
                self._mod_structures(model, model_mod, structure, structure_mod)
                model, structure = model_mod, structure_mod

            # TODO: Spawn the jobs across a number of CPUs. ample_util.workers_util.run_scripts() maybe?
            log = os.path.join(self.work_dir, model_name + "_{0}.log".format(self.method))
            cmd = [self.executable, model, structure]
            ample_util.run_command(cmd, logfile=log, directory=self.work_dir)

            # Delete the modified structures if not wanted
            if not (keep_modified_structures and identical_sequences and identical_structures and self.method != "tmalign"):
                os.remove(model_mod)
                os.remove(structure_mod)

            try:
                # Reset the TM log parser to default values
                pt.reset()
                # Parse the TM method logfile to extract the data
                pt.parse(log)
            except Exception:
                msg = "Issues processing the {0} log file: {1}".format(self.method, log)
                LOGGER.critical(msg)
                log = "None"

            _entry = self._store(model_name, model, log, structure, pt)
            entries.append(_entry)

        self.entries = entries
        return entries

    def _identical_structures(self, model, structure):
        """
        Description
        -----------
        Check for an all-vs-all comparison to not go through overhead of modifying
        identical structure files. Important when making large comparisons.
        Comparison itself important for statistics like nr_common_residues

        Arguments
        ---------
        model: str
           The path to a PDB file
        structure: str
           The path to a PDB file

        Returns
        -------
        identical_structures: bool

        """
        return True if filecmp.cmp(model, structure) else False

    def _mod_structures(self, model, model_mod, structure, structure_mod):
        """
        Modify the two structure files to match each other

        Description
        -----------
        Structure files often contain unequal residue numberings, missing residues in the
        chain or other mal-formatted data. This function aims to remove such discrepancies
        to allow for the most accurate comparisons possible.

        Parameters
        ----------
        model : str
           Path to the model pdb structure file
        model_mod : str
           Path to the modified model pdb structure file [does not need to exist]
        structure : str
           Path to the reference pdb structure file
        structure_mod : str
           Path to the modified reference pdb structure file [does not need to exist]
        """

        # Disable the info logger to not spam the user with which chain of native extracted.
        # Happens for every model + native below
        # http://stackoverflow.com/questions/2266646/how-to-i-disable-and-re-enable-console-logging-in-python
        logging.disable(logging.INFO)

        model_seq = pdb_edit.sequence(model).values()[0]
        structure_seq = pdb_edit.sequence(structure).values()[0]

        # Align the sequences to see how much of the predicted decoys are in the xtal
        aligned_seq_list = alignment_parser.AlignmentParser().align_sequences(model_seq,
                                                                              structure_seq)
        model_seq_ali = aligned_seq_list[0]
        structure_seq_ali = aligned_seq_list[1]

        # Get the gaps in both sequences
        model_gaps = self._find_gaps(model_seq_ali)
        structure_gaps = self._find_gaps(structure_seq_ali)

        ## STAGE 1 - REMOVE RESIDUES ##
        model_stage1 = ample_util.tmp_file_name(delete=False, directory=self.work_dir, suffix=".pdb")
        structure_stage1 = ample_util.tmp_file_name(delete=False, directory=self.work_dir, suffix=".pdb")

        # Get first residue number to adjust list of residues to remove
        model_res1 = self._residue_one(model)
        structure_res1 = self._residue_one(structure)

        # Match the residue lists to fit the residue 1 number
        model_gaps = [i + structure_res1 - 1 for i in model_gaps]
        structure_gaps = [i + model_res1 - 1 for i in structure_gaps]

        # Use gaps of other sequence to even out
        pdb_edit.select_residues(model, model_stage1, delete=structure_gaps)
        pdb_edit.select_residues(structure, structure_stage1, delete=model_gaps)

        ## STAGE 2 - RENUMBER RESIDUES ##
        pdb_edit.renumber_residues(model_stage1, model_mod)
        pdb_edit.renumber_residues(structure_stage1, structure_mod)

        os.unlink(model_stage1)
        os.unlink(structure_stage1)

        logging.disable(logging.NOTSET)

        return

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
           List of indices that contain gaps

        """
        return [i + 1 for i, c in enumerate(seq) if c == "-"]

    def _store(self, name, model, logfile, structure, pt):
        # Generic data that either both parsers have or are defined independently of the parser
        data_storage = {'name': name, 'model': model, 'TM_log': logfile, 'structure': structure,
                        'tmscore': pt.tm, 'rmsd': pt.rmsd, 'nr_residues_common': pt.nr_residues_common}

        # Specific attributes that either but not both parsers have
        if hasattr(pt, 'gdtts'):
            data_storage['gdtts'] = pt.gdtts
        if hasattr(pt, 'gdtha'):
            data_storage['gdtha'] = pt.gdtha
        if hasattr(pt, 'maxsub'):
            data_storage['maxsub'] = pt.maxsub
        if hasattr(pt, 'seq_id'):
            data_storage['seq_id'] = pt.seq_id

        return data_storage


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

    def compare_structures(self, *args, **kwargs):
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
        return self.comparison(*args, **kwargs)
        

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

    def compare_structures(self, *args, **kwargs):
        """
        Compare a list of model structures to a second list of reference structures

        Parameters
        ----------
        models : list
           List containing the paths to the model structure files
        structures : list
           List containing the paths to the reference structure files
        fasta : str
           The path to a FASTA file for improved alignment
        all_vs_all : bool
           Flag to compare all models against all structures
        keep_modified_structures : bool
           Flag to delete intermediate, modified structure files
        identical_sequences : bool
           Flag to avoid any modification of files due to sequence identity

        Returns
        -------
        entries : list
           List of TMscore data entries on a per-model basis
        """
        return self.comparison(*args, **kwargs)


def tm_available(app):
    """
    Check if TM binary is available

    Returns
    -------
    bool
    """
    try:
        ample_util.find_exe(app)
    except:
        return False
    return True

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--allvall', action="store_true", help="All vs all comparison")
    parser.add_argument('-m', '-models', dest="models", nargs="+", required=True)
    parser.add_argument('-s', '-structures', dest="structures", nargs="+", required=True)
    tmapp = parser.add_mutually_exclusive_group()
    tmapp.add_argument('-tm', '-tmscore', dest="tmscore", type=str)
    tmapp.add_argument('-ta', '-tmalign', dest="tmalign", type=str)
    args = parser.parse_args()

    if args.tmalign:
        tm = TMalign(args.tmalign)
    elif args.tmscore:
        tm = TMscore(args.tmscore)
    return tm.compare_structures(args.models, args.structures, all_vs_all=args.allvall)

if __name__ == "__main__":
    print main()
