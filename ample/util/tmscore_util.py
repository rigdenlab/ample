#!/usr/bin/env ccp4-python

import csv
import collections
import logging
import os
import warnings

from ample.parsers import alignment_parser
from ample.parsers import tmscore_parser
from ample.util import ample_util
from ample.util import pdb_edit

__author__ = "Felix Simkovic"
__date__ = "09.11.2015"
__version__ = "1.0"

LOGGER = logging.getLogger(__name__)

TMScoreModel = collections.namedtuple("TMScoreModel",
                                      ["name", "pdbin", "TMSCORE_log", "structure",
                                       "tm", "maxsub", "gdtts", "gdtha", "rmsd",
                                       "nr_residues_common"])


def tmscore_available():
    """
    Check if TMscore binary is available

    Returns
    -------
    bool
    """
    try:
        ample_util.find_exe("TMscore")
    except:
        return False
    return True


class TMscorer(object):
    """
    Wrapper to handle TMscoring for one or more structures

    Attributes
    ----------
    entries : list
       List containing the TMscore entries on a per-model basis
    pickle_file : str
       Path to a pickled data file containing AMPLE data
    structure : str
       Path to the reference structure
    tmscore_exe : str
       Path to the TMscore executable
    work_dir : str
       Path to the working directory
    """

    def __init__(self, structure, tmscore_exe, wdir=None):
        """
        Parameters
        ----------
        structure : str
           Path to the reference structure
        tmscore_exe : str
           Path to the TMscore executable
        work_dir : str
           Path to the working directory
        """
        self.entries = []
        self.pickle_file = None
        self.structure = structure
        self.tmscore_exe = tmscore_exe
        self.work_dir = wdir
        return

    def main(self, pdb_list_file, keep_modified_structures=False, identical_sequences=False):
        """
        Wrapper function for compare_to_structure

        Description
        -----------
        This function is primarily a wrapper around the actual comparison functions should be
        used in combination with a file containing the paths to protein models. Ultimately,
        this function makes use of the ``compare_to_structure()`` function.

        Parameters
        ----------
        pdb_list_file : str
           File containing a list of models
        keep_modified_structures : bool
           Flag to delete intermediate, modified structure files
        identical_sequences : bool
           Flag to avoid any modification of files due to sequence identity
        """
        # Read all pdb files from list
        pdbs = self._read_list(pdb_list_file)
        # Compare each pdb in list against structure defined
        self.entries = self.compare_to_structure(pdbs, 
                                                 keep_modified_structures=keep_modified_structures, 
                                                 identical_sequences=identical_sequences)
        return

    def compare_to_structure(self, pdb_list, keep_modified_structures=False, identical_sequences=False):
        """
        Compare a list of structures to a reference structure

        Parameters
        ----------
        pdb_list : list
           list containing paths to model files
        keep_modified_structures : bool
           Flag to delete intermediate, modified structure files
        identical_sequences : bool
           Flag to avoid any modification of files due to sequence identity
        
        Returns
        -------
        entries : list
           List of TMscore data entries on a per-model basis
        """
        if not identical_sequences and not alignment_parser.BIOPYTHON_AVAILABLE:
            msg = 'Cannot compare non-identical sequences as Biopython is not available'
            LOGGER.critical(msg)
            raise RuntimeError(msg)

        LOGGER.info('-------Evaluating decoys/models-------')        
        entries = []                                                       # For data storage

        structure = os.path.abspath(self.structure)                        # Path to reference structure
        structure_name = os.path.splitext(os.path.basename(structure))[0]  # Filename
        pdb_list_abs = [os.path.abspath(model) for model in pdb_list]      # Full paths to models

        # Create a TMscore logfile parser
        pt = tmscore_parser.TMscoreLogParser()

        for pdbin in pdb_list_abs:
            LOGGER.debug("Working on %s" % pdbin)
            pdbin_name = os.path.splitext(os.path.basename(pdbin))[0]  # Filename        
            if not os.path.exists(pdbin): 
                LOGGER.warning("Cannot find {0}".format(pdbin))
                continue
            
            # Modify structures to be identical as required by TMscore binary
            if not identical_sequences:
                pdbin_mod = os.path.join(self.work_dir, pdbin_name + "_mod.pdb")
                structure_mod = os.path.join(self.work_dir, pdbin_name + "_" + 
                                             structure_name + "_mod.pdb")
                self.mod_structures(pdbin, pdbin_mod, structure, structure_mod)    
            model = pdbin_mod if not identical_sequences else pdbin
            reference = structure_mod if not identical_sequences else structure
    
            log = os.path.join(self.work_dir, pdbin_name + "_tmscore.log")
            self.execute_comparison(model, reference, log)

            # Delete the modified structures if not wanted       
            if not keep_modified_structures and not identical_sequences:
                os.remove(pdbin_mod)
                os.remove(structure_mod)

            try:
                # Reset the TMscoreLogParser to default values of 0.0 for every score.
                pt.reset()
                # Parse the TMscore logfile to extract the scores
                pt.parse(log)
            except Exception as e:
                LOGGER.critical(e.msg)
                log = "None"

            entry = self._store(pdbin_name, pdbin, log, self.structure, pt)
            entries.append(entry)

        return entries

    def execute_comparison(self, model, reference, log=None):
        """
        Wrapper to execute the TMscore comparison command

        Paramters
        ---------
        model : str
           Path to the model structure file
        reference : str
           Path to the reference structure file
        log : str
           Path to the log file

        Returns
        -------
        return_code : int
           Return code of the process
        """
        # Create a command list and execute TMscore
        cmd = [self.tmscore_exe, model, reference]
        return_code = ample_util.run_command(cmd, logfile=log, directory=self.work_dir)
        return return_code
    
    def dump_csv(self, csv_file):
        """
        Dump the entry data to a csv file

        Parameters
        ----------
        csv_file : str
           Path to a file to write the data to

        Warnings
        --------
        This function was deprecated and will be removed in future releases
        """
        msg = "This function was deprecated and will be removed in a future release"
        warnings.warn(msg, DeprecationWarning, stacklevel=2)

        if not len(self.entries):
            return

        with open(csv_file, 'w') as f:
            fieldnames = self.entries[0]._asdict().keys()
            dw = csv.DictWriter(f, fieldnames=fieldnames)
            dw.writeheader()
            for e in self.entries:
                dw.writerow(e._asdict())
        LOGGER.info("Wrote csvfile: {0}".format(os.path.abspath(csv_file)))
        return

    def mod_structures(self, pdbin, pdbin_mod, structure, structure_mod):
        """
        Modify the two structure files to match each other

        Description
        -----------
        Structure files often contain unequal residue numberings, missing residues in the
        chain or other mal-formatted data. This function aims to remove such discrepancies
        to allow for the most accurate comparisons possible.

        Parameters
        ----------
        pdbin : str
           Path to the model pdb structure file
        pdbin_mod : str
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

        pdbin_seq = pdb_edit.sequence(pdbin).values()[0]
        structure_seq = pdb_edit.sequence(structure).values()[0]

        # Align the sequences to see how much of the predicted decoys are in the xtal
        aligned_seq_list = alignment_parser.AlignmentParser().align_sequences(pdbin_seq, 
                                                                              structure_seq)
        pdbin_seq_ali = aligned_seq_list[0]
        structure_seq_ali = aligned_seq_list[1]
        
        # Get the gaps in both sequences
        pdbin_gaps = self.find_gaps(pdbin_seq_ali)
        structure_gaps = self.find_gaps(structure_seq_ali)

        ## STAGE 1 - REMOVE RESIDUES ##
        pdbin_stage1 = ample_util.tmp_file_name(delete=False)
        structure_stage1 = ample_util.tmp_file_name(delete=False)

        # Get first residue number to adjust list of residues to remove
        pdbin_res1 = self.residue_one(pdbin)
        structure_res1 = self.residue_one(structure)
        
        # Match the residue lists to fit the residue 1 number
        pdbin_gaps = [i + structure_res1-1 for i in pdbin_gaps]
        structure_gaps = [i + pdbin_res1-1 for i in structure_gaps]

        # Use gaps of other sequence to even out
        pdb_edit.select_residues(pdbin, pdbin_stage1.name, delete=structure_gaps)
        pdb_edit.select_residues(structure, structure_stage1.name, delete=pdbin_gaps)

        ## STAGE 2 - RENUMBER RESIDUES ##
        pdb_edit.renumber_residues(pdbin_stage1.name, pdbin_mod)
        pdb_edit.renumber_residues(structure_stage1.name, structure_mod)

        os.unlink(pdbin_stage1.name)
        os.unlink(structure_stage1.name)

        logging.disable(logging.NOTSET)

        return

    def residue_one(self, pdb):
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
                index = int(line[5])
        return index

    def find_gaps(self, seq):
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
        return [i+1 for i, c in enumerate(seq) if c == "-"]
    
    def read_sequence(self, seq):
        """
        Determine the sequence offset

        Parameters
        ----------
        seq : str
           String of amino acids

        Returns
        -------
        offset : int
           Offset of sequence

        Warnings
        --------
        This function was deprecated and will be removed in future releases
        """

        msg = "This function was deprecated and will be removed in future release"
        warnings.warn(msg, DeprecationWarning, stacklevel=2)

        offset = 0
        for char in seq:
            if char == "-":
                offset += 1
            if char != "-":
                break
        return offset
    
    def _read_list(self, list_file):
        return [l.strip() for l in open(list_file, 'r')]
    
    def _store(self, name, pdbin, logfile, structure, pt):
        return TMScoreModel(name=name, pdbin=pdbin,
                            TMSCORE_log=logfile, structure=structure,
                            tm=pt.tm, maxsub=pt.maxsub, gdtts=pt.gdtts,
                            gdtha=pt.gdtha, rmsd=pt.rmsd,
                            nr_residues_common=pt.nr_residues_common)


