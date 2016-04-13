#!/usr/bin/env ccp4-python

import collections
import logging
import os
import tempfile

from ample.parsers import alignment_parser
from ample.parsers import tmscore_parser
from ample.util import ample_util
from ample.util import pdb_edit

__author__ = "Felix Simkovic"
__date__ = "09.11.2015"

LOGGER = logging.getLogger(__name__)

TMScoreModel = collections.namedtuple("TMScoreModel", 
                                      ["name", "pdbin", "TMSCORE_log", "structure", 
                                      "tm", "maxsub", "gdtts", "gdtha", "rmsd", 
                                      "nrResiduesCommon"])


def tmscoreAvail():
    """Check if TMscore binary is available"""
    try: ample_util.find_exe("TMscore")
    except: return False
    return True


class TMscorer(object):
    
    def __init__(self, structure, tmscore_exe, wdir=None):
        self.pickle_file = None
        self.structure = structure
        self.tmscore_exe = tmscore_exe
        self.work_dir = wdir
        return

    def main(self, pdb_list_file, keep_modified_structures=False, identical_sequences=False):
        """Wrapper function for compare_to_structure
        
        :pdb_list_file:            file containing a list of models
        :keep_modified_structures: do not delete any intermediate, modified structure files
        :sequence_identical:       avoid any modification of files due to sequence identity
        """
        # Read all pdb files from list
        pdbs = self._read_list(pdb_list_file)
        # Compare each pdb in list against structure defined
        self.entries = self.compare_to_structure(pdbs, 
                                                 keep_modified_structures=keep_modified_structures, 
                                                 identical_sequences=identical_sequences)
        return

    def compare_to_structure(self, pdb_list, keep_modified_structures=False, identical_sequences=False):
        """Compare a list of structures to a reference structure
        
        :pdb_list:                 list containing paths to model files
        :keep_modified_structures: do not delete any intermediate, modified structure files
        :sequence_identical:       avoid any modification of files due to sequence identity
        
        :returns:                  list of TMscore data per model
        """
        if not identical_sequences and not alignment_parser.BIOPYTHON_AVAILABLE:
            msg = 'Cannot compare non-identical sequences as Biopython is not available'
            LOGGER.critical(msg)
            raise RuntimeError(msg)

        LOGGER.info('-------Evaluating decoys/models-------')        
        entries = []                                                    # For data storage

        structure = os.path.abspath(self.structure)                        # Path to reference structure
        structure_name = os.path.splitext(os.path.basename(structure))[0]  # Filename
        pdb_list_abs = [ os.path.abspath(model) for model in pdb_list ]    # Full paths to models

        for pdbin in pdb_list_abs:
            LOGGER.debug("Working on %s" % pdbin)
            pdbin_name = os.path.splitext(os.path.basename(pdbin))[0]  # Filename        
            if not os.path.exists(pdbin): 
                LOGGER.warning("Cannot find {0}".format(pdbin))
                continue
            
            # Modify structures to be identical as required by TMscore binary
            if not identical_sequences:
                pdbin_mod     = os.path.join(self.work_dir, pdbin_name + "_mod.pdb")
                structure_mod = os.path.join(self.work_dir, pdbin_name + "_" + 
                                             structure_name + "_mod.pdb")
                self.mod_structures(pdbin, pdbin_mod, structure, structure_mod)    
            model     = pdbin_mod     if not identical_sequences else pdbin
            reference = structure_mod if not identical_sequences else structure
    
            log = os.path.join(self.work_dir, pdbin_name + "_tmscore.log")
            self.execute_comparison(model, reference, log)

            # Delete the modified structures if not wanted       
            if not keep_modified_structures and not identical_sequences:
                os.remove(pdbin_mod)
                os.remove(structure_mod)
            
            # Do the try clause here to allow anything that is required from here to throw
            # exceptions. In that case we revert to the TMscoreLogParser default values of
            # 0.0 for every score.
            pt = tmscore_parser.TMscoreLogParser()
            
            try:
                pt.parse(log)
            except Exception as e:
                LOGGER.critical(e.msg)
                log = "None"
                
            entry = self._store(pdbin_name, pdbin, log, self.structure, pt)
            entries.append(entry)
                
        return entries

    def execute_comparison(self, model, reference, log=None):
        # Create a command list and execute TMscore
        cmd = [ self.tmscore_exe, model, reference ]
        p = ample_util.run_command(cmd, logfile=log, directory=self.work_dir)
        return p
    
    def dump_csv(self, csv_file):
        if not len(self.entries): return
        import csv
        with open(csv_file, 'w') as f:
            fieldnames = self.entries[0]._asdict().keys()
            dw = csv.DictWriter(f, fieldnames=fieldnames)
            dw.writeheader()
            for e in self.entries:
                dw.writerow(e._asdict())
        print "Wrote csvfile: {0}".format(os.path.abspath(csv_file))
        return

    def mod_structures(self, pdbin, pdbin_mod, structure, structure_mod):
        """Make sure the decoy and the xtal pdb align to get an accurate TM-score"""
        
        # Disable the info logger to not spam the user with which chain of native extracted.
        # Happens for every model + native below
        # http://stackoverflow.com/questions/2266646/how-to-i-disable-and-re-enable-console-logging-in-python
        logging.disable( logging.INFO )

        pdbin_seq     = pdb_edit.sequence(pdbin).values()[0]
        structure_seq = pdb_edit.sequence(structure).values()[0]

        # Align the sequences to see how much of the predicted decoys are in the xtal
        aligned_seq_list = alignment_parser.AlignmentParser().align_sequences(pdbin_seq, 
                                                                              structure_seq)
        pdbin_seq_ali     = aligned_seq_list[0]
        structure_seq_ali = aligned_seq_list[1]
        
        # Get the gaps in both sequences
        pdbin_gaps     = self.find_gaps(pdbin_seq_ali)
        structure_gaps = self.find_gaps(structure_seq_ali)

        ## STAGE 1 - REMOVE RESIDUES ##
        pdbin_stage1 = tempfile.NamedTemporaryFile(delete=False)
        pdbin_stage1.close()
        structure_stage1 = tempfile.NamedTemporaryFile(delete=False)
        structure_stage1.close()

        # Get first residue number to adjust list of residues to remove
        pdbin_res1     = self.residue_one(pdbin)
        structure_res1 = self.residue_one(structure)
        
        # Match the residue lists to fit the residue 1 number
        pdbin_gaps     = [ i + structure_res1-1 for i in pdbin_gaps ]
        structure_gaps = [ i + pdbin_res1-1     for i in structure_gaps ]

        # Use gaps of other sequence to even out
        pdb_edit.select_residues(pdbin, pdbin_stage1.name, delete=structure_gaps)
        pdb_edit.select_residues(structure, structure_stage1.name, delete=pdbin_gaps)

        ## STAGE 2 - RENUMBER RESIDUES ##
        pdb_edit.renumber_residues(pdbin_stage1.name, pdbin_mod)
        pdb_edit.renumber_residues(structure_stage1.name, structure_mod)

        os.unlink(pdbin_stage1.name)
        os.unlink(structure_stage1.name)

        return

    def residue_one(self, pdb):
        for line in open(pdb, 'r'):
            if line.startswith("ATOM"):
                line = line.split()
                return int(line[5])

    def find_gaps(self, seq):
        return [ i+1 for i, c in enumerate(seq) if c=="-" ]
    
    def read_sequence(self, seq):
        offset = 0
        for char in seq:
            if char=="-": offset += 1
            if char!="-": break
        return offset
    
    def _read_list(self, list_file):
        return [ l.strip() for l in open(list_file, 'r') ]
    
    def _store(self, name, pdbin, logfile, structure, pt):
        return TMScoreModel(name=name, pdbin=pdbin, 
                           TMSCORE_log=logfile, structure=structure, 
                           tm=pt.tm, maxsub=pt.maxsub, gdtts=pt.gdtts,
                           gdtha=pt.gdtha, rmsd=pt.rmsd,
                           nrResiduesCommon=pt.nrResiduesCommon)


