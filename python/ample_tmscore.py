#!/usr/bin/env ccp4-python

'''
09.11.2015

@author: hlfsimko
'''

# System
import argparse
import collections
import logging
import os
import pickle
import sys
import tempfile
import unittest

if "CCP4_AMPLE_ROOT" in os.environ.keys() and "CCP4" in os.environ.keys():
    root = os.environ["CCP4_AMPLE_ROOT"]
elif "CCP4" in os.environ.keys():
    root = os.path.join(os.environ["CCP4"], "share", "ample")
else:
    raise RuntimeError('CCP4 not found')

sys.path.insert(0, os.path.join(root, "parsers"))

# Custom
import ample_util
import parse_tmscore
import pdb_edit

try:
    import parse_alignment
    _BIOPYTHON = True
except ImportError:
    _BIOPYTHON = False


class TMscorer(object):
    
    _TMSCORE_DATA_FACTORY = collections.namedtuple("Model", 
                                                   ["name", "pdbin", "TMSCORE_log", "structure", 
                                                    "tm", "maxsub", "gdtts", "gdtha", "rmsd", 
                                                    "nrResiduesCommon"])
    
    def __init__(self, structure, tmscore_exe, wdir=None):
        self.logger = logging.getLogger()
        self.pickle_file = None
        self.structure = structure
        self.tmscore_exe = tmscore_exe
        self.workingDIR = wdir
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
        
        entries = []                                                    # For data storage
        
        self.logger.info('-------Evaluating decoys/models-------')

        structure = os.path.abspath(self.structure)                     # Path to reference structure
        structure_name = os.path.basename(structure).rsplit(".", 1)[0]  # Filename
        pdb_list_abs = [ os.path.abspath(model) for model in pdb_list ] # Full paths to models

        for pdbin in pdb_list_abs:
            self.logger.debug("Working on %s" % pdbin)
            pdbin_name     = os.path.basename(pdbin).rsplit(".", 1)[0]  # Filename        
            if not os.path.exists(pdbin): 
                self.logger.warning("Cannot find {0}".format(pdbin))
                continue
            
            # Modify structures to be identical as required by TMscore binary
            if not identical_sequences:
                pdbin_mod     = os.path.join(self.workingDIR, pdbin_name + "_mod.pdb")
                structure_mod = os.path.join(self.workingDIR, pdbin_name + "_" + structure_name + "_mod.pdb")
                self.mod_structures(pdbin, pdbin_mod, structure, structure_mod)    
            model     = pdbin_mod     if not identical_sequences else pdbin
            reference = structure_mod if not identical_sequences else structure
    
            # Create a command list and execute TMscore
            log = os.path.join(self.workingDIR, pdbin_name + "_tmscore.log")
            cmd = [ self.tmscore_exe, model, reference ]
            p = ample_util.run_command(cmd, logfile=log, directory=self.workingDIR)

            # Delete the modified structures if not wanted       
            if not keep_modified_structures and not identical_sequences:
                os.remove(pdbin_mod)
                os.remove(structure_mod)
            
            # Do the try clause here to allow anything that is required from here to throw
            # exceptions. In that case we revert to the TMscoreLogParser default values of
            # 0.0 for every score.
            pt = parse_tmscore.TMscoreLogParser()
            
            try:
                pt.parse(log)
            except Exception as e:
                self.logger.critical(e.msg)
                log = "None"
                
            entry = self._store(pdbin_name, pdbin, log, self.structure, pt)
            entries.append(entry)
                
        return entries

    def mod_structures(self, pdbin, pdbin_mod, structure, structure_mod):
        """Make sure the decoy and the xtal pdb align to get an accurate TM-score"""
        
        if not _BIOPYTHON:
            raise ImportError
        
        # Disable the info logger to not spam the user with which chain of native extracted.
        # Happens for every model + native below
        # http://stackoverflow.com/questions/2266646/how-to-i-disable-and-re-enable-console-logging-in-python
        logging.disable( logging.INFO )

        pdbin_seq     = pdb_edit.sequence(pdbin).values()[0]
        structure_seq = pdb_edit.sequence(structure).values()[0]

        # Align the sequences to see how much of the predicted decoys are in the xtal
        aligned_seq_list = parse_alignment.AlignmentParser().align_sequences(pdbin_seq, structure_seq)
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
        return self._TMSCORE_DATA_FACTORY(name=name, pdbin=pdbin, 
                                          TMSCORE_log=logfile, structure=structure, 
                                          tm=pt.tm, maxsub=pt.maxsub, gdtts=pt.gdtts,
                                          gdtha=pt.gdtha, rmsd=pt.rmsd,
                                          nrResiduesCommon=pt.nrResiduesCommon)


class TestTMScore(unittest.TestCase):
    def setUp(self):
        self.TM = TMscorer("foo", "bar", "cho")

    def testRead(self):
        seq1 = "-----AAAA---"
        ref1_offset = 5

        seq2 = "AAA---"
        ref2_offset = 0

        seq3 = "AAAA"
        ref3_offset = 0

        offset1 = self.TM.read_sequence(seq1)
        offset2 = self.TM.read_sequence(seq2)
        offset3 = self.TM.read_sequence(seq3)

        self.assertEqual(ref1_offset, offset1)
        self.assertEqual(ref3_offset, offset2)
        self.assertEqual(ref3_offset, offset3)

    def testGaps(self):
        seq1 = "AAAA---AA--AA"
        ref_gaps1 = [5, 6, 7, 10, 11]

        seq2 = "---AA-AA"
        ref_gaps2 = [1, 2, 3, 6]

        seq3 = "-AAA--"
        ref_gaps3 = [1, 5, 6]

        gaps1 = self.TM.find_gaps(seq1)
        gaps2 = self.TM.find_gaps(seq2)
        gaps3 = self.TM.find_gaps(seq3)

        self.assertEqual(ref_gaps1, gaps1)
        self.assertEqual(ref_gaps2, gaps2)
        self.assertEqual(ref_gaps3, gaps3)


if __name__ == "__main__":
    # Check whether we can import BioPython
    if not _BIOPYTHON: sys.exit("Upgrade to CCP4 version 7.0 or greater to use this script")
    import ample_statistics     # Only need this here

    parser = argparse.ArgumentParser()
    parser.add_argument("--identical", dest="identical", action="store_true", default=False,
                        help="reference structure and models have identical sequences (default: False)")
    parser.add_argument("--keep", dest="keep", action="store_true", default=False,
                        help="keep intermediate structures (default: False)")
    parser.add_argument("structure",
                        help="reference structure")
    parser.add_argument("pdb_list_file",
                        help="list containing model structures")
    parser.add_argument("tmscore",
                        help="TMscore binary")
    parser.add_argument("-d", dest="wdir", type=str, default=os.getcwd(),
                        help="working directory")
    args = parser.parse_args()

    t = TMscorer(os.path.abspath(args.structure), os.path.abspath(args.tmscore),
                os.path.abspath(args.wdir))
    t.main(args.pdb_list_file, keep_modified_structures=args.keep, identical_sequences=args.identical)
    
    tmscores = [ i.tm for i in t.entries ] 
    print "Median TM-score: {0}".format(round(ample_statistics.median(tmscores), 3))
    print "Mean   TM-score: {0}".format(round(ample_statistics.mean(tmscores), 3))
