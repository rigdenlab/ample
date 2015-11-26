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

# 3rd party
import numpy

if not "CCP4" in os.environ.keys(): raise RuntimeError('CCP4 not found')
root = os.path.join(os.environ['CCP4'], "share", "ample")
#root = os.sep.join( os.path.abspath(__file__).split( os.sep )[:-2] )
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
    def __init__(self, structure, tmscore_exe, wdir=None):
        self.logger = logging.getLogger()
        self.pickle_file = None
        self.structure = structure
        self.tmscore_exe = tmscore_exe
        self.workingDIR = wdir
        return

    def main(self, pdb_list_file):
        # Read all pdb files from list
        pdbs = self._read_list(pdb_list_file)

        # Compare each pdb in list against structure defined
        self.entries = self.compare_to_structure(pdbs)

        # Dumpe all data in a pickle file
        self.pickle_file = os.path.join(self.workingDIR, "tmresults.pkl")
        pickle.dump(self.entries, open(self.pickle_file, 'w'))
        return

    def compare_to_structure(self, pdb_list):
        self.logger.info('-------Evaluating decoys/models-------')

        entries = []

        for pdb in pdb_list:
            if not os.path.exists(pdb): continue

            name = os.path.basename(pdb).rsplit(".", 1)[0]
            pdbin = os.path.abspath(pdb)

            self.logger.info("Working on %s" % name)

            # Write a temporary file with the aligned residues from the native structure
            structure_mod = tempfile.NamedTemporaryFile(delete=False)
            structure_mod.close()
            pdbin_mod = tempfile.NamedTemporaryFile(delete=False)
            pdbin_mod.close()
            
            # Initialise the logparser that stores all the scores
            pt = parse_tmscore.TMscoreLogParser()
            
            # Do the try clause here to allow anything that is required from here to throw
            # exceptions. In that case we revert to the TMscoreLogParser default values of
            # 0.0 for every score. 
            try:
                self.mod_structures(pdbin, pdbin_mod.name, self.structure, structure_mod.name)
            
                # Create a command list and execute TMscore
                log = os.path.join(self.workingDIR, name+".tmscore.log")
                cmd = [ self.tmscore_exe, pdbin_mod.name, structure_mod.name ]
                p = ample_util.run_command(cmd, logfile=log, directory=self.workingDIR)
    
                os.unlink(structure_mod.name)
                os.unlink(pdbin_mod.name)
                
                # Parse the log file to the parser
                pt.parse(log)
                
            except Exception:
                log = "None"
                
            entry = self._store(name, pdbin, log, self.structure, pt)
            entries.append(entry)
                
        return entries

    def mod_structures(self, pdbin, pdbin_mod, structure, structure_mod):
        ''' Make sure the decoy and the xtal pdb align to get an accurate Tm score '''
        
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
        pdbin_gaps     = [i + structure_res1-1 for i in pdbin_gaps]
        structure_gaps = [i + pdbin_res1-1     for i in structure_gaps]

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
        gaps = [i+1 for i, c in enumerate(seq) if c=="-"]
        return gaps
    
    def read_sequence(self, seq):
        offset = 0
        for char in seq:
            if char=="-": offset += 1
            if char!="-": break
        return offset
    
    def _read_list(self, list_file):
        return [l.strip() for l in open(list_file, 'r')]
    
    def _store(self, name, pdbin, logfile, structure, pt):
        entry = {"name": name,
                 "pdbin": pdbin,
                 "TMSCORE_log": logfile,
                 "structure": structure,
                 "tm": pt.tm,
                 "maxsub": pt.maxsub,
                 "gdtts": pt.gdtts,
                 "gdtha": pt.gdtha,
                 "rmsd": pt.rmsd,
                 "nrResiduesCommon": pt.nrResiduesCommon}
        return entry
    
    
class Statistics(object):
    def __init__(self, pickle_file, pdb_list_file):
        self.pdb_list_file = pdb_list_file
        self.pickle_file = pickle_file
        self.mean = 0.0
        self.median = 0.0
        return

    def main(self, score):
        # Load the required data
        pdb_names = self._read_list(self.pdb_list_file)
        data = pickle.load(open(self.pickle_file, 'r'))
        
        # Obtain a list of all defined scores for defined pdbs
        score_list = self.extract(data, pdb_names, score)
        score_list_sorted = sorted(score_list)
        
        assert len(score_list_sorted) > 0, \
                "No scores in score list file"
        
        assert len(score_list_sorted)==len(pdb_names), \
                "Divergent counts between scores and pdb names"

        self.mean   = numpy.mean(score_list_sorted)
        self.median = numpy.median(score_list_sorted)
        return

    def extract(self, data, pdb_names, score):
        return [entry[score] for entry in data \
                        if entry['name'] in pdb_names]

    def _read_list(self, list_file):
        return [os.path.basename(l.strip()).rsplit('.', 1)[0] \
                    for l in open(list_file, 'r')]
        

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


class TestStatistics(unittest.TestCase):
    def testExtract(self):
        data = [{"name": "test1", "pdbin": "test1",
                 "TMSCORE_log": "test1", "structure": "test1",
                 "tm": 0.523, "maxsub": 0.333, "gdtts": 0.355,
                 "gdtha": 0.424, "rmsd": 1.345, "nrResiduesCommon": 10},
                {"name": "test2", "pdbin": "test2",
                 "TMSCORE_log": "test2", "structure": "test2",
                 "tm": 0.140, "maxsub": 0.222, "gdtts": 0.234,
                 "gdtha": 0.100, "rmsd": 10.444, "nrResiduesCommon": 10},
                {"name": "test3", "pdbin": "test3",
                 "TMSCORE_log": "test3", "structure": "test3",
                 "tm": 0.3, "maxsub": 0.5, "gdtts": 0.1,
                 "gdtha": 0.6, "rmsd": 1.0, "nrResiduesCommon": 10},
                {"name": "test4", "pdbin": "test4",
                 "TMSCORE_log": "test4", "structure": "test4",
                 "tm": 0.6, "maxsub": 0.8, "gdtts": 0.9,
                 "gdtha": 0.4, "rmsd": 9.345, "nrResiduesCommon": 10}]

        s = Statistics("foo", "bar")


        pdb_list = ["test1", "test3", "test4"]

        ref_data_tm = [0.523, 0.3, 0.6]
        out_data_tm = s.extract(data, pdb_list, "tm")
        self.assertItemsEqual(ref_data_tm, out_data_tm)

        ref_data_rmsd = [1.345, 1.0, 9.345]
        out_data_rmsd = s.extract(data, pdb_list, "rmsd")
        self.assertItemsEqual(ref_data_rmsd, out_data_rmsd)

        ref_data_gdtts = [0.355, 0.1, 0.9]
        out_data_gdtts = s.extract(data, pdb_list, "gdtts")
        self.assertItemsEqual(ref_data_gdtts, out_data_gdtts)

        ref_data_gdtha = [0.424, 0.6, 0.4]
        out_data_gdtha = s.extract(data, pdb_list, "gdtha")
        self.assertItemsEqual(ref_data_gdtha, out_data_gdtha)

        ref_data_maxsub = [0.333, 0.5, 0.8]
        out_data_maxsub = s.extract(data, pdb_list, "maxsub")
        self.assertItemsEqual(ref_data_maxsub, out_data_maxsub)


        pdb_list_2 = ["test2", "test4"]
        ref_data_maxsub_2 = [0.222, 0.8]
        out_data_maxsub_2 = s.extract(data, pdb_list_2, "maxsub")
        self.assertItemsEqual(ref_data_maxsub_2, out_data_maxsub_2)


if __name__ == "__main__":
    # Check whether we can import Biopythin
    if not _BIOPYTHON: sys.exit("Upgrade to a CCP4 version with the new interface to use this script")

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-score', metavar='[ tm | maxsub | gdtts | gdtha ]',
                        type=str, default='tm', help='The score to extract [default: tm]')
    parser.add_argument("structure")
    parser.add_argument("pdb_list_file")
    parser.add_argument("tmscore")
    parser.add_argument("-wdir", type=str, default=os.getcwd())
    args = parser.parse_args()

    assert args.score in ["tm", "maxsub", "gdtts", "gdtha", "rmsd"], "Unknown score %s" % args.score

    t = TMscorer(os.path.abspath(args.structure), os.path.abspath(args.tmscore),
                os.path.abspath(args.wdir))
    t.main(args.pdb_list_file)

    pkl_file = t.pickle_file

    s = Statistics(pkl_file, os.path.abspath(args.pdb_list_file))
    s.main(args.score)
    print "Mean", s.mean
    print "Median", s.median
