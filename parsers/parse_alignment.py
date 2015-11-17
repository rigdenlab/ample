#!/usr/bin/env ccp4-python

import Bio.Alphabet
import Bio.Seq
import Bio.SeqIO
import unittest

class AlignmentParser(object):
    """ Parser for manipulation of MSAs """

    def remove_gaps(self, infile, outfile):
        """ Method to remove gaps in an alignment """
        # SeqIO instead of AlignIO to avoid error with diff seq lengths
        # Especially a problem when removing gaps
        records = list(Bio.SeqIO.parse(open(infile, 'r'), 'fasta'))
        
        records_formatted = self._remove_gaps(records)
        
        Bio.SeqIO.write(records_formatted, open(outfile, 'w'), "fasta")
        return

    def _remove_gaps(self, records):
        """ Format the string to remove the gaps
                and make every residue uppercase
        """
        for record in records:
            seq = self.format_seq(str(record.seq))
            record.seq = Bio.Seq.Seq(seq, Bio.Alphabet.single_letter_alphabet)
        return records

    def format_seq(self, seq):
        return str(seq).replace("-", "").upper()
##End AlignmentParser


class Test(unittest.TestCase):
    def setUp(self):
        self.ap = AlignmentParser()

    def testFormat(self):
        seq1 = "aaaaaa"
        seq2 = "A-A--A"
        seq3 = "a---a"

        ref_seq1 = "AAAAAA"
        ref_seq2 = "AAA"
        ref_seq3 = "AA"

        out_seq1 = self.ap.format_seq(seq1)
        out_seq2 = self.ap.format_seq(seq2)
        out_seq3 = self.ap.format_seq(seq3)
        
        self.assertEqual(ref_seq1, out_seq1)
        self.assertEqual(ref_seq2, out_seq2)
        self.assertEqual(ref_seq3, out_seq3)
##End Test

