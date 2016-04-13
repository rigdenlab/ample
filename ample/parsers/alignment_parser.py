
import os
import sys

BIOPYTHON_AVAILABLE=False
try:
    import Bio.AlignIO
    import Bio.Alphabet
    import Bio.pairwise2
    import Bio.Seq
    import Bio.SeqIO
    BIOPYTHON_AVAILABLE=True
except ImportError:
    msg = "Cannot import Biopython - some functionality will not be available." + os.linesep
    sys.stderr.write(msg)

class AlignmentParser(object):
    """ Parser for manipulation of MSAs """
    
    def a3mToTrimmed(self, alnFile, outFile):
        """Convert an .a3m HH-suite file to a FASTA format
        
        :alnFile: A3M format alignment file
        :outFile: MSA output file
        """
        # A3m is similar to FASTA except indels and lower case letters
        records = list(Bio.SeqIO.parse(open(alnFile, 'r'), 'fasta'))
        records_formatted = self._a3mToTrimmed(records)
        
        Bio.SeqIO.write(records_formatted, open(outFile, 'w'), "fasta")
        return
        
    def _a3mToTrimmed(self, records):
        """Takes two file handlers"""
        for record in records:
            record.seq = Bio.Seq.Seq("".join([c for c in str(record.seq) if not c.islower()]),
                                     Bio.Alphabet.single_letter_alphabet)
        return records
        
    def align_sequences(self, seq1, seq2):
        """Global pairwise alignment of two sequences
        
        :returns: aligned sequences as tuple
        """
        alignment = Bio.pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)
        return (alignment[-1][0], alignment[-1][1])
    
    def resetA3M(self, alnFile, outFile):
        """ Reset a A3M alignment to FASTA sequences
        
        Wrapper function for a3mToTrimmed() and removeGaps()
        
        :alnFile: A3M format alignment file
        :outFile: FASTA sequence file
        """
        records = list(Bio.SeqIO.parse(open(alnFile, 'r'), 'fasta'))
        
        records_trimmed = self._a3mToTrimmed(records)
        records_trimmed_gapped = self._removeGaps(records)
        
        Bio.SeqIO.write(records_trimmed_gapped, open(outFile, 'w'), "fasta")
        return
        
    def read(self, alnFile, alnFormat):
        """Read a multiple sequence alignment
        
        :alnFile: MSA file
        :alnFormat: MSA format
        
        :returns: Biopython MSA generator
        """
        return Bio.AlignIO.parse(open(alnFile, 'r'), alnFormat)
        
    def removeGaps(self, alnFile, outFile):
        """ Remove all gaps in a MSA
        
        :alnFile: FASTA alignment file
        :outFile: FASTA sequence file
        """
        # Instead of reading an AlignIO object read as SeqIO object to treat them
        # straight away as individual sequences
        
        records = list(Bio.SeqIO.parse(open(alnFile, 'r'), 'fasta'))
        records_formatted = self._removeGaps(records)
        
        Bio.SeqIO.write(records_formatted, open(outFile, 'w'), "fasta")
        return

    def _removeGaps(self, records):
        """ Format the string to remove the gaps
        """
        for record in records:  
            record.seq = Bio.Seq.Seq(str(record.seq).replace("-", ""), 
                                     Bio.Alphabet.single_letter_alphabet)
        return records

    def write(self, alignment, file, format):
        """Write a MSA object to a file"""
        Bio.AlignIO.write(alignment, open(file, 'w'), format)
        return
    

