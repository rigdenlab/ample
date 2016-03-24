

from Bio import SeqIO
from Bio import SeqUtils

AA_ALPHABET = SeqUtils.IUPACData.protein_letters

def read(file, format):
    return SeqIO.parse(open(file, 'r'), format)

def canonicalise(seq):
    sequence = ""
    for char in seq:
        char = char.upper()
        if char not in AA_ALPHABET:
            msg = "There appear to be non-standard AA in your sequence: '{0}'\nPlease format to only use standard AA.".format(char)
            raise RuntimeError(msg)
        sequence += char  
    return sequence