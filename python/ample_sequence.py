
import os
import unittest

import pdb_edit

class Sequence(object):
    """A class to handle a fasta file"""
    
    def __init__(self,fasta=None,pdb=None,canonicalise=False):
        """Initialise the object"""
        
        self.MAXWIDTH=80 # maximum width of any line 
        self.name='unknown'
        self.headers = [] # title lines
        self.sequences = [] # The fasta sequences (just AA) as a single string
        self.pdbs = [] # pdb files that sequences were read from
        self.fasta_files = [] # The fasta files the sequence were ready from
        
        if fasta:
            self.from_fasta(fasta,canonicalise=canonicalise)
        elif pdb:
            self.from_pdb(pdbin=pdb)
        return
    
    def from_fasta(self,fasta_file,canonicalise=True):
        name=os.path.splitext(os.path.basename(fasta_file))[0]
        self.name=name
        with open(fasta_file, "r") as f:
            self._parse_fasta(f,fasta_file=fasta_file,canonicalise=canonicalise)
        return
    
    def from_pdb(self,pdbin):
        name=os.path.splitext(os.path.basename(pdbin))[0]
        self.name=name
        chain2seq = pdb_edit.sequence(pdbin)
        assert len(chain2seq),"Could not read sequence from pdb: {0}".format(pdbin)
        self.headers = []
        self.sequences = []
        self.pdbs = []
        self.fasta_files = []
        for chain in sorted(chain2seq.keys()):
            seq=chain2seq[chain]
            self.headers.append(">From pdb: {0} chain {1} length {2}".format(name,chain,len(seq)))
            self.sequences.append(seq)
            self.pdbs.append(pdbin)
            self.fasta_files.append(None) # Need to make sure pdb abd fasta arrays have same length
        return

    def _parse_fasta(self, fasta, fasta_file=None, canonicalise=True):
        """Parse the fasta file int our data structures & check for consistency 
        Args:
        fasta -- list of strings or open filehandle to read from the fasta file
        """
        self.headers = []
        self.sequences = []
        self.fasta_files = [] 
        self.pdbs = [] 
        sequence = None
        first=True
        for line in fasta:
            line = line.strip().rstrip(os.linesep)
            if first and not line.startswith(">"):
                raise RuntimeError,"FASTA files must start with a > character!"
            else:
                first=False
            if not line: continue # skip blank lines
            if line.startswith(">"):
                self.headers.append(line)
                if sequence:
                    self.sequences.append(sequence)
                    if fasta_file:
                        self.fasta_files.append(fasta_file)
                        self.pdbs.append(None) # Need to make sure pdb abd fasta arrays have same length
                sequence=""
                continue
            sequence += line
        
        # add final sequence
        self.sequences.append(sequence)
        assert len(self.sequences)==len(self.headers)
        if canonicalise: self.canonicalise()
        return
    
    def canonicalise(self):
        """
        Reformat the fasta file
        Needed because Rosetta has problems reading fastas. For it to be read, it has to have no spaces in the sequence,
        a name that is 4 characters and an underscore (ABCD_), everything has to be uppercase, and there has to be a
        return carriage at the end - this has to be linux formatted as when I make a fasta in windows, it doesnt recognize
        the return carriage.
        Rosetta has a lot of problems with fastas so we put in this script to deal with it.
        """
        aa = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        for i,seq in enumerate(self.sequences):
            cs=""
            for char in seq:
                char = char.upper()
                if char not in aa:
                    raise RuntimeError,"There appear to be non-standard AA in your sequence: '{0}'\nPlease format to only use standard AA.".format(char)
                cs+=char
            self.sequences[i]=cs
        return
    
    def cat(self,seq2,fasta_file):
        """Concatentate two (single chain) sequences into a single fasta file"""
        headers = [self.headers[0], seq2.headers[0]]
        sequences = [self.sequences[0], seq2.sequences[0]]
        s = self._fasta_str(headers,sequences)
        with open(fasta_file,'w') as f:
            for line in s: f.write(line)
        return 
    
    def fasta_str(self):
        if not len(self.sequences): raise RuntimeError,"No sequences have been read!"
        headers=[]
        for i, seq in enumerate(self.sequences):
            try:
                h=self.headers[i]
            except IndexError:
                h=">Sequence: {0} Length: {1}".format(i,len(seq))
            headers.append(h)
        return self._fasta_str(headers, self.sequences)
    
    def _fasta_str(self,headers,sequences):
        s=""
        for i, seq in enumerate(sequences):
            s += headers[i]+'\n'
            for chunk in range(0, len(seq), self.MAXWIDTH):
                s += seq[chunk:chunk+self.MAXWIDTH]+"\n"
            s += "\n"
        s+="\n" # Add last newline
        return s

    def length(self,seq_no=0):
        return len(self.sequences[seq_no])
    
    def numSequences(self):
        return len(self.sequences)
    
    def sequence(self,seq_no=0):
        return self.sequences[seq_no]
    
    def pirStr(self,seqNo=0):
        """Return a canonical MAXWIDTH PIR representation of the file as a line-separated string"""
    
        pirStr = [ self.title + "\n\n" ] # Title plus the obligatory blank line
        
        # Now reformat to MAXWIDTH chars
        for chunk in range(0, len(self.sequances[seqNo]), self.MAXWIDTH):
            pirStr.append(self.sequences[seqNo][ chunk:chunk+self.MAXWIDTH ]+"\n")
            
        # Add last newline
        pirStr.append("\n")
        
        return pirStr

    def toPir(self, input_fasta, output_pir=None ):
        """Take a fasta file and output the corresponding PIR file"""
        self.parseFasta(input_fasta)
        if not output_pir:
            dir, ffile = os.path.split(input_fasta)
            fname = os.path.splitext(ffile)[0]
            output_pir = os.path.join( dir, fname+".pir" )
        with open( output_pir, "w") as pirout:
            for line in self.pirStr():
                pirout.write(line)
        return

    def write_fasta(self,fasta_file):
        with open(fasta_file,'w') as f:
            for s in self.fasta_str():
                f.write(s)
        return
    
    def __add__(self,other):
        self.headers += other.headers
        self.sequences += other.sequences
        self.pdbs += other.pdbs
        self.fasta_files += other.fasta_files
        return self
    
class Test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = cls.thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -1 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        return

    def testOK(self):
        """Reformat a fasta"""
 
        infasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTI
LAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD"""
        
        fp = Sequence()
        fp._parse_fasta(infasta.split(os.linesep))

        outfasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYW
ARYADWLFTTPLLLLDLALLVDADQGTILAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKA
ESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD

"""

        self.assertEqual(outfasta, "".join(fp.fasta_str()))
        self.assertEqual(fp.length(), 249)
        return
              
    def testFailChar(self):
        
        infasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVXAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTI
LAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD"""

        fp = Sequence()
        self.assertRaises( RuntimeError, fp._parse_fasta, infasta.split(os.linesep))    
        return
    
    def testAdd(self):
        s1 = Sequence(pdb=os.path.join(self.testfiles_dir,'1GU8.pdb'))
        s2 = Sequence(fasta=os.path.join(self.testfiles_dir,'2uui.fasta'))
        s1 += s2
        
        self.assertTrue(len(s1.sequences),2)
        self.assertTrue(len(s1.headers),2)
        self.assertTrue(len(s1.pdbs),2)
        self.assertTrue(len(s1.fasta_files),2)
        
        # Test write for the hell of it
        out_fasta = 'seq_add.fasta'
        s1.write_fasta(out_fasta)
        os.unlink(out_fasta)
        
        return
    
#
# Run unit tests
if __name__ == "__main__":
    unittest.main(verbosity=2)
