
import os
import unittest

class Sequence(object):
    """A class to handle a fasta file"""
    
    def __init__( self ):
        """Initialise the object"""
        
        self.MAXWIDTH=80 # maximum width of any line 
        self._reset()
        return
    
    def parse_fasta(self,fastaFile):
        with open( fastaFile, "r") as f:
            self._parse_fasta( f )
        return

    def _reset( self ):
        """Reset the object"""
        self.headers = [] # title lines
        self.sequences = [] # The fasta sequences (just AA) as a single string
        return

    def _parse_fasta(self, fasta):
        """Parse the fasta file int our data structures & check for consistency 
        Args:
        fasta -- list of strings or open filehandle to read from the fasta file
        """

        self._reset() 
        sequence = None
        first=True
        for line in fasta:
            line = line.strip().rstrip(os.linesep)
            if first and not line.startswith( ">" ):
                raise RuntimeError,"FASTA files must start with a > character!"
            else:
                first=False
            if not line: continue # skip blank lines
            if line.startswith( ">" ):
                self.headers.append(line)
                if sequence: self.sequences.append(sequence)
                sequence=""
                continue
            sequence += line
        
        # add final sequence
        self.sequences.append(sequence)
        assert len(self.sequences)==len(self.headers)
        self.canonicalise()
        return
    
    def canonicalise(self):
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
    
    def length(self,seqNo=0):
        return len(self.sequences[seqNo])
    
    def numSequences(self):
        return len(self.sequences)
    
    def sequence(self):
        assert len(self.sequences) == 1,"> 1 sequence"
        return self.sequences[0]
    
    def pirStr(self,seqNo=0):
        """Return a canonical MAXWIDTH PIR representation of the file as a line-separated string"""
    
        pirStr = [ self.title + "\n\n" ] # Title plus the obligatory blank line
        
        # Now reformat to MAXWIDTH chars
        for chunk in range(0, len(self.sequances[seqNo]), self.MAXWIDTH):
            pirStr.append(self.sequences[seqNo][ chunk:chunk+self.MAXWIDTH ]+"\n")
            
        # Add last newline
        pirStr.append("\n")
        
        return pirStr
    
    def reformat_fasta( self, input_fasta, output_fasta):
        """
        Reformat the fasta file
        Needed because Rosetta has problems reading fastas. For it to be read, it has to have no spaces in the sequence,
        a name that is 4 characters and an underscore (ABCD_), everything has to be uppercase, and there has to be a
        return carriage at the end - this has to be linux formatted as when I make a fasta in windows, it doesnt recognize
        the return carriage.
        Rosetta has a lot of problems with fastas so we put in this script to deal with it.
        """
        
        self.parseFasta(input_fasta)
        with open( output_fasta, "w") as f:
            for line in self.fastaStr():
                f.write( line )
        return

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
    
    def writeFasta(self,fastaFile):
        if not len(self.sequences): raise RuntimeError,"No sequences have been read!"
        with open(fastaFile,'w') as f:
            for i, seq in enumerate(self.sequences):
                try:
                    h=self.headers[i]
                except IndexError:
                    h=">Sequence: {0}".format(i)
                f.write(h+'\n')
                for chunk in range(0, len(seq), self.MAXWIDTH):
                    f.write(seq[chunk:chunk+self.MAXWIDTH]+"\n")
            f.write("\n") # Add last newline
        return

class Test(unittest.TestCase):
    def testOK(self):
        """Reformat a fasta"""
 
        infasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTI
LAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD"""
        
        fp = Sequence()
        fp._parse_fasta( infasta.split( os.linesep ) )

        outfasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYW
ARYADWLFTTPLLLLDLALLVDADQGTILAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKA
ESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD

"""

        self.assertEqual( outfasta, "".join( fp.fastaStr() ) )
        self.assertEqual( fp.length(), 249)
              
    def testFailChar(self):
        
        infasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVXAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTI
LAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD"""

        fp = Sequence()
        self.assertRaises( RuntimeError, fp._parse_fasta, infasta.split( os.linesep ) )       

def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testOK'))
    suite.addTest(Test('testFailMulti'))
    suite.addTest(Test('testFailChar'))
    return suite
    
#
# Run unit tests
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(testSuite())
