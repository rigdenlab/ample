
import os
import re
import unittest


class FastaParser(object):
    """A class to handle a fasta file"""
    
    def __init__( self ):
        """Initialise the object"""
        
        self.length = None # The length of the parse fasta

    def _parse_fasta( self, fasta ):
        """Return the reformatted fasta as a list of strings 
        Args:
        fasta -- list of strings or open filehandle to read from the fasta file
        
        Return:
        Reformatted fasta as a list of strings
        """
        
        newfasta = []
        sequence = ""
        header=None
        for line in fasta:
            
            # remove whitespace and line-endings
            line = line.strip()
            line = line.rstrip( os.linesep )
            
            # Deal with header
            if line.startswith( ">" ):
                if header:
                    raise RuntimeError,"There appears to be more than one sequence in your fasta.\nPlease remove all but the first."
    
                header = line[0:80]+"\n"
                newfasta.append( header )
                continue
            
            sequence += line
        
        aa = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        
        newsequence = ""
        # Check for unwanted characters
        
        self.length=0
        for char in sequence:
            char = char.upper()
            if char not in aa:
                raise RuntimeError,"There appear to be non-standard AA in your sequence: '{0}'\nPlease format to only use standard AA.".format(char)
            newsequence+=char
            self.length+=1
        
        # Now reformat to 80 chars
        llen=80
        for chunk in range( 0, len(newsequence), llen ):
            newfasta.append( newsequence[ chunk:chunk+llen ]+"\n"  )
        
        # Add last newline
        newfasta.append("\n")
        
        return newfasta
        
    def reformat_fasta( self, input_fasta, output_fasta):
        """
        Reformat the fasta file
        Needed because Rosetta has problems reading fastas. For it to be read, it has to have no spaces in the sequence,
        a name that is 4 characters and an underscore (ABCD_), everything has to be uppercase, and there has to be a
        return carriage at the end - this has to be linux formatted as when I make a fasta in windows, it doesnt recognize
        the return carriage.
        Rosetta has a lot of problems with fastas so we put in this script to deal with it.
        """
        
        f = open( input_fasta, "r")
        newfasta =  self._parse_fasta( f )
        f.close()
        
        fasout=open( output_fasta, "w")
        for line in newfasta:
            fasout.write( line )
        fasout.close()
        
        return


class Test(unittest.TestCase):


    def testOK(self):
        """Reformat a fasta"""
 
        infasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTI
LAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD"""
        
        fp = FastaParser()
        newfasta = fp._parse_fasta( infasta.split( os.linesep ) )

        outfasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYW
ARYADWLFTTPLLLLDLALLVDADQGTILAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKA
ESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD

"""

        self.assertEqual( outfasta, "".join(newfasta) )
        self.assertEqual( fp.length, 249)
  
              
    def testFailMulti(self):
        
        infasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTI
LAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD
>3HAP:A|PDBID|CHAIN|SEQUENCE
"""

        fp = FastaParser()
        self.assertRaises( RuntimeError, fp._parse_fasta, infasta.split( os.linesep )  )
        
 
    def testFailChar(self):
        
        infasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVXAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTI
LAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD"""

        fp = FastaParser()
        self.assertRaises( RuntimeError, fp._parse_fasta, infasta.split( os.linesep ) )       
    
if __name__ == "__main__":
    unittest.main()
    #parse_fasta("/home/jmht/in.fasta", "/home/jmht/out.fasta")
