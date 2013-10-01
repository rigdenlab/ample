
import os
import re
import unittest


class FastaParser(object):
    """A class to handle a fasta file"""
    
    def __init__( self ):
        """Initialise the object"""
        
        self.MAXWIDTH=80 # maximum width of any line 
        self._reset()

    def _reset( self ):
        """Reset the object"""
        
        self.title = None # title line
        self.length = None # The length of the parse fasta
        self.seqStr = None # The fasta sequence (just AA) as a single string
        self.formattedFasta = None # The reformatted fasta

    def _parse_fasta( self, fasta ):
        """Parse the fasta file int our data structures & check for consistency 
        Args:
        fasta -- list of strings or open filehandle to read from the fasta file
        """

        self._reset() 

        sequence = ""
        for line in fasta:
            
            # remove whitespace and line-endings
            line = line.strip()
            line = line.rstrip( os.linesep )
            
            # Deal with header
            if line.startswith( ">" ):
                if self.title:
                    raise RuntimeError,"There appears to be more than one sequence in your fasta.\nPlease remove all but the first."
                
                if len(line) > self.MAXWIDTH:
                    self.title = line[0:self.MAXWIDTH]
                else:
                    self.title = line
                continue
            
            sequence += line
        
        aa = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        self.seqStr = ""
        
        # Check for unwanted characters and get canonical string and length
        self.length=0
        for char in sequence:
            char = char.upper()
            if char not in aa:
                raise RuntimeError,"There appear to be non-standard AA in your sequence: '{0}'\nPlease format to only use standard AA.".format(char)
            self.seqStr +=char
            self.length +=1
        
        return
    
    def fastaStr(self):
        """Return a canonical MAXWIDTH fasta representation of the file as a line-separated string"""
    
        formattedFasta = [ self.title + "\n" ]
        # Now reformat to MAXWIDTH chars
        for chunk in range( 0, len(self.seqStr), self.MAXWIDTH ):
            formattedFasta.append( self.seqStr[ chunk:chunk+self.MAXWIDTH ]+"\n"  )
            
        # Add last newline
        formattedFasta.append("\n")
        
        return formattedFasta
    
    def pirStr(self):
        """Return a canonical MAXWIDTH PIR representation of the file as a line-separated string"""
    
        pirStr = [ self.title + "\n\n" ] # Title plus the obligatory blank line
        
        # Now reformat to MAXWIDTH chars
        for chunk in range( 0, len(self.seqStr), self.MAXWIDTH ):
            pirStr.append( self.seqStr[ chunk:chunk+self.MAXWIDTH ]+"\n"  )
            
        # Add last newline
        pirStr.append("\n")
        
        return pirStr
    
#     def getSeqStr(self, fastaFile ):
#         """Return the fasta string from the given file"""
#         
#         f = open( fastaFile, "r")
#         self._parse_fasta( f )
#         f.close()
#         
#         return self.seqStr
        
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
        self._parse_fasta( f )
        f.close()
        
        fasout=open( output_fasta, "w")
        for line in self.fastaStr():
            fasout.write( line )
        fasout.close()
        
        return

    def toPir(self, input_fasta, output_pir=None ):
        """Take a fasta file and output the corresponding PIR file"""

        f = open( input_fasta, "r")
        self._parse_fasta( f )
        f.close()
        
        if not output_pir:
            dir, ffile = os.path.split( input_fasta )
            fname = os.path.splitext(ffile)[0]
            output_pir = os.path.join( dir, fname+".pir" )
        
        pirout=open( output_pir, "w")
        for line in self.pirStr():
            pirout.write( line )
        pirout.close()
        
        return   

class Test(unittest.TestCase):


    def testOK(self):
        """Reformat a fasta"""
 
        infasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTI
LAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD"""
        
        fp = FastaParser()
        fp._parse_fasta( infasta.split( os.linesep ) )

        outfasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYW
ARYADWLFTTPLLLLDLALLVDADQGTILAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKA
ESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD

"""

        self.assertEqual( outfasta, "".join( fp.fastaStr() ) )
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
