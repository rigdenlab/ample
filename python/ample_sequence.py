#!/usr/bin/env ccp4-python
import os
import unittest

import pdb_edit

class Sequence(object):
    """A class to handle a fasta file"""
    
    def __init__(self, fasta=None, pdb=None, canonicalise=False):
        """Initialise the object"""
        
        self.MAXWIDTH=80 # maximum width of any line 
        self.headers = [] # title lines
        self.sequences = [] # The fasta sequences (just AA) as a single string
        self.resseqs = [] # The fasta sequences (just AA) as a single string
        self.pdbs = [] # pdb files that sequences were derived from
        self.chains = [] # The chain that the sequence belongs too (if applicable)
        self.fasta_files = [] # The fasta files the sequence were read from
        
        if fasta:
            self.from_fasta(fasta, canonicalise=canonicalise)
        elif pdb:
            self.from_pdb(pdbin=pdb)
        return
    
    def add_pdb_data(self, pdbin):
        """Add the resseq information from a pdb to ourselves when we already have the sequence information from a fasta file - such as from an alignment
        We assume that there will be gaps (-) in the sequence and the letters may be upper or lower case
        Currently this only supports adding data for single-chain pdbs
        """
        assert len(self.headers) and len(self.sequences)
        
        assert os.path.isfile(pdbin),"Cannot find pdb file: {0}".format(pdbin)
        
        # Assume the name of the pdb is the first part of the header
        fname = os.path.basename(pdbin)
        name, ext = os.path.splitext(fname)
        assert ext == '.pdb', "PDB files must have extension .pdb"
        
        # Find where in the list of data this sequence is
        got=False
        for idx, h in enumerate(self.headers):
            n = h[1:].split('.')[0]
            if n == name:
                got=True
                break
        
        if not got: raise RuntimeError,"Could not find matching pdb name for {0} in headers {1}".format(fname,self.headers)
        seqd = pdb_edit.sequence_data(pdbin)
        assert len(seqd) == 1,'Currently only support adding data for single chain pdbs'
        
        chain = seqd.keys()[0]
        
        # Add the pdb and chain data
        self.pdbs[idx] = fname
        self.chains[idx] = chain
        
        self.resseqs[idx] = [] # Clear out any existing resseqs
        sequence = seqd[chain][0]
        resseqs = seqd[chain][1]
        # Loop through both sequences
        needle = 0
        GAP='-'
        for res1 in self.sequences[idx]:
            if res1 == GAP:
                self.resseqs[idx].append(None)
            else:
                assert res1.upper() == sequence[needle]
                self.resseqs[idx].append( resseqs[needle] )
                needle += 1
        return True

    def fasta_str(self,pdbname=False):
        if not len(self.sequences): raise RuntimeError,"No sequences have been read!"
        headers=[]
        for i, header in enumerate(self.headers):
            #if pdbname: header = os.path.basename(self.pdbs[i])
            if pdbname: header = ">{0}".format(os.path.basename(self.pdbs[i]))
            headers.append(header)
        return self._fasta_str(headers, self.sequences)
    
    def _fasta_str(self,headers,sequences):
        s=""
        for i, seq in enumerate(sequences):
            s += headers[i]+'\n'
            for chunk in range(0, len(seq), self.MAXWIDTH):
                s += seq[chunk:chunk+self.MAXWIDTH]+"\n"
            s += "\n"
        return s
   
    def from_fasta(self, fasta_file, canonicalise=True, resseq=True):
        with open(fasta_file, "r") as f:
            self._parse_fasta(f, fasta_file=fasta_file, canonicalise=canonicalise)
        if resseq:
            # Add automatically calculated ressegs starting from 1
            #assert len(self.resseqs) == 0,"Altering existing resseqs!"
            for seq in self.sequences:
                self.resseqs.append([])
                for i in range(len(seq)):
                    self.resseqs[-1].append(i+1)
        return
    
    def from_pdb(self, pdbin):
        pdbin_name = os.path.basename(pdbin)
        chain2data = pdb_edit.sequence_data(pdbin)
        assert len(chain2data),"Could not read sequence from pdb: {0}".format(pdbin)
        self.headers = []
        self.sequences = []
        self.pdbs = []
        self.resseqs = []
        self.fasta_files = []
        for chain in sorted(chain2data.keys()):
            seq = chain2data[chain][0]
            resseq = chain2data[chain][1]
            self.headers.append(">From pdb: {0} chain={1} length={2}".format(pdbin_name,chain,len(seq)))
            self.sequences.append(seq)
            self.resseqs.append(resseq)
            self.pdbs.append(pdbin_name)
            self.chains.append(chain)
            self.fasta_files.append(None) # Need to make sure pdb abd fasta arrays have same length
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
    


    def length(self,seq_no=0):
        return len(self.sequences[seq_no])
    
    def numSequences(self):
        return len(self.sequences)
    
    def _parse_fasta(self, fasta, fasta_file=None, canonicalise=True):
        """Parse the fasta file int our data structures & check for consistency 
        Args:
        fasta -- list of strings or open filehandle to read from the fasta file
        """
        headers = []
        sequences = []
        sequence = ""
        header = None
        for line in fasta:
            line = line.strip().rstrip(os.linesep)
            if not line: continue # skip blank lines
            
            if not header:
                if not line.startswith(">"):
                    raise RuntimeError,"FASTA sequencs must be prefixed with a > character: {0}".format(line)
                header = line
                continue
            
            if header and line.startswith(">"):
                headers.append(header)
                sequences.append(sequence)
                header = line
                sequence = ""
            else:
                sequence += line
        
        # Add the last header and sequence
        headers.append(header)
        sequences.append(sequence)
        
        # Now add all the collected data        
        self.headers = []
        self.sequences = []
        self.resseqs = []
        self.fasta_files = []
        self.pdbs = []
        self.chains = []
        for h, s in zip(headers,sequences):
            self.headers.append(h)
            self.sequences.append(s)
            self.fasta_files.append(fasta_file)
            self.resseqs.append(None)
            self.pdbs.append(None)
            self.chains.append(None)
            
        if canonicalise: self.canonicalise()
        return

    def pirStr(self,seqNo=0):
        """Return a canonical MAXWIDTH PIR representation of the file as a line-separated string"""
    
        pirStr = [ self.title + "\n\n" ] # Title plus the obligatory blank line
        
        # Now reformat to MAXWIDTH chars
        for chunk in range(0, len(self.sequances[seqNo]), self.MAXWIDTH):
            pirStr.append(self.sequences[seqNo][ chunk:chunk+self.MAXWIDTH ]+"\n")
            
        # Add last newline
        pirStr.append("\n")
        
        return pirStr

    def sequence(self,seq_no=0):
        return self.sequences[seq_no]

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

    def write_fasta(self,fasta_file,pdbname=False):
        with open(fasta_file,'w') as f:
            for s in self.fasta_str(pdbname=pdbname):
                f.write(s)
        return
    
    def __add__(self,other):
        self.headers += other.headers
        self.sequences += other.sequences
        self.resseqs += other.resseqs
        self.pdbs += other.pdbs
        self.chains += other.chains
        self.fasta_files += other.fasta_files
        return self
    
    def __str__(self):
        return self.__repr__() + "\n" + self.fasta_str()
        
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

    def testAdd(self):
        s1 = Sequence(pdb=os.path.join(self.testfiles_dir,'1GU8.pdb'))
        s2 = Sequence(fasta=os.path.join(self.testfiles_dir,'2uui.fasta'))
        s1 += s2
        
        self.assertTrue(len(s1.sequences),2)
        self.assertTrue(len(s1.resseqs),2)
        self.assertTrue(len(s1.headers),2)
        self.assertTrue(len(s1.pdbs),2)
        self.assertTrue(len(s1.chains),2)
        self.assertTrue(len(s1.fasta_files),2)
        
        # Test write for the hell of it
        out_fasta = 'seq_add.fasta'
        s1.write_fasta(out_fasta)
        os.unlink(out_fasta)
        return

    def testAddPdbData(self):
        fasta1 = os.path.join(self.testfiles_dir,'1ujb_2a6pA_3c7tA.afasta')
        pdbin1 = os.path.join(self.ample_dir, 'examples', 'homologs','1ujb.pdb')
        pdbin2 = os.path.join(self.ample_dir, 'examples', 'homologs','2a6pA.pdb')
        pdbin3 = os.path.join(self.ample_dir, 'examples', 'homologs','3c7tA.pdb')
        s1 = Sequence(fasta=fasta1)
        s1.add_pdb_data(pdbin1)
        s1.add_pdb_data(pdbin2)
        s1.add_pdb_data(pdbin3)
        
        self.assertEqual(s1.pdbs[0], os.path.basename(pdbin1))
        self.assertEqual(s1.chains[0],'A')
        self.assertEqual(s1.pdbs[1], os.path.basename(pdbin2))
        self.assertEqual(s1.chains[1],'A')
        self.assertEqual(s1.pdbs[2], os.path.basename(pdbin3))
        self.assertEqual(s1.chains[2],'A')
        
        p1r = [None, None, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, 
               None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 
               35, 36, 37, 38, 39, 40, 41, 42, 43, 44, None, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, None, None, 70, 71, 72, 73, 
               74, 75, 76, 77, 78, 79, 80, 81, 82, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, 
               None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, 83, None, None, 84, 85, 86, 87, 88, 89, 90, 91, None, 
               92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, None, None, None, None, None, None, 
               None, None, 121, 122, 123, None, None, None, 124, 125, None, None, 126, None, None, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, None, None, 139, 140, 141, 142, 
               143, 144, 145, 146, 147, 148, 149, None, None, 150, 151, None, 152, 153, 154, 155, 156, None, None, None, None, None, None, None, None, None, None]
        
        p2r = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, 
               None, None, None, None, None, None, None, None, None, None, None, 23, 24, 25, 26, 27, 28, None, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 
               50, None, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, None, None, None, None, 74, None, None, 75, None, 76, 77, 78, 79, 80, 81, 82, 
               83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, None, None, None, None, None, None, None, None, 105, 106, 107, 108, 109, 110, 111, 112, 
               113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, None, 133, 134, 135, 136, 137, None, None, 138, 139, 140, 141, 142, 143, 144, 145,
                146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, None, None, None, None, None, None, None, None, 160, 161, 162, 163, None, 164, 165, 166, 167, 168, 169, None, 
                None, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, None, None, None, 185, 186, 187, 188, 189, 190, 191, 192, None, 193, 194, 195, 196, None, None, 
                None, None, None, None, None, None, None, None, None, None, None, None]
 
        
        p3r = [None, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 
               112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, None, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 
               146, 147, 148, 149, None, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 
               180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 
               215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 
               250, None, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 
               284, 285, 286, 287, 288, None, None, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, None, None, 304, None, None, None, 305, 306, 307, 308, 309, 310, 
               311, 312, 313, 314, 315, None, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330]
        
        self.assertEqual(s1.resseqs[0], p1r)
        self.assertEqual(s1.resseqs[1], p2r)
        self.assertEqual(s1.resseqs[2], p3r)
        
        return

    def testFailChar(self):
        
        infasta=""">3HAP:A|PDBID|CHAIN|SEQUENCE
QAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVXAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTI
LAAVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSA
GDGAAATSD"""

        fp = Sequence()
        self.assertRaises(RuntimeError, fp._parse_fasta, infasta.split(os.linesep))    
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

    def testFromPdb(self):
        s1 = Sequence(pdb=os.path.join(self.testfiles_dir,'4DZN.pdb'))
        self.assertEqual(s1.pdbs, ['4DZN.pdb', '4DZN.pdb', '4DZN.pdb'])
        self.assertEqual(s1.chains, ['A', 'B', 'C'])
    
        outfasta=""">From pdb: 4DZN.pdb chain=A length=31
GEIAALKQEIAALKKEIAALKEIAALKQGYY

>From pdb: 4DZN.pdb chain=B length=31
GEIAALKQEIAALKKEIAALKEIAALKQGYY

>From pdb: 4DZN.pdb chain=C length=31
GEIAALKQEIAALKKEIAALKEIAALKQGYY

"""
        self.assertEqual(outfasta, "".join(s1.fasta_str()))
        
        return

    def testResSeq(self):
        pdbin = os.path.join(self.testfiles_dir,'1D7M.pdb')
        s1 = Sequence(pdb=pdbin)
        self.assertTrue(len(s1.sequences),2)
        self.assertTrue(len(s1.headers),2)
        self.assertTrue(len(s1.pdbs),2)
        self.assertEqual(s1.pdbs[0],os.path.basename(pdbin))
        self.assertTrue(s1.resseqs[0][-1],343)
        return
    
    def testAlignFile(self):
        pdbin1 = os.path.join(self.testfiles_dir,'1D7M.pdb')
        pdbin2 = os.path.join(self.testfiles_dir,'1GU8.pdb')
        pdbin3 = os.path.join(self.testfiles_dir,'2UUI.pdb')
        s1 = Sequence(pdb=pdbin1)
        s1 += Sequence(pdb=pdbin2)
        s1 += Sequence(pdb=pdbin3)
        
        ref = """>1D7M.pdb
EMANRLAGLENSLESEKVSREQLIKQKDQLNSLLASLESEGAEREKRLRELEAKLDETLKNLELEKLARMELEARLAKTE
KDRAILELKLAEAIDEKSKLE

>1D7M.pdb
EMANRLAGLENSLESEKVSREQLIKQKDQLNSLLASLESEGAEREKRLRELEAKLDETLKNLELEKLARMELEARLAKTE
KDRAILELKLAEAIDEKSKLE

>1GU8.pdb
VGLTTLFWLGAIGMLVGTLAFAWAGRDAGSGERRYYVTLVGISGIAAVAYVVMALGVGWVPVAERTVFAPRYIDWILTTP
LIVYFLGLLAGLDSREFGIVITLNTVVMLAGFAGAMVPGIERYALFGMGAVAFLGLVYYLVGPMTESASQRSSGIKSLYV
RLRNLTVILWAIYPFIWLLGPPGVALLTPTVDVALIVYLDLVTKVGFGFIALDAAATL

>2UUI.pdb
MHHHHHHKDEVALLAAVTLLGVLLQAYFSLQVISARRAFRVSPPLTTGPPEFERVYRAQVNCSEYFPLFLATLWVAGIFF
HEGAAALCGLVYLFARLRYFQGYARSAQLRLAPLYASARALWLLVALAALGLLAHFLPAALRAALLGRLRTLLPWA

"""
        self.assertEqual(s1.fasta_str(pdbname=True),ref)
        return

#
# Run unit tests
if __name__ == "__main__":
    unittest.main(verbosity=2)
