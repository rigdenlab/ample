#!/usr/bin/env ccp4-python
import logging
import os

import exit_util
import pdb_edit

class Sequence(object):
    """A class to handle a fasta file"""
    
    def __init__(self, fasta=None, pdb=None, canonicalise=False):
        """Initialise the object"""
        
        self.MAXWIDTH=80 # maximum width of any line 
        self.name = None # identifier for this sequence object
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
        self.name = os.path.splitext(os.path.basename(fasta_file))[0]
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
        self.name = os.path.splitext(pdbin_name)[0]
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

def process_fasta(amoptd):
    # Check we can find the input fasta
    if not os.path.exists(str(amoptd['fasta'])):
        msg = 'Cannot find fasta file: {0}'.format(amoptd['fasta'])
        exit_util.exit_error(msg)
    
    # Reformat to what we need
    logging.debug('Parsing FASTA file')
    try: fp = Sequence(fasta=amoptd['fasta'])
    except Exception as e:
        msg = "Error parsing FASTA file: {0}\n\n{1}".format(amoptd['fasta'],e.message)
        exit_util.exit_error(msg)
    if fp.numSequences() != 1:
        msg = "ERROR! Fasta file {0} has > 1 sequence in it.".format(amoptd['fasta'])
        exit_util.exit_error(msg)
    
    # Length checks
    amoptd['fasta_length'] = fp.length()
    logging.info("Fasta is {0} amino acids long".format(amoptd['fasta_length']))
    
    # Check we have a decent length
    if amoptd['fasta_length'] < 9:
        msg = "ERROR! Fasta is of length {0}. This is much too short!".format(amoptd['fasta_length'])
        exit_util.exit_error(msg)
    
    # Check we will be able to truncate at this level
    if (float(amoptd['fasta_length']) / 100) * float(amoptd['percent']) < 1:
        msg = "Cannot truncate a fasta sequence of length {0} with {1} percent intervals. Please select a larger interval.".format(amoptd['fasta_length'], amoptd['percent'])
        exit_util.exit_error(msg)
        
    # Check that the sequence doesn't have a his-tag in it
    if not amoptd['allow_his_tag']:
        his_tag = 'HHHHHH'
        i = fp.sequence().find(his_tag)
        l = fp.length()
        if (0 <= i <= 20) or (l-20 <= i <= l):
            msg = 'The fasta sequence contains a his tag sequence {0} at position {1}. If you wish to use ample with this sequence, please use the \"-allow_his_tag True\" option'.format(his_tag,i)
            exit_util.exit_error(msg)
    
    # Fasta is ok, so write out a canonical fasta in the work directory
    outfasta = os.path.join(amoptd['work_dir'], amoptd['name'] + '_.fasta')
    fp.write_fasta(outfasta)
    amoptd['fasta'] = outfasta
    amoptd['sequence'] = fp.sequence()

    return
        
