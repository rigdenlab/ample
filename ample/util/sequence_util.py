#!/usr/bin/env ccp4-python
from collections import OrderedDict
import logging
import os

import iotbx.pdb

from ample.util import ample_util
from ample.util import exit_util


def sequence(pdbin):
    return _sequence(iotbx.pdb.pdb_input(pdbin).construct_hierarchy())


def _sequence(hierarchy):
    """Extract the sequence of residues from a pdb file."""
    chain2data = _sequence_data(hierarchy)
    return dict((k, chain2data[k][0]) for k in chain2data.keys())


def _sequence1(hierarchy):
    """Return sequence of the first chain"""
    d = _sequence(hierarchy)
    return d[sorted(d.keys())[0]]


def sequence_data(pdbin):
    return _sequence_data(iotbx.pdb.pdb_input(pdbin).construct_hierarchy())


def _sequence_data(hierarchy):
    """Extract the sequence of residues and resseqs from a pdb file."""
    chain2data = OrderedDict()
    for chain in hierarchy.models()[0].chains():  # only the first model
        if not chain.is_protein():
            continue
        seq, resseq = chain_data(chain)
        chain2data[chain.id] = (seq, resseq)
    return chain2data


def chain_data(chain):
    seq = ""
    resseq = []
    for residue in chain.conformers()[0].residues():  # Just look at the first conformer
        # See if any of the atoms are non-hetero - if so we add this residue
        if any([not atom.hetero for atom in residue.atoms()]):
            seq += ample_util.three2one[residue.resname]
            resseq.append(residue.resseq_as_int())
    return seq, resseq


def chain_sequence(chain):
    if not chain.is_protein():
        return None
    return chain_data(chain)[0]


class Sequence(object):
    """A class to handle a fasta file"""

    def __init__(self, fasta=None, pdb=None, canonicalise=False):
        """Initialise the object"""

        self.MAXWIDTH = 80  # maximum width of any line
        self.name = None  # identifier for this sequence object
        self.headers = []  # title lines
        self.sequences = []  # The fasta sequences (just AA) as a single string
        self.resseqs = []  # The fasta sequences (just AA) as a single string
        self.pdbs = []  # pdb files that sequences were derived from
        self.chains = []  # The chain that the sequence belongs too (if applicable)
        self.fasta_files = []  # The fasta files the sequence were read from

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
        if not os.path.isfile(pdbin):
            raise RuntimeError("Cannot find provided pdb file: {}".format(pdbin))

        if len(self.headers) < 1 and len(self.sequences) < 1:
            raise ValueError("Require at least one one sequence.")

        fname = os.path.basename(pdbin)
        name, ext = os.path.splitext(fname)
        if ext != ".pdb":
            raise TypeError("PDB files must have extension .pdb")

        got = False
        for idx, h in enumerate(self.headers):
            n = h[1:].split('.')[0]
            if n == name:
                got = True
                break

        if not got:
            raise RuntimeError("Could not find matching pdb name for {0} in headers {1}".format(fname, self.headers))
        seqd = sequence_data(pdbin)
        if len(seqd) > 1:
            raise RuntimeError("Currently only support adding data for single chain pdbs")

        chain = seqd.keys()[0]

        self.pdbs[idx] = fname
        self.chains[idx] = chain

        self.resseqs[idx] = []
        sequence = seqd[chain][0]
        resseqs = seqd[chain][1]

        needle = 0
        GAP = '-'
        for res1 in self.sequences[idx]:
            if res1 == GAP:
                self.resseqs[idx].append(None)
            else:
                assert res1.upper() == sequence[needle]
                self.resseqs[idx].append(resseqs[needle])
                needle += 1
        return True

    def fasta_str(self, pdbname=False):
        if len(self.sequences) < 1:
            raise RuntimeError("No sequences have been read!")
        headers = []
        for i, header in enumerate(self.headers):
            if pdbname:
                header = ">{0}".format(os.path.basename(self.pdbs[i]))
            headers.append(header)
        return self._fasta_str(headers, self.sequences)

    def _fasta_str(self, headers, sequences):
        s = ""
        for i, seq in enumerate(sequences):
            s += headers[i] + os.linesep
            for chunk in range(0, len(seq), self.MAXWIDTH):
                s += seq[chunk : chunk + self.MAXWIDTH] + os.linesep
            s += os.linesep
        return s

    def from_fasta(self, fasta_file, canonicalise=True, resseq=True):
        self.name = os.path.splitext(os.path.basename(fasta_file))[0]
        with open(fasta_file, "r") as f:
            self._parse_fasta(f, fasta_file=fasta_file, canonicalise=canonicalise)
        if resseq:
            for i, seq in enumerate(self.sequences):
                if len(self.resseqs) >= i + 1 and self.resseqs[i] is None:
                    self.resseqs[i] = []
                    for j in range(len(seq)):
                        self.resseqs[i].append(j + 1)

    def from_pdb(self, pdbin):
        pdbin_name = os.path.basename(pdbin)
        self.name = os.path.splitext(pdbin_name)[0]
        chain2data = sequence_data(pdbin)
        if len(chain2data) < 1:
            raise RuntimeError("Could not read sequence from pdb: {0}".format(pdbin))
        self.headers = []
        self.sequences = []
        self.pdbs = []
        self.resseqs = []
        self.fasta_files = []
        for chain in sorted(chain2data.keys()):
            seq = chain2data[chain][0]
            resseq = chain2data[chain][1]
            self.headers.append(">From pdb: {0} chain={1} length={2}".format(pdbin_name, chain, len(seq)))
            self.sequences.append(seq)
            self.resseqs.append(resseq)
            self.pdbs.append(pdbin_name)
            self.chains.append(chain)
            self.fasta_files.append(None)

    def canonicalise(self):
        """Reformat the fasta file
        
        Description
        -----------
        Needed because Rosetta has problems reading fastas. For it to be read, it has to have no spaces in the sequence,
        a name that is 4 characters and an underscore (ABCD_), everything has to be uppercase, and there has to be a
        return carriage at the end - this has to be linux formatted as when I make a fasta in windows, it doesnt recognize
        the return carriage.

        Rosetta has a lot of problems with fastas so we put in this script to deal with it.

        """
        aa = list("ACDEFGHIKLMNPQRSTVWY")
        for i, seq in enumerate(self.sequences):
            cs = ""
            for char in seq:
                char = char.upper()
                if char in aa:
                    cs += char
                else:
                    msg = "Non-standard amino acid in your sequence: '{}'. Please format to standard amino acids!"
                    raise RuntimeError(msg.format(char))
            self.sequences[i] = cs

    def length(self, seq_no=0):
        return len(self.sequences[seq_no])

    def mutate_residue(self, from_aa, res_seq, to_aa, seq_id=0):
        """Change residue type from_aa at position res_seq to be of type to_aa
        
        Note: res_seq positions start from 1 not zero!

        Parameters
        ----------
        from_aa : str
           Single-letter amino acid
        res_seq : int
           Residue sequence number
        to_aa : str
           Single-letter amino acid
        seq_id : int
           The index of the sequence to operate on (counting from zero)
        
        """
        seq = self.sequences[seq_id]
        assert seq[res_seq - 1] == from_aa, "Amino acid at position {0} is {1} not {2}".format(
            res_seq, seq[res_seq - 1], from_aa
        )
        self.sequences[seq_id] = seq[: res_seq - 1] + to_aa + seq[res_seq:]
        return

    def numSequences(self):
        return len(self.sequences)

    def _parse_fasta(self, fasta, fasta_file=None, canonicalise=True):
        """Parse the fasta file int our data structures & check for consistency"""
        headers = []
        sequences = []
        sequence = ""
        header = None
        for line in fasta:
            line = line.strip().rstrip(os.linesep)
            if not line:
                continue

            if header is None:
                if line.startswith(">"):
                    header = line
                    continue
                else:
                    raise RuntimeError("FASTA sequencs must be prefixed with a > character: {}".format(line))

            if header and line.startswith(">"):
                headers.append(header)
                sequences.append(sequence)
                header = line
                sequence = ""
            else:
                sequence += line

        headers.append(header)
        sequences.append(sequence)

        logging.debug("Stripping whitespace characters from sequence")
        logging.debug("Stripping '*' character of end of sequence")
        for i in range(len(sequences)):
            seq = sequences[i].replace(" ", "")
            if seq.endswith('*'):
                seq = seq[:-1]
            sequences[i] = seq

        self.headers, self.sequences, self.resseqs, self.fasta_files, self.pdbs, self.chains = [], [], [], [], [], []
        for h, s in zip(headers, sequences):
            self.headers.append(h)
            self.sequences.append(s)
            self.fasta_files.append(fasta_file)
            self.resseqs.append(None)
            self.pdbs.append(None)
            self.chains.append(None)

        if canonicalise:
            self.canonicalise()

    def pirStr(self, seqNo=0):
        """Return a canonical MAXWIDTH PIR representation of the file as a line-separated string"""
        pirStr = [self.title + os.linesep + os.linesep]
        for chunk in range(0, len(self.sequances[seqNo]), self.MAXWIDTH):
            pirStr.append(self.sequences[seqNo][chunk : chunk + self.MAXWIDTH] + os.linesep)
        pirStr.append(os.linesep)
        return pirStr

    def sequence(self, seq_no=0):
        return self.sequences[seq_no]

    def toPir(self, input_fasta, output_pir=None):
        """Take a fasta file and output the corresponding PIR file"""
        self.parseFasta(input_fasta)
        if not output_pir:
            dir, ffile = os.path.split(input_fasta)
            fname = os.path.splitext(ffile)[0]
            output_pir = os.path.join(dir, fname + ".pir")
        with open(output_pir, "w") as pirout:
            for line in self.pirStr():
                pirout.write(line)
        return

    def write_fasta(self, fasta_file, pdbname=False):
        with open(fasta_file, 'w') as f:
            for s in self.fasta_str(pdbname=pdbname):
                f.write(s)
        return

    def __add__(self, other):
        self.headers += other.headers
        self.sequences += other.sequences
        self.resseqs += other.resseqs
        self.pdbs += other.pdbs
        self.chains += other.chains
        self.fasta_files += other.fasta_files
        return self

    def __str__(self):
        return self.__repr__() + os.linesep + self.fasta_str()


def process_fasta(amoptd, canonicalise=False):
    # Check we can find the input fasta
    if not os.path.exists(str(amoptd['fasta'])):
        msg = 'Cannot find fasta file: {0}'.format(amoptd['fasta'])
        exit_util.exit_error(msg)

    # Reformat to what we need
    logging.debug('Parsing FASTA file')
    try:
        fp = Sequence(fasta=amoptd['fasta'], canonicalise=canonicalise)
    except Exception as e:
        msg = "Error parsing FASTA file: {0}\n\n{1}".format(amoptd['fasta'], e.message)
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
        msg = "Cannot truncate a fasta sequence of length {0} with {1} percent intervals. Please select a larger interval.".format(
            amoptd['fasta_length'], amoptd['percent']
        )
        exit_util.exit_error(msg)

    # Check that the sequence doesn't have a his-tag in it
    if not amoptd['allow_his_tag']:
        his_tag = 'HHHHHH'
        i = fp.sequence().find(his_tag)
        l = fp.length()
        if (0 <= i <= 20) or (l - 20 <= i <= l):
            msg = 'The fasta sequence contains a his tag sequence {0} at position {1}. If you wish to use ample with this sequence, please use the \"-allow_his_tag True\" option'.format(
                his_tag, i
            )
            exit_util.exit_error(msg)

    # Fasta is ok, so write out a canonical fasta in the work directory
    outfasta = os.path.join(amoptd['work_dir'], amoptd['name'] + '_.fasta')
    fp.write_fasta(outfasta)
    amoptd['fasta'] = outfasta
    amoptd['sequence'] = fp.sequence()

    return
