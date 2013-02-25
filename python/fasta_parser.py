#!/usr/bin/python2.6

import re

def parse_fasta(fasta, outfasta):
    """
    Reformat the fasta file
    """
    aa=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    
    fastin=open(fasta)
    seq = ''
    fasout=open(outfasta, "w")
    
    for line in fastin:
        if re.match('^>', line):
            fasout.write(line)
        if not re.match('^>', line):
            seq+=line
    
    seqout = ''
    for x in seq:
        if x.upper() in aa:
            seqout+=x.upper()
        if re.match('\n', x):
            seqout+=x
          
    fasout.write(seqout+'\n')
    fasout.close()

#fasta = '/home/jaclyn/DOMAINS/ample/TEST/1al6/1al6_.fasta'
#outfasta = '/home/jaclyn/DOMAINS/ample/TEST/ROSETTA_MR_1/fas'
#parse_fasta(fasta, outfasta)




