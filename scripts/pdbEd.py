#!/usr/bin/env python
'''
Created on 30 May 2013

@author: jmht

Useful stuff for PDBs - currently just remove HETATM lines
'''

import os
import sys

import pdb_edit

# inpdb = sys.argv[1]
# 
# outpdb=None
# if len(sys.argv) == 3:
#     outpdb = sys.argv[2]
#     
# if not outpdb:
#     name = os.path.splitext( os.path.basename(inpdb) )[0]
#     dirname = os.path.dirname( os.path.abspath( inpdb ) )
#     outpdb = os.path.join( dirname, name + "_clean.pdb" )


def strip_hetatm( inpdb, outpdb):
    """Remove all hetatoms from pdbfile"""
    o = open( outpdb, 'w' )
    
    hremoved=-1
    for i, line in enumerate( open(inpdb) ):
        
        # Remove EOL
        line = line.rstrip( "\n" )
        
        # Remove any HETATOM lines and following ANISOU lines
        if line.startswith("HETATM"):
            hremoved = i
            continue
        
        if line.startswith("ANISOU") and i == hremoved+1:
            continue
        
        o.write( line + "\n" )
        
        
    o.close()
    

aaDict = {
    'ALA' : 'A',    
    'ARG' : 'R',    
    'ASN' : 'N',    
    'ASP' : 'D',    
    'CYS' : 'C',    
    'GLU' : 'E',    
    'GLN' : 'Q',    
    'GLY' : 'G',    
    'HIS' : 'H',    
    'ILE' : 'I',    
    'LEU' : 'L',    
    'LYS' : 'K',    
    'MET' : 'M',    
    'PHE' : 'F',    
    'PRO' : 'P',    
    'SER' : 'S',    
    'THR' : 'T',    
    'TRP' : 'W',    
    'TYR' : 'Y',   
    'VAL' : 'V,'
}

# http://stackoverflow.com/questions/3318625/efficient-bidirectional-hash-table-in-python
aaDict.update( dict((v, k) for (k, v) in aaDict.items()) )




sys.exit()
    
    
    
o = open(inpdb)

sequence = []
last = None
for line in o:
    if line.startswith("ATOM"):
        a = pdb_edit.PdbAtom( line )
        if a.resSeq != last:
            last = a.resSeq
            sequence.append( a.resName )
            
print sequence
