#!/usr/bin/env python
'''
Created on 30 May 2013

@author: jmht

Useful stuff for PDBs - currently just remove HETATM lines
'''


import os
import sys

import pdb_edit

three2one = {
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
#aaDict.update( dict((v, k) for (k, v) in aaDict.items()) )
one2three =  dict((v, k) for (k, v) in three2one.items()) 


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
    
    return
    


def keep_matching( refpdb, targetpdb, outpdb):
    """Create a new pdb file that only contains that atoms in targetpdb that are
    also in refpdb
    
    It renumbers the atoms as we go
    
    NB: Assumes that both pdb files only contain one model
    
    Args:
    refpdb: path to pdb that contains the minimal set of atoms we want to keep
    targetpdb: path to the pdb that will be stripped of non-matching atoms
    outpdb: output path for the stripped pdb
    """


    # Go through refpdb and find which residues are present
    f = open(refpdb, 'r')
    
    # map of resSeq to list of PdbAtom objects
    residues = {}
    
    last = None
    for line in f:
        if line.startswith("MODEL"):
            raise RuntimeError, "Multi-model file!"
        
        if line.startswith("ATOM"):
            a = pdb_edit.PdbAtom( line )
            if a.resSeq != last:
                if a.resSeq in residues:
                    raise RuntimeError,"Multiple chains in pdb - found residue #: {0} again.".format(a.resSeq)
                last = a.resSeq
                residues[ last ] = [ a ]
            else:
                residues[ last ].append( a )
                
    f.close()
    
    # Now read in target pdb and output everything bar the atoms in this file that
    # don't match those in the refpdb
    t = open(targetpdb,'r')
    out = open(outpdb,'w')
    
    for line in t:
        if line.startswith("MODEL"):
            raise RuntimeError, "Multi-model file!"
        
        # Stop at TER
        if line.startswith("TER"):
            break
        
        if line.startswith("ATOM"):
            
            atom = pdb_edit.PdbAtom( line )
            
            # Skip any residues that don't match
            if atom.resSeq not in residues:
                continue
            
            # Skip any atoms that aren't in the reference PDB
            got=False
            for i, ratom in enumerate(residues[ atom.resSeq ]):
                if atom.name == ratom.name:
                    residues[ atom.resSeq ].pop(i)
                    got=True
                    break
                
            if not got:
                continue
            
        # Output everything else
        out.write(line)
        
    t.close()
    out.close()
    
    return


refpdb = "/Users/jmht/Documents/AMPLE/data/test/3PCV/refmac_phaser_loc0_ALL_poly_ala_trunc_2.822761_rad_1_UNMOD.pdb"
refpdb = "/Users/jmht/Documents/AMPLE/data/test/3PCV/refmac_phaser_loc0_ALL_poly_ala_trunc_2.822761_rad_1_UNMOD_1CHAIN_2.pdb"
#refpdb = "/Users/jmht/Documents/AMPLE/data/test/3PCV/ref.pdb"
targetpdb = "/Users/jmht/Documents/AMPLE/data/test/3PCV/3PCV_clean.pdb"
#targetpdb = "/Users/jmht/Documents/AMPLE/data/test/3PCV/target.pdb"
outpdb = "/Users/jmht/Documents/AMPLE/data/test/3PCV/jens.pdb"

keep_matching(refpdb, targetpdb, outpdb)
#strip_hetatm("/Users/jmht/Documents/AMPLE/data/test/3PCV/3PCV_clean.pdb","/Users/jmht/Documents/AMPLE/data/test/3PCV/3PCV_clean_htm.pdb")

