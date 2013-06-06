#!/usr/bin/env python
'''
Created on 30 May 2013

@author: jmht

Useful stuff for PDBs - currently just remove HETATM lines
'''

import os
import sys

inpdb = sys.argv[1]

outpdb=None
if len(sys.argv) == 3:
    outpdb = sys.argv[2]
    
if not outpdb:
    name = os.path.splitext( os.path.basename(inpdb) )[0]
    dirname = os.path.dirname( os.path.abspath( inpdb ) )
    outpdb = os.path.join( dirname, name + "_clean.pdb" )
    
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
    
    
    
