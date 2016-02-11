#!/usr/bin/env ccp4-python

import glob
import logging
import os
import sys
sys.path.append('/opt/ample-dev1/python')

import ample_sequence
import ample_scwrl
import ample_util

import iotbx.pdb
from iotbx.pdb import amino_acid_codes

def add_sidechains_pulchra(pdbin, pdbout):
    pulchra_exe = '/media/data/shared/TM/quark_models/pulchra_3.06/pulchra.pl'
    
    # The pdb file that pulchra will produce
    #pulchra_pdb = 'center_' + os.path.splitext(os.path.basename(pdbout))[0] + '.rebuilt.pdb'
    pulchra_pdb = 'pul_' + pdbin
    cmd = [ pulchra_exe, pdbin ]
    logfile = 'pulchra.log'
    rtn = ample_util.run_command(cmd, logfile=logfile)
    if rtn != 0:
        #raise RuntimeError("{0} failed".format(cmd))
        if os.path.isfile(pulchra_pdb): os.unlink(pulchra_pdb)
        return None
    
    # No dodgy return code, but no output file
    if not os.path.isfile(pulchra_pdb): return None
    
    # Tt worked!
    os.rename(pulchra_pdb, pdbout)
    return pdbout

def add_sidechains_maxsprout(pdbin, pdbout):
    maxsprout_exe = '/media/data/shared/TM/quark_models/maxsprout_2006.10/maxsprout.pl'

    cmd = [ maxsprout_exe, pdbin, pdbout ]
    logfile = 'maxsprout.log'
    rtn = ample_util.run_command(cmd, logfile=logfile)
    if rtn != 0:
        #raise RuntimeError("{0} failed".format(cmd))
        if os.path.isfile(pdbout): os.unlink(pdbout)
        return None
    
    # No dodgy return code, but no output file
    if not os.path.isfile(pdbout): return None
    
    return pdbout

def coord_to_pdb(coords, sequence, coord_id, clean=True):
    assert len(sequence) == len(coords)
    pdbca = "{0}.pdb".format(coord_id)
    if not pdb_from_ca(coords, sequence, pdbca):
        raise RuntimeError("Error processing coords {0}".format(coord_id))
    #print("Wrote file {0}".format(pdbca))
    
    pdb_side = ample_util.filename_append(pdbca,"pulchra")
    if not add_sidechains_pulchra(pdbca, pdb_side):
        print "Failed to add pulchra sidechains for coord {0}".format(coord_id)
        pdb_side = ample_util.filename_append(pdbca,"maxsprout")
        if not add_sidechains_maxsprout(pdbca, pdb_side):
            print "*** Failed to add maxsprout sidechains for coord {0}".format(coord_id)
            return None
    
    # Now add sidechains with SCWRL
    pdbout = ample_util.filename_append(pdb_side,"scwrl")
    SCWRL.add_sidechains(pdbin=pdb_side, pdbout=pdbout, hydrogens=False, strip_oxt=False)
    
    if clean: map(os.unlink,[pdbca, pdb_side])
    
    return pdbout

def pdb_from_ca(coords, sequence, pdbout, chain):
    root = iotbx.pdb.hierarchy.root()
    model = iotbx.pdb.hierarchy.model()
    root.append_model(model)
    chain = iotbx.pdb.hierarchy.chain(id=chain)
    model.append_chain(chain)
    
    for i, aa1 in enumerate(sequence):
        resseq = "{0: 4}".format(i+1)
        # Needs to be 4 characters
        rg = iotbx.pdb.hierarchy.residue_group(resseq=resseq)
        chain.append_residue_group(residue_group=rg)
        aa3 = amino_acid_codes.three_letter_given_one_letter[aa1]
        ag = iotbx.pdb.hierarchy.atom_group(resname=aa3)
        rg.append_atom_group(atom_group=ag)
        
        atom = iotbx.pdb.hierarchy.atom()
        atom.set_element("C")
        atom.set_name(" CA ")
        atom.set_serial(resseq)
        atom.set_xyz(coords[i])
        ag.append_atom(atom)
    
    root.write_pdb_file(pdbout)
    return pdbout
      
def read_tra(tra_file):
    coords = []
    with open(tra_file) as f:
        for i, line in enumerate(f.readlines()):
            line = line.strip()
            if not line: continue
            fields = line.split()
            if len(fields) == 3:
                try:
                    coords[-1].append(map(float, fields))
                except Exception as e:
                    raise RuntimeError("Error reading trajectory: {0} line {1}: {2}".format(tra_file, i, e))
            else:
                coords.append([])
    return coords

#logging.basicConfig(level=logging.DEBUG)
scwrl_exe = '/opt/scwrl4/Scwrl4'
SCWRL = ample_scwrl.Scwrl(scwrl_exe)

NUM_MODELS = 5000
start_dir = os.path.abspath(os.getcwd())
root = '/media/data/shared/TM/quark_models/QUARK_orig'
for pdb_code in  ['1GU8', '2BHW' ,'2BL2', '2EVU', '2O9G', '2UUI', '2WIE', '2X2V',
                  '2XOV', '3GD8', '3HAP', '3LDC', '3OUF', '3PCV', '3RLB', '3U2F', '4DVE']:
    
    print("Processing pdb {0}".format(pdb_code))
    directory = os.path.join(start_dir,pdb_code)
    if not os.path.isdir(directory): os.mkdir(directory)
    os.chdir(directory)
    
    chain = 'A'
    if pdb_code == '3U2F': chain = 'K'
    
    count = 0
    fasta_file = os.path.join('/media/data/shared/TM/tm_data',pdb_code,"{0}.fasta".format(pdb_code))
    AS = ample_sequence.Sequence(fasta=fasta_file)
    sequence = AS.sequences[0]
    
    # Loop through all trajectory files and get a list of all coordinates
    coords = []
    for idx_tra, tra_file in enumerate(glob.glob(os.path.join(root,"{0}{1}".format(pdb_code,chain),'rep*.tra*'))):
        print("Processing tra_file {0}".format(tra_file))
        coords += read_tra(tra_file)
    
    # Work out stride so that we end up with 5000 models
    ncoords = len(coords)
    #stride = int(round( float(ncoords) / float(5000) ))
    stride = ncoords / NUM_MODELS
    
    count = 0
    for i in range(0, ncoords, stride):
        coord_id = "{0}{1}_{2}".format(pdb_code, chain, i)
        print "Making pdb ",coord_id
        coord_to_pdb(coords[i], sequence, coord_id, clean=True)
        count += 1
        if count >= NUM_MODELS: break

    
            



"""
Old code from unpack_quark.py

headers = []
coord_list = []
with open(fn) as f:
    while True:
        # read header
        header = f.readline().strip()
        if not header: break
        headers.append(header)
        fields = header.strip().split()
        nca = int(fields[0])
        ca = []
        for _ in range(nca):
            x, y, z = f.readline().strip().split()
            ca.append((float(x),float(y),float(z)))
        coord_list.append(ca)
        
assert len(headers) == len(coord_list)


with open(fasta) as f:
    f.readline()
    seq=""
    while True:
        s = f.readline().strip()
        if not s: break
        seq += s

seq3 = [ one2three[s] for s in seq ]
assert len(seq3) == len(coord_list[0])

bname = os.path.splitext(os.path.basename(fn))[0]

for i in range(len(headers)):
    fname = "{0}_{1}.pdb".format(bname,i)
    with open(fname,'w') as f:
        f.write("REMARK   3 QUARK MODEL {0} FROM FILE: {1}\n".format(i,fn))
        f.write("REMARK   3 {0}\n".format(headers[i]))
        f.write("MODEL 1\n")
        for j, coord in enumerate(coord_list[i]):
            x,y,z = coord
            f.write("ATOM  {0:5d}  CA  {1} A {2:3d}    {3:8.3F}{4:8.3F}{5:8.3F}  0.50 30.00           C\n".format(int(j),seq3[j],int(j),x,y,z))
        #           "ATOM      1  N   VAL A   2      35.075  18.239 -14.019  0.50 41.75           N"
        f.write("TER\n")
        f.write("ENDMDL\n\n")

"""
