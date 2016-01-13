#!/usr/bin/env ccp4-python

import glob
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
    pulchra_pdb = 'center_' + os.path.splitext(os.path.basename(pdbout))[0] + '.rebuilt.pdb'
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

def pdb_from_ca(coords, sequence, pdbout):
    root = iotbx.pdb.hierarchy.root()
    model = iotbx.pdb.hierarchy.model()
    root.append_model(model)
    chain = iotbx.pdb.hierarchy.chain(id="A")
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

scwrl_exe = '/opt/scwrl4/Scwrl4'
SCWRL = ample_scwrl.Scwrl(scwrl_exe)

root = '/media/data/shared/TM/quark_models/QUARK_orig'
for pdb_code in  ['1GU8', '2BHW' ,'2BL2', '2EVU', '2O9G', '2UUI', '2WIE', '2X2V',
                  '2XOV', '3GD8', '3HAP', '3LDC', '3OUF', '3PCV', '3RLB', '3U2F', '4DVE']:
    
    if not os.path.isdir(pdb_code): os.mkdir(pdb_code)
    os.chdir(pdb_code)
    
    if pdb_code == ' 3U2F': chain = 'K'
    else: chain = 'A'
    
    count = 0
    fasta_file = os.path.join('/media/data/shared/TM/tm_data',pdb_code,"{0}.fasta".format(pdb_code))
    AS = ample_sequence.Sequence(fasta=fasta_file)
    sequence = AS.sequences[0]
    for idx_tra, tra_file in enumerate(glob.glob(os.path.join(root,"{0}{1}".format(pdb_code,chain),'rep*.tra*'))):
        coord_list = read_tra(tra_file)
        for idx_coord, coords in enumerate(coord_list):
            assert len(sequence) == len(coords)
            count += 1
            pdbca = "{0}{1}_{2}.pdb".format(pdb_code, chain, count)
            if not pdb_from_ca(coords, sequence, pdbca):
                raise RuntimeError("Error processing tra_file {0} count {1}".format(tra_file, count))
        
            pdb_side = ample_util.filename_append(pdbca,"pulchra")
            if not add_sidechains_pulchra(pdbca, pdb_side):
                print "Failed to add pulchra sidechains for tra_file {0} count {1}".format(tra_file, count)
                pdb_side = ample_util.filename_append(pdbca,"maxsprout")
                if not add_sidechains_maxsprout(pdbca, pdb_side):
                    print "Failed to add maxsprout sidechains for tra_file {0} count {1}".format(tra_file, count)
                    continue
            
            # Now add sidechains with SCWRL
            pdbout = ample_util.filename_append(pdb_side,"scwrl")
            SCWRL.add_sidechains(pdbin=pdb_side, pdbout=pdbout, hydrogens=False, strip_oxt=False)
            
            
            


