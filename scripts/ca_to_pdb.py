#!/usr/bin/env ccp4-python

import os, sys
sys.path.append('/opt/ample-dev1/python')

import ample_sequence
import ample_scwrl
import ample_util
import iotbx.pdb
from iotbx.pdb import amino_acid_codes

cfile = '1GU8A.rep30.tra10D.1'
fasta_file = '1GU8.fasta'
cfile = '2UUI.rep30.tra1D.1'
fasta_file = '2UUI.fasta'

coords = []
with open(cfile) as f:
    for i, line in enumerate(f.readlines()):
        if i == 0 : continue
        line = line.strip()
        if not line: continue
        coords.append(map(float, line.split()))
        
# with open('jens.xyz', 'w') as f:
#     f.write('{0}\n'.format(len(coords)))
#     f.write('test\n')
#     for c in coords:
#         f.write('C    {0}   {1}   {2}\n'.format(c[0], c[1], c[2]))


AS = ample_sequence.Sequence(fasta=fasta_file)

assert len(AS.sequences[0]) == len(coords)

root = iotbx.pdb.hierarchy.root()
model = iotbx.pdb.hierarchy.model()
root.append_model(model)
chain = iotbx.pdb.hierarchy.chain(id="A")
model.append_chain(chain)

for i, aa1 in enumerate(AS.sequences[0]):
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
    

pdbout = cfile+'.pdb'
root.write_pdb_file(pdbout)


pulchra_exe = '/media/data/shared/TM/quark_models/pulchra_3.06/pulchra.pl'

cmd = [ pulchra_exe, pdbout ]
logfile = 'pulchra.log'
rtn = ample_util.run_command(cmd, logfile=logfile)
if rtn != 0: raise RuntimeError("{0} failed".format(cmd))

pulchra_pdb = 'center_' + os.path.splitext(os.path.basename(pdbout))[0] + '.rebuilt.pdb'
if not os.path.isfile(pulchra_pdb): raise RuntimeError("Missing pulchra output file")

maxsprout_exe = '/media/data/shared/TM/quark_models/maxsprout_2006.10/maxsprout.pl'
maxsprout_pdb = ample_util.filename_append(pdbout,"maxsprout")
cmd = [ maxsprout_exe, pdbout, maxsprout_pdb ]
logfile = 'maxsprout.log'
rtn = ample_util.run_command(cmd, logfile=logfile)
if rtn != 0: raise RuntimeError("{0} failed".format(cmd))

if not os.path.isfile(maxsprout_pdb): raise RuntimeError("Missing maxsprout output file")


scwrl_exe = '/opt/scwrl4/Scwrl4'
SCWRL = ample_scwrl.Scwrl(scwrl_exe)
SCWRL.add_sidechains(pdbin=None, pdbout=None, hydrogens=False, strip_oxt=True)

