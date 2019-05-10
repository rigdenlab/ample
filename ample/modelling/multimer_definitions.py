DIMER = 'dimer'
TRIMER = 'trimer'
TETRAMER = 'tetramer'

SYMFILE_DIMER = """symmetry_name c3
subunits 3
recenter
number_of_interfaces  1
E = 3*VRT0001 + 3*(VRT0001:VRT0002)
anchor_residue COM
virtual_transforms_start
start -1,0,0 0,1,0 0,0,0
rot Rz 3
virtual_transforms_stop
connect_virtual JUMP1 VRT0001 VRT0002
connect_virtual JUMP2 VRT0002 VRT0003
set_dof BASEJUMP x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)
"""


SYMFILE_TRIMER = """symmetry_name c3
subunits 3
recenter
number_of_interfaces  1
E = 3*VRT0001 + 3*(VRT0001:VRT0002)
anchor_residue COM
virtual_transforms_start
start -1,0,0 0,1,0 0,0,0
rot Rz 3
virtual_transforms_stop
connect_virtual JUMP1 VRT0001 VRT0002
connect_virtual JUMP2 VRT0002 VRT0003
set_dof BASEJUMP x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)
"""

BROKER_SETUP_STR = """CLAIMER FoldandDockClaimer
END_CLAIMER
"""

FLAGSFILE_STR = """-run:protocol broker 
-broker:setup {broker_file} 
-symmetry:symmetry_definition {symdef_file} 
-symmetry:initialize_rigid_body_dofs 
-fold_and_dock:rigid_body_cycles 1 
-fold_and_dock:rigid_body_frequency 5 
-fold_and_dock:rotate_anchor_to_x 
-run:reinitialize_mover_for_each_job 
-score:weights score13_env_hb 
-abinitio::recover_low_in_stages 0 
-abinitio:rg_reweight 0.001 
-abinitio:use_filters false 
-packing:ex1 
-packing:ex1:level 1 
-packing:ex2 
-packing:ex2:level 1 
-packing:extrachi_cutoff 0 
-in:file:fasta {fasta_file} 
-in:file:frag3 {frags3} 
-in:file:frag9 {frags9}
-evaluation:rmsd_column _ 
-evaluation:symmetric_rmsd 
#-out:nstruct 1000 
-relax:quick 
-relax:jump_move 
#-out:path:pdb ./ 
-out:pdb 
-mute core.io.database 
-out:file:scorefile score.sc 
-constraints:cst_file {constraint_file} 
-constraints:cst_fa_file {constraint_file} 
-constraints:cst_weight 10 
-constraints:cst_fa_weight 10
"""

