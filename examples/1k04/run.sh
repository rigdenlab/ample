#!/bin/bash

# 
# Script to run missing domains test case
#
# Change the below if the path to these external dependencies of AMPLE
# cannot be found in your standard PATH
ROSETTA_DIR=/opt/rosetta3.4

export PATH=\
/opt/spicker:\
/opt/theseus_src:\
/opt/scwrl4:\
/opt/maxcluster:\
/opt/shelx:\
$PATH

# This is only required if you are running ample from outside
# the standard CCP4 directory
AMPLEDIR=/opt/ample-dev1
export PYTHONPATH=${AMPLEDIR}/python:$PYTHONPATH

${AMPLEDIR}/bin/ample \
      -mtz ${PWD}/1k04_cad.mtz \
      -domain_all_chains_fasta ${PWD}/1k04_.fasta  \
      -domain_all_chains_pdb  ${PWD}/Known_40.pdb \
      -frags3mers ${PWD}/aa1k04_03_05.200_v1_3 \
      -frags9mers ${PWD}/aa1k04_09_05.200_v1_3 \
      -rosetta_dir $ROSETTA_DIR \
      -make_frags False \
      -make_models True \
      -CLUSTER False \
      -shelx_cycles 10 \
      -nmodels 30 \
      -early_terminate True  \
      -percent 5 \
      -use_shelxe True \
      -use_arpwarp False \
      -use_buccaneer False \
      -molreponly True \
      -nproc 1 \
      -num_clusters 1 \
      -allatom True \
