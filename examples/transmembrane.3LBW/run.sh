#!/bin/bash
# 
# Script to run Transmembrane test case
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
      -name 3LBW \
      -mtz ${PWD}/3LBW.sf.mtz   \
      -fasta ${PWD}/3LBW.fasta  \
      -rosetta_dir $ROSETTA_DIR \
      -transmembrane True \
      -nr "/opt/nr/nr" \
      -blast_dir "/opt/blast-2.2.26" \
      -usehoms False \
      -make_frags False \
      -make_models False \
      -shelx_cycles 10 \
      -nmodels 50 \
      -early_terminate True  \
      -percent 5 \
      -use_shelxe True \
      -use_arpwarp False \
      -use_buccaneer False \
      -molreponly True \
      -nproc 5 \
      -num_clusters 1 \
      -allatom True \
      -ensembles  ${PWD}/ROSETTA_MR_2/ensembles_1

#      -frags3mers ${PWD}/ROSETTA_MR_0/rosetta_fragments/3LBW_.200.3mers \
#     -frags9mers ${PWD}/ROSETTA_MR_0/rosetta_fragments/3LBW_.200.9mers \
#     -transmembrane_spanfile ${PWD}/ROSETTA_MR_0/rosetta_fragments/3LBW_.span\
#     -transmembrane_lipofile ${PWD}/ROSETTA_MR_0/rosetta_fragments/3LBW.lips4 \
