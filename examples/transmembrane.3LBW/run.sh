#!/bin/bash
# 
# Script to run Transmembrane test case
#
# Change the below if the path to these external dependencies of AMPLE
# cannot be found in your standard PATH
ROSETTA_DIR=/opt/rosetta3.4
ROSETTA_DIR=/Users/jmht/Documents/AMPLE/programs/rosetta-3.5

#export PATH=\
#/opt/spicker:\
#/opt/theseus_src:\
#/opt/scwrl4:\
#/opt/maxcluster:\
#/opt/shelx:\
#$PATH

. /Users/jmht/Documents/AMPLE/programs/PATHS
BLAST_DIR=/Applications/Blast/blast-2.2.26
NR_DIR=/opt/nr/nr
NR_DIR=/Users/jmht/Documents/AMPLE/ample-dev1/examples/transmembrane.3LBW/nr

# This is only required if you are running ample from outside
# the standard CCP4 directory
AMPLEDIR=/opt/ample-dev1
AMPLEDIR=/Users/jmht/Documents/AMPLE/ample-dev1
export PYTHONPATH=${AMPLEDIR}/python:$PYTHONPATH

${AMPLEDIR}/bin/ample \
      -name 3LBW \
      -mtz 3LBW.sf.mtz   \
      -fasta 3LBW.fasta  \
      -rosetta_dir $ROSETTA_DIR \
      -transmembrane True \
      -nr $NR_DIR \
      -blast_dir $BLAST_DIR \
      -use_homs False \
      -frags_3mers 3LBW.200.3mers  \
      -frags_9mers 3LBW.200.9mers  \
      -make_models True \
      -shelx_cycles 10 \
      -nmodels 50 \
      -early_terminate True  \
      -percent 5 \
      -use_shelxe True \
      -use_arpwarp False \
      -use_buccaneer False \
      -molrep_only True \
      -nproc 2 \

#      -frags3mers ${PWD}/ROSETTA_MR_0/rosetta_fragments/3LBW_.200.3mers \
#      -num_clusters 1 \
#      -ensembles  ${PWD}/ROSETTA_MR_2/ensembles_1
#     -frags9mers ${PWD}/ROSETTA_MR_0/rosetta_fragments/3LBW_.200.9mers \
#     -transmembrane_spanfile ${PWD}/ROSETTA_MR_0/rosetta_fragments/3LBW_.span\
#     -transmembrane_lipofile ${PWD}/ROSETTA_MR_0/rosetta_fragments/3LBW.lips4 \
