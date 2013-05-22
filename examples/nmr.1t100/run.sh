#!/bin/bash
# 
# Script to run NMR test case
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
      -name 1t00 \
      -mtz ${PWD}/1t00.mtz \
      -fasta ${PWD}/1T00_.fasta  \
      -NMR_model_in ${PWD}/1OKD.pdb \
      -frags3mers ${PWD}/1t00_.200.3mers \
      -frags9mers ${PWD}/1t00_.200.9mers \
      -rosetta_dir $ROSETTA_DIR \
      -usehoms True \
      -make_frags False \
      -make_models True \
      -shelx_cycles 10 \
      -nmodels 100 \
      -early_terminate True  \
      -percent 5 \
      -use_shelxe True \
      -use_arpwarp False \
      -use_buccaneer False \
      -molreponly True \
      -nproc 5 \
      -num_clusters 1 \
      -allatom False \
      -F FP \
      -SIGF SIGFP \
      -FREE FREE  \

