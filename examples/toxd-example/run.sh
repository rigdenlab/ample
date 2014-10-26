#!/bin/bash

# 
# Script to run toxd test case
#

# Always need to supply the path to the Rosetta directory
ROSETTA_DIR=/opt/rosetta-3.5

# Change the below if the path to these external dependencies of AMPLE
# cannot be found in your standard PATH
export PATH=\
/opt/shelx:\
$PATH

# Set which ample binary we are using
AMPLE=$CCP4/bin/ample
#AMPLE=/opt/ample-dev1/bin/ample

$AMPLE \
-rosetta_dir $ROSETTA_DIR \
-fasta toxd_.fasta \
-mtz 1dtx.mtz \
-frags_3mers aat000_03_05.200_v1_3 \
-frags_9mers aat000_09_05.200_v1_3 \
-nmodels 30 \
-percent 50 \
-use_shelxe True \
-use_arpwarp False \
-use_buccaneer False \
-molrep_only True \
-nproc 2 \

#-models_dir ${AMPLEDIR}/examples/toxd-example/ROSETTA_MR_0/models \
#-ensembles_dir ${AMPLEDIR}/examples/toxd-example/ROSETTA_MR_0/ensembles_1 \
