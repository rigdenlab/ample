#!/bin/bash

# 
# Script to run Transmembrane test case
#

# Set path to include where shelxe is located
export PATH=\
/opt/shelx:\
$PATH

# Need to specify where rosetta is located
rosetta_dir=/opt/rosetta_2014.35.57232_bundle

# If you are using rosetta <= 3.5 then you need to specify the
# location of the blast root directory and the nr database directory.
# (later versions of rosetta will install these for you)
# If so, add the following flags to the script
# -blast_dir /Applications/Blast/blast-2.2.26 \
# -nr /opt/nr/nr \

${CCP4}/bin/ample.py \
-name 3LBW \
-mtz 3LBW.sf.mtz   \
-fasta 3LBW.fasta  \
-rosetta_dir $rosetta_dir \
-transmembrane True \
-use_homs False \
-frags_3mers 3LBW.200.3mers  \
-frags_9mers 3LBW.200.9mers  \
-nproc 5 \

