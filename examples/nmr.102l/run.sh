#!/bin/sh

# NMR ensembling example - 2LC9 is an ensemble model of a minor and transiently formed state of a T4 lysozyme mutant. The target 102l is X-ray data

ccp4-python -u ${CCP4}/bin/ample -mtz 102l.mtz -fasta 102L.fasta -NMR_model_in 2LC9.pdb \
      -name 102l \
      -NMR_Truncate_only True \
      -use_homs True \
      -early_terminate True  \
      -use_arpwarp False \
      -molrep_only True \
      -nproc 1 \
      -F F -SIGF SIGF -FREE FreeR_flag \
      -use_arpwarp False \
      -run_dir ${PWD}
