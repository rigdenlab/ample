#!/bin/sh

# NMR ensembling example - 2LC9 is an ensemble model of a minor and transiently formed state of a T4 lysozyme mutant. The target 102l is X-ray data



AMPLE=${CCP4}/bin/ample
AMPLE=/opt/ample-dev1/bin/ample.py

ccp4-python -u $AMPLE \
-mtz 102l.mtz \
-fasta 102L.fasta \
-name 102l \
-NMR_model_in 2LC9.pdb \
-NMR_Truncate_only True \
-quick_mode True
