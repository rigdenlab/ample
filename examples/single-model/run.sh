#!/bin/bash

$CCP4/bin/ample \
-fasta input/1ujb.fasta \
-mtz input/1ujb-sf.mtz \
-single_model input/3c7t.pdb \
-do_mr False \
-truncation_method scores \
-truncation_scorefile input/3c7t_scores.csv \
-truncation_scorefile_header residue Concoord AMN  \

