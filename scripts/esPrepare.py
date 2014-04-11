#!/usr/bin/env ccp4-python

import cPickle
import os
import shutil
import sys

sys.path.insert(0, "/opt/ample-dev1/python")
sys.path.insert(0, "opt/ample-dev1/scripts")

import ample_util
import split_models

ensembleDir = "/media/data/shared/coiled-coils/ensemble/ensemble.run1"
workdir = "/home/jmht/Documents/work/CC/ensemble_single_chadwick"
os.chdir(workdir)

sfail = "/home/jmht/Documents/work/CC/ensemble_results/single.fail"


d = [ tuple(line.strip().split()) for line in open(sfail) ]

for pdbCode, ensembleName in d:
    
    # Need to extract the column labels for the mtz file - use the old script
    s = os.path.join( ensembleDir,
                      pdbCode,
                      "ROSETTA_MR_0/MRBUMP/cluster_1",
                      ensembleName+".sub" 
                      )
    
    SIGF = None
    F    = None
    FREE = None
    with open(s) as f:
        for line in f:
            if line.startswith("LABIN"):
                #LABIN SIGF=SIGFP F=FP FreeR_flag=FREE
                f = line.strip().split()
                l,v = f[1].split("=")
                assert l =="SIGF"
                SIGF = v
                
                l,v = f[2].split("=")
                assert l =="F"
                F = v
                
                l,v = f[3].split("=")
                assert l =="FreeR_flag"
                FREE = v
                break
            
    print F,SIGF,FREE
    
    
    # Create the directory for the files
    jobDir = os.path.join( workdir,"{0}_{1}",format(pdbCode,ensembleName))
    os.mkdir(jobDir)
    sys.exit()
    
    # Copy files in
    ofasta = os.path.join( ensembleDir, pdbCode, "{0}_1.fasta".format( pdbCode ) )
    fasta = os.path.join( jobDir, "{0}_1.fasta".format( pdbCode ) )
    shutil.copy2( ofasta, fasta )
    
    omtz = os.path.join( ensembleDir, pdbCode, "{0}-cad.mtz".format( pdbCode ) )
    mtz = os.path.join( jobDir, "{0}-cad.mtz".format( pdbCode ) )
    shutil.copy2( omtz, mtz )
    
    print "processing {0}-{1} in dir {2}".format( pdbCode, ensembleName, jobDir )
    
    # Split the ensemble into models
    ensembleFile = os.path.join( ensembleDir,
                                 pdbCode,
                                 "ROSETTA_MR_0/ensembles_1",
                                 ensembleName+".pdb" 
                                 )
    modelIds = split_models.split_pdb( ensembleFile )
    
    #
    # Create scripts to run them all
    #
    for ensemblePdb in modelIds:
        
        ename = os.path.splitext( os.path.split( ensemblePdb )[1] )[0]
    
        newScript = os.path.join( jobDir, "{0}.sub".format( ename ) )
        stext = """#!/bin/bash


# Set up SGE variables
#$ -j y
#$ -cwd
#$ -w e
#$ -V
#$ -o /volatile/jmht42/coiled-coils/ensemble_single/{0}/{2}.log
#$ -pe smp 16

pushd /volatile/jmht42/coiled-coils/ensemble_single/{0}


mrbump HKLIN /volatile/jmht42/coiled-coils/ensemble_single/{0}/{0}-cad.mtz SEQIN /volatile/jmht42/coiled-coils/ensemble_single/{0}/{0}_1.fasta HKLOUT OUT.mtz  XYZOUT OUT.pdb << eof
LABIN SIGF={3} F={4} FreeR_flag={5}
JOBID {0}_{2}
MRPROGRAM phaser
LOCALFILE /volatile/jmht42/coiled-coils/ensemble_single/{0}/{2}.pdb CHAIN ALL RMS 0.1
SCOPSEARCH False
PQSSEARCH False
SSMSEARCH False
FAST False
DOFASTA False
MDLD False
MDLC False
MDLM False
MDLP False
MDLS False
MDLU True
UPDATE False
BUCC  False
BCYCLES  5
ARPWARP  False
ACYCLES  10
SHELXE  True
SHLXEXE  /home/jmht42/bin/shelxe
SCYCLES  15
FIXSG True
SXREBUILD TRUE
CHECK False
LITE True
PICKLE False
TRYALL True
USEACORN False
USEENSEM False
CLEAN False
DEBUG True
PJOBS  16
PKEY  KILL  TIME  360
END
eof

popd

""".format( pdbCode,
            ensembleName,
            ename,
            SIGF,
            F,
            FREE )
    
        with open( newScript, 'w') as f:
            f.write(stext )

    break

    
