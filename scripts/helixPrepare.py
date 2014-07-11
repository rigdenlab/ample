#!/usr/bin/env ccp4-python

import cPickle
import os
import shutil
import subprocess
import sys

sys.path.insert(0, "/opt/ample-dev1/python")
sys.path.insert(0, "opt/ample-dev1/scripts")

import add_sidechains_SCWRL
import ample_util
from analyse_run import AmpleResult
import pdb_edit

#workdir = "/home/jmht/coiled-coils/ideal_helices"
workdir = "/home/jmht/Documents/work/CC/ideal_helices"
#smDir = "/home/jmht/coiled-coils/single_model"
ensembleDir = "/media/data/shared/coiled-coils/ensemble/ensemble.run2"
chainMakerExe="/home/jmht/Documents/avogadro/chainMaker/chainMaker"
scwrlExe="/opt/scwrl4/Scwrl4"
pfile = "/home/jmht/Documents/work/CC/ensemble_results/ar_results.pkl"


os.chdir(workdir)

# We've generated multiple sequences but many models will map onto each sequence
# For each model we are only interested if the sequence solves the structure - we can
# link back to the original ensemble later using a dictionary

# Maps the sequence to the list of ensembles 
sequence2ensemble = {}
pdbCode2Col = {}

# Unpickle the results
with open( pfile ) as f:
    allResults = cPickle.load( f )

scwrl = add_sidechains_SCWRL.Scwrl( scwrlExe=scwrlExe )
pdbEdit = pdb_edit.PDBEdit()

for result in allResults:

    # We only care about succeess
    if not ( result.shelxeCC >= 25 and result.shelxeAvgChainLength >= 10 ):
        continue

    # Similar to the single_models we only process the first radius ensembles
    if result.ensembleRadiusThreshold != 1:
        continue
    
    # See if we already have a dict for this target
    if result.pdbCode not in sequence2ensemble:
        sequence2ensemble[ result.pdbCode ] = {}
        
        # Need to extract the column labels for the mtz file - use the old script
        s = os.path.join( ensembleDir,
                          result.pdbCode,
                          result.ensembleName+".sub" 
                          )
        with open(s) as f:
            for line in f:
                if line.startswith("LABIN"):
                    #LABIN SIGF=SIGFP F=FP FreeR_flag=FREE
                    f = line.strip().split()
                    pdbCode2Col[ result.pdbCode ] = {}
                    l,v = f[1].split("=")
                    assert l =="SIGF"
                    pdbCode2Col[ result.pdbCode ]['SIGF'] = v
                    
                    l,v = f[2].split("=")
                    assert l =="F"
                    pdbCode2Col[ result.pdbCode ]['F'] = v
                    
                    l,v = f[3].split("=")
                    assert l =="FreeR_flag"
                    pdbCode2Col[ result.pdbCode ]['FREE'] = v
                    break
    
    # We need a helix
    if result.rioHelixSequence is None or len( result.rioHelixSequence ) < 3:
        continue
    
    # Got sequence so see if we already have it for this target
    if result.rioHelixSequence in sequence2ensemble[ result.pdbCode ]:
        # Add it to the list
        sequence2ensemble[ result.pdbCode ][ result.rioHelixSequence ].append( result.ensembleName )
    else:
        # New sequence so create a fresh list
        sequence2ensemble[ result.pdbCode ][ result.rioHelixSequence ] = [ result.ensembleName ]
        
    scriptDir = os.path.join( workdir, result.pdbCode )
    if not os.path.isdir( scriptDir ):
        os.mkdir( scriptDir )
    os.chdir( scriptDir )
    
    # Copy files in
    ofasta = os.path.join( ensembleDir, result.pdbCode, "{0}_1.fasta".format( result.pdbCode ) )
    fasta = os.path.join( scriptDir, "{0}_1.fasta".format( result.pdbCode ) )
    shutil.copy2( ofasta, fasta )
    
    omtz = os.path.join( ensembleDir, result.pdbCode, "{0}-cad.mtz".format( result.pdbCode ) )
    mtz = os.path.join( scriptDir, "{0}-cad.mtz".format( result.pdbCode ) )
    shutil.copy2( omtz, mtz )
    
    print "processing {0}-{1} in dir {2}".format( result.pdbCode, result.rioHelixSequence, scriptDir )
    
    # Generate full pdb from helix sequence
    #HACK REPLACE ANY X with A
    sequence = result.rioHelixSequence.replace('X','A')
    helixPdb1   = os.path.join( scriptDir, result.rioHelixSequence + ".pdb" )
    cmd = [ chainMakerExe, "-s", sequence, "-o", helixPdb1 ]
    #print "RUNNING "," ".join(cmd)
    retcode = ample_util.run_command(cmd, dolog=False)
    if not os.path.isfile( helixPdb1 ):
        raise RuntimeError,"ERROR creating PDB: {0}".format( helixPdb1 )
    
    # Optimise sidechain configuration with Scwrl
    helixPdbAA   = os.path.join( scriptDir, result.rioHelixSequence + "_All_atom.pdb" )
    scwrl.addSidechains( pdbin=helixPdb1, pdbout= helixPdbAA )
    
    helixPdbScwrl   = os.path.join( scriptDir, result.rioHelixSequence + "_SCWRL_reliable_sidechains.pdb" )
    helixPdbPolya   = os.path.join( scriptDir, result.rioHelixSequence + "_poly_ala.pdb" )
    
    # Now strip back as we did for the full ensembles
    pdbEdit.reliable_sidechains( helixPdbAA, helixPdbScwrl )
    pdbEdit.backbone( helixPdbAA, helixPdbPolya )
    
    #
    # Create scripts to run them all
    #
    for t in [ 'All_atom', 'SCWRL_reliable_sidechains', 'poly_ala']:
    
        newScript = os.path.join( scriptDir, result.rioHelixSequence+"_{0}.sub".format( t ) )
        stext = """#!/bin/bash

#$ -j y
#$ -cwd
#$ -w e
#$ -V
#$ -o /data1/jmht/ideal_helices/{0}/{1}_{2}.log
#$ -N J_{0}
#$ -pe threaded 16

pushd /data1/jmht/ideal_helices/{0}

setenv CCP4_SCR $TMPDIR


mrbump HKLIN /data1/jmht/ideal_helices/{0}/{0}-cad.mtz SEQIN /data1/jmht/ideal_helices/{0}/{0}_1.fasta HKLOUT OUT.mtz  XYZOUT OUT.pdb << eof
LABIN SIGF={3} F={4} FreeR_flag={5}
JOBID {0}_{1}_{2}
MRPROGRAM phaser
LOCALFILE /data1/jmht/ideal_helices/{0}/{1}_{2}.pdb CHAIN ALL RMS 0.1
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
SHLXEXE  /home/jmht/bin/shelxe
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

""".format( result.pdbCode,
                result.rioHelixSequence,
                t,
                pdbCode2Col[ result.pdbCode ]['SIGF'],
                pdbCode2Col[ result.pdbCode ]['F'],
                pdbCode2Col[ result.pdbCode ]['FREE'],
               )
    
        with open( newScript, 'w') as f:
            f.write(stext )

    #break

# Now write out mapping file
os.chdir(workdir)
with open("ensemble2Helix.txt", 'w') as f:
    for pdbCode in sequence2ensemble.keys():
        print sequence2ensemble[ pdbCode ]
        for sequence, elist in sequence2ensemble[ pdbCode ].iteritems():
            for ensemble in elist:
                f.write( "{0}    {1}    {2}\n".format( pdbCode, ensemble, sequence  ) 
                        )


with open("ensemble2Helix.pkl", 'w') as f:
    cPickle.dump( sequence2ensemble, f)
    
