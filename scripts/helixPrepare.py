#!/usr/bin/env python


"""


"""


import cPickle
import os
import shutil
import subprocess
import sys

sys.path.append("/home/jmht/ample-dev1/python")
sys.path.append("/home/jmht/ample-dev1/scripts")

import add_sidechains_SCWRL
import ample_util
from analyse_run import AmpleResult
import pdb_edit

workdir = "/home/jmht/coiled-coils/ideal_helices"
smDir = "/home/jmht/coiled-coils/single_model"
chainMakerExe="/home/jmht/chainMaker/chainMaker"
scwrlExe="/home/jmht/scwrl/Scwrl4"
pfile = "/home/jmht/analysis/full_ensemble/ar_results.pkl"

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
    
    # We need a helix
    if result.helixSequence is None or len( result.helixSequence ) < 3:
        continue
    
    #runDir  = os.path.join( scriptDir1, "search_" + result.ensembleName+"_mrbump" )
    #if os.path.isdir( runDir ):
    #    if os.path.isfile( os.path.join( runDir, "results", "finished.txt" )):
    #        #print "FINISHED ",runDir
    #        continue
    #    else:
    #        shutil.rmtree( runDir )
    
    scriptDir = os.path.join( workdir, result.pdbCode )
    if not os.path.isdir( scriptDir ):
        os.mkdir( scriptDir )
    os.chdir( scriptDir )
    newScript = os.path.join( scriptDir, result.ensembleName+".sub" )

    print "processing {0} in dir {1}".format( result.ensembleName, scriptDir )
    
    # Generate full pdb from helix sequence
    #HACK REPLACE ANY X with A
    sequence = result.helixSequence.replace('X','A')
    helixPdb   = os.path.join( scriptDir, result.ensembleName + ".pdb" )
    helixPdbCM = ample_util.filename_append( helixPdb, "chainMaker")
    cmd = [ chainMakerExe, "-s", sequence, "-o", helixPdbCM ]
    #print "RUNNING "," ".join(cmd)
    retcode = ample_util.run_command(cmd, dolog=False)
    if not os.path.isfile( helixPdbCM ):
        raise RuntimeError,"ERROR creating PDB: {0}".format( helixPdbCM )
    
    # Optimise sidechain configuration with Scwrl
    helixPdbScwrl = ample_util.filename_append( helixPdbCM, "scwrl")
    scwrl.addSidechains( pdbin=helixPdbCM, pdbout= helixPdbScwrl )
    
    # Now strip back as we did for the full ensembles
    if result.ensembleSideChainTreatment == "All_atom":
        # Copy scrwl file to be named ensemble
        shutil.copy2( helixPdbScwrl, helixPdb )
    elif result.ensembleSideChainTreatment == "SCWRL_reliable_sidechains":
        pdbEdit.reliable_sidechains( helixPdbScwrl, helixPdb )
    elif result.ensembleSideChainTreatment == "poly_ala":
        pdbEdit.backbone( helixPdbScwrl, helixPdb )
    else:
        assert False
    
    # Generate a fasta file for the helix
    #fasta = os.path.join( scriptDir, result.ensembleName+".fasta" )
    #with open( fasta, 'w') as w:
    #    w.write( ">helix for:{0}\n".format( result.ensembleName ) )
    #    # Max length is 62 so only one line
    #    w.write( result.helixSequence +"\n" )
    
    # Now copy the original script and update it with the new pdb
    oldScript = os.path.join( smDir, result.pdbCode, result.ensembleName+".sub" )
    newScript = os.path.join( scriptDir, result.ensembleName+".sub" )
    fasta = os.path.join( smDir, result.pdbCode, "{0}_1.fasta".format( result.pdbCode ) )
    mtz = os.path.join( smDir, result.pdbCode, "{0}-cad.mtz".format( result.pdbCode ) )
    
    with open( oldScript, 'r') as o, open( newScript, 'w' ) as n:
        for line in o:
            if line.startswith("#$ -o"):
                line = "#$ -o {0}\n".format( os.path.join(scriptDir,result.ensembleName+".log" ) )
            elif line.startswith("pushd"):
                line = "pushd {0}\n".format( scriptDir )
            elif line.startswith("mrbump"):
                line = "mrbump HKLIN {0} SEQIN {1} HKLOUT OUT.mtz  XYZOUT OUT.pdb << eof\n".format( mtz, fasta )
            elif line.startswith("LOCALFILE"):
                line = "LOCALFILE {0} CHAIN ALL RMS 0.1\n".format( helixPdb )
            
            # Write out the line to the new script
            n.write( line )
    
    # Make the script executble
    os.chmod(newScript, 0o777)
    
    # Submit the job
    subprocess.call( ["qsub",newScript] )


