#!/usr/bin/env python


"""


"""


import cPickle
import os
import shutil
import subprocess

#sys.path.append("/opt/ample-dev1/python")
#sys.path.append("/opt/ample-dev1/scripts")
sys.path.append("/Users/jmht/Documents/AMPLE/ample-dev1/python")
sys.path.append("/Users/jmht/Documents/AMPLE/ample-dev1/scripts")

import add_sidechains_SCWRL
import ample_util
import pdb_edit

# Unpickle the results
pfile = os.path.join( dataRoot,"results.pkl" )
with open( pfile ) as f:
    allResults = cPickle.load( f )

workdir = "wdir"
smDir = ""
chainMakerExe=""
scwrlExe=""

scwrl = add_sidechains_SCWRL.Scwrl( scwrlExe=scwrlExe )
pdbEdit = pdb_edit.PDBEdit()

for result in allResults:
    
    # We only care about succeess
    if not ( result.shelxeCC >= 25 and result.shelxeAvgChainLength >= 10 ):
        continue
    
    # We need a helix
    if result.helixSequence is None or len( result.helixSequence ) < 1:
        continue
    
    scriptDir = os.path.join( workdir, result.pdbCode )
    if not os.path.isdir( scriptDir ):
        os.mkdir( scriptDir )
    
    os.chdir( scriptDir )
    
    # Generate full pdb from helix sequence
    helixPdbScwrl = os.path.join( workdir, result.ensembleName + ".pdb" )
    helixPdbCM = ample_util.filename_append( helixPdb, "chainMaker")
    cmd = [ chainMaker, "-s", result.helixSequence, "-o", helixPdbCM ]
    retcode = ample_util.run_command(cmd, dolog=False)
    
    # Optimise sidechain configuration with Scwrl
    helixPdbScwrl = ample_util.filename_append( helixPdbCM, "scwrl")
    scwrl.addSidechains( pdbin=helixPdbCM, pdbout= helixPdbScwrl )
    
    # Now strip back as we did for the full ensembles
    if result.ensembleSideChainTreatment == "all_atom":
        # Copy scrwl file to be named ensemble
        shutil.copy2( helixPdbScwrl, helixPdb )
    elif result.ensembleSideChainTreatment == "reliable_sidechains":
        pdbEdit.reliable_sidechains( helixPdbScwrl, helixPdb )
    elif result.ensembleSideChainTreatment == "poly_ala":
        pdbEdit.backbone( helixPdbScwrl, helixPdb )
    
    # Generate a fasta file for the helix
    fasta = os.path.join( scriptDir, result.ensembleName+".fasta" )
    with open( fasta, 'w') as w:
        w.write( ">helix for:{0}\n".format( result.ensembleName ) )
        # Max length is 62
        w.write( result.helixSequence +"\n" )
    
    # Now copy the original script and update it with the new pdb
    oldScript = os.path.join( smDir, result.pdbCode, result.ensembleName+".sub" )
    newScript = os.path.join( scriptDir, result.pdbCode, result.ensembleName+".sub" )
    mtz = os.path.join( smDir, result.pdbCode, "{0}-cad.mtz".format( result.pdbCode ) )
    
    with open( oldScript, 'r') as o, open( newScript, 'w' ) as n:
        for line in o:
            if line.startswith("#BSUB -o"):
                line = "#BSUB -o {0}\n".format( os.path.join(scriptDir,result.ensembleName+".log" ) )
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
    #subprocess.call( ["qsub",newScript] )
    
    