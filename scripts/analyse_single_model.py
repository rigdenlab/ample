import copy
import csv
import cPickle
import os
import shutil
import sys
import traceback

from analyse_run import AmpleResult, MrbumpLogParser, MolrepLogParser, processMrbump, analyseSolution
from mrbump_results import MrBumpResult
import ample_util
import contacts
import csymmatch
import maxcluster
import pdb_edit
import pdb_model
import phaser_parser
import reforigin
import shelxe_log
import residue_map

def clearResult( ampleResult ):
    """Reset all things that may change"""
    
    ampleResult.solution = None
    ampleResult.resultDir = None
    ampleResult.rfact =  None
    ampleResult.rfree =  None
    ampleResult.mrProgram =  None
    
    ampleResult.phaserLLG = None
    ampleResult.phaserTFZ = None
    ampleResult.phaserPdb = None
    ampleResult.phaserTime = None
    
    ampleResult.molrepScore = None
    ampleResult.molrepTime = None
    ampleResult.molrepPdb = None
    
    ampleResult.shelxePdb = None
    ampleResult.shelxeCC = None
    ampleResult.shelxeAvgChainLength = None
    ampleResult.shelxeMaxChainLength = None
    ampleResult.shelxeNumChains = None
    
    
    ampleResult.reforiginRMSD = None
    ampleResult.floatingOrigin = None
    ampleResult.csymmatchOriginOk = None
    ampleResult.csymmatchOriginOk   = None
    ampleResult.shelxeCsymmatchShelxeScore  = None
    ampleResult.shelxeTM = None
    ampleResult.shelxeTMPairs = None
    ampleResult.shelxeRMSD =None
    
    ampleResult.contactData        = None
    ampleResult.numContacts        = None
    ampleResult.inregisterContacts = None
    ampleResult.ooregisterContacts = None
    ampleResult.backwardsContacts  = None
    ampleResult.contactOrigin      = None
    ampleResult.goodContacts       = None
    ampleResult.nocatContacts      = None
    ampleResult.helixSequence      = None
    ampleResult.lenHelix = None

    return




#
# MAIN STARTS HERE
#

rundir = "/home/jmht/Documents/test/CC/single_model"
singleModelDir="/media/data/shared/coiled-coils/single_model"
dataRoot="/media/data/shared/coiled-coils/ensemble/"

results = []

# Unpickle the original results
pfile = "/home/jmht/Documents/test/CC/run2/ar_results.pkl"
with open( pfile ) as f:
    ensembleResults = cPickle.load( f )

# Unpickle the mrbump results for this job
pfile = os.path.join( singleModelDir,"results.pkl" )
with open( pfile ) as f:
    mrbumpResults = cPickle.load( f )


pdbedit = pdb_edit.PDBEdit()
for pdbCode in sorted( mrbumpResults.keys() ):
    #if pdbCode != "1M5I":
    #    continue
    
    workdir = os.path.join( rundir, pdbCode )
    if not os.path.isdir( workdir ):
        os.mkdir( workdir )
    os.chdir( workdir )
        
    print "\nResults for ",pdbCode
    
    # Directory where all the data for this run live
    dataDir = os.path.join( dataRoot, pdbCode )
    
    # Get path to native Extract all the nativePdbInfo from it
    nativePdb = os.path.join( dataDir, "{0}.pdb".format( pdbCode ) )
    pdbedit = pdb_edit.PDBEdit()
    nativePdbInfo = pdbedit.get_info( nativePdb )
    
    # First check if the native has > 1 model and extract the first if so
    if len( nativePdbInfo.models ) > 1:
        print "nativePdb has > 1 model - using first"
        nativePdb1 = ample_util.filename_append( filename=nativePdb, astr="model1", directory=workdir )
        pdbedit.extract_model( nativePdb, nativePdb1, modelID=nativePdbInfo.models[0].serial )
        nativePdb = nativePdb1
        
    # Standardise the PDB to rename any non-standard AA, remove solvent etc
    nativePdbStd = ample_util.filename_append( filename=nativePdb, astr="std", directory=workdir )
    pdbedit.standardise( nativePdb, nativePdbStd )
    nativePdb = nativePdbStd
    
    # Get the new Info about the native
    nativePdbInfo = pdbedit.get_info( nativePdb )
    
    # Get information on the origins for this spaceGroup
    originInfo = pdb_model.OriginInfo( spaceGroupLabel=nativePdbInfo.crystalInfo.spaceGroup )
    
    # For maxcluster comparsion of shelxe model we need a single chain from the native so we get this here
    if len( nativePdbInfo.models[0].chains ) > 1:
        chainID = nativePdbInfo.models[0].chains[0]
        nativeAs1Chain  = ample_util.filename_append( filename=nativePdbInfo.pdb,
                                                       astr="1chain".format( chainID ), 
                                                       directory=workdir )
        pdbedit.to_single_chain( nativePdbInfo.pdb, nativeAs1Chain )
    else:
        nativeAs1Chain = nativePdbInfo.pdb
    
    
    # Get hold of a full model so we can do the mapping of residues
    refModelPdb = os.path.join( dataDir, "models/S_00000001.pdb".format( pdbCode ) )
    resSeqMap = residue_map.residueSequenceMap()
    refModelPdbInfo = pdbedit.get_info( refModelPdb )
    resSeqMap.fromInfo( refInfo=refModelPdbInfo,
                        refChainID=refModelPdbInfo.models[0].chains[0], # Only 1 chain in model
                        targetInfo=nativePdbInfo,
                        targetChainID=nativePdbInfo.models[0].chains[0]
                      )
    
    dsspLog = os.path.join( dataDir, "{0}.dssp".format( pdbCode  )  )

    
    # Loop through results for that job
    for mrbumpResult in mrbumpResults[ pdbCode ]:

        print "processing: {0}".format( mrbumpResult.name[9:-6] )
        #if mrbumpResult.name[9:-6] != "SCWRL_reliable_sidechains_trunc_24.348237_rad_1":
        #    continue
        
        # Need to specify the log
        mrbumpResult.mrbumpLog = os.path.join( singleModelDir, pdbCode, 
                                               "{0}.log".format( mrbumpResult.name[9:-6] ) )
    
        # Update the Mrbump result object and set all values in the Ample Result
        processMrbump( mrbumpResult )

        # Find the original result object
        ampleResult=None
        for e in ensembleResults:
            if e.pdbCode == pdbCode and e.ensembleName == mrbumpResult.ensembleName:
                ampleResult = copy.deepcopy( e )
        
        assert ampleResult is not None,"Could not find result: {0} {1}".format( pdbCode, mrbumpResult.ensembleName )
        clearResult( ampleResult )
        results.append( ampleResult )
        
        # Added here
        ampleResult.spaceGroup = originInfo.spaceGroup()
        
        ampleResult.solution =  mrbumpResult.solution
        ampleResult.resultDir = mrbumpResult.mrDir
        ampleResult.rfact =  mrbumpResult.rfact
        ampleResult.rfree =  mrbumpResult.rfree
        ampleResult.mrProgram =  mrbumpResult.program
    
        if mrbumpResult.program == "phaser":
            #ampleResult.phaserLog = mrbumpResult.phaserLog
            ampleResult.phaserLLG = mrbumpResult.phaserLLG
            ampleResult.phaserTFZ = mrbumpResult.phaserTFZ
            ampleResult.phaserPdb = mrbumpResult.phaserPdb
            ampleResult.phaserTime = mrbumpResult.phaserTime
        elif mrbumpResult.program == "molrep":
            #ampleResult.molrepLog = mrbumpResult.molrepLog
            ampleResult.molrepScore = mrbumpResult.molrepScore
            ampleResult.molrepTime = mrbumpResult.molrepTime
            ampleResult.molrepPdb = mrbumpResult.molrepPdb
        else:
            raise RuntimeError,"Unrecognised program!"
        
        #ampleResult.shelxeLog = mrbumpResult.shelxeLog
        ampleResult.shelxePdb = mrbumpResult.shelxePdb
        ampleResult.shelxeCC = mrbumpResult.shelxeCC
        ampleResult.shelxeAvgChainLength = mrbumpResult.shelxeAvgChainLength
        ampleResult.shelxeMaxChainLength = mrbumpResult.shelxeMaxChainLength
        ampleResult.shelxeNumChains = mrbumpResult.shelxeNumChains
        ampleResult.estChainsASU = mrbumpResult.estChainsASU
        
        try:
            analyseSolution( ampleResult=ampleResult,
                             nativePdbInfo=nativePdbInfo,
                             nativePdbAs1Chain=nativePdbAs1Chain,
                             refModelPdbInfo=refModelPdbInfo,
                             resSeqMap=resSeqMap,
                             originInfo=originInfo,
                             dsspLog=dsspLog,
                             workdir=workdir )
        except Exception,e:
            print "ERROR ANALYSING SOLUTION: {0} {1}".format( pdbCode, mrbumpResult.ensembleName )
            print traceback.format_exc()

for r in results:
    print r

pfile = os.path.join( rundir, "ar_results.pkl")
f = open( pfile, 'w' )
ampleDict = cPickle.dump( results, f  )

cpath = os.path.join( rundir, 'results.csv' )
csvfile =  open( cpath, 'wb')
csvwriter = csv.writer(csvfile, delimiter=',',
                        quotechar='"', quoting=csv.QUOTE_MINIMAL)

header=False
for r in results:
    if not header:
        csvwriter.writerow( r.titlesAsList() )
        header=True
    csvwriter.writerow( r.valuesAsList() )

csvfile.close()
