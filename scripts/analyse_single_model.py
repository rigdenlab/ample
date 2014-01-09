import copy
import csv
import cPickle
import os
import shutil
import sys
import traceback

from analyse_run import AmpleResult, MrbumpLogParser, MolrepLogParser
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


def processMrbump( mrbumpResult ):
    
    # Add attributes to object
    mrbumpResult.phaserLLG = None
    mrbumpResult.phaserTFZ = None
    mrbumpResult.phaserPdb = None
    mrbumpResult.phaserLog = None
    mrbumpResult.phaserTime = None
    mrbumpResult.molrepLog = None
    mrbumpResult.molrepScore = None
    mrbumpResult.molrepTime = None
    mrbumpResult.molrepPdb = None
    mrbumpResult.shelxePdb = None
    mrbumpResult.shelxeLog = None
    mrbumpResult.shelxeCC = None
    mrbumpResult.shelxeAvgChainLength = None
    mrbumpResult.shelxeMaxChainLength = None
    mrbumpResult.shelxeNumChains = None

    # Need to remove last component as we recored the refmac directory
    mrDir = os.sep.join( mrbumpResult.resultDir.split(os.sep)[:-1] )
    # HACK - we run the processing on cytosine so differnt place
    #/data2/jmht/coiled-coils/single_ensemble
    mrDir = mrDir.replace( "/data2/jmht/coiled-coils/single_ensemble","/media/data/shared/coiled-coils/single_model" )
    mrbumpResult.mrDir = mrDir
    
    mrbumpResult.ensembleName = mrbumpResult.name[9:-6]
    
    mrbumpP = MrbumpLogParser( mrbumpResult.mrbumpLog )
    mrbumpResult.estChainsASU = mrbumpP.noChainsTarget
    
    if mrbumpResult.program == "phaser":
        
        phaserPdb = os.path.join( mrDir,"refine","{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(mrbumpResult.program, mrbumpResult.ensembleName) )
        if os.path.isfile( phaserPdb ):
            phaserP = phaser_parser.PhaserPdbParser( phaserPdb )
            mrbumpResult.phaserLLG = phaserP.LLG
            mrbumpResult.phaserTFZ = phaserP.TFZ
            mrbumpResult.phaserPdb = phaserPdb
            
            phaserLog = os.path.join( mrDir, "{0}_loc0_ALL_{1}_UNMOD.log".format(mrbumpResult.program, mrbumpResult.ensembleName) )
            mrbumpResult.phaserLog = phaserLog
            phaserP = phaser_parser.PhaserLogParser( phaserLog )
            mrbumpResult.phaserTime = phaserP.time
        
    elif mrbumpResult.program == "molrep":
        molrepLog = os.path.join( mrDir, "molrep.log" )
        mrbumpResult.molrepLog = molrepLog
        molrepP = MolrepLogParser( molrepLog )
        mrbumpResult.molrepScore = molrepP.score
        mrbumpResult.molrepTime = molrepP.time
        
        molrepPdb = os.path.join( mrDir,"refine","{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(mrbumpResult.program, mrbumpResult.ensembleName) )
        if os.path.isfile( molrepPdb ):
            mrbumpResult.molrepPdb = molrepPdb
    else:
        assert False
        
    #
    # SHELXE PROCESSING
    #
    # Now read the shelxe log to see how we did
    shelxePdb = os.path.join( mrDir, "build/shelxe", "shelxe_{0}_loc0_ALL_{1}_UNMOD.pdb".format( mrbumpResult.program, mrbumpResult.ensembleName ) )
    if os.path.isfile( shelxePdb):
        mrbumpResult.shelxePdb = shelxePdb
        
    shelxeLog = os.path.join( mrDir, "build/shelxe/shelxe_run.log" )
    if os.path.isfile( shelxeLog ):
        mrbumpResult.shelxeLog = shelxeLog
        shelxeP = shelxe_log.ShelxeLogParser( shelxeLog )
        mrbumpResult.shelxeCC = shelxeP.CC
        mrbumpResult.shelxeAvgChainLength = shelxeP.avgChainLength
        mrbumpResult.shelxeMaxChainLength = shelxeP.maxChainLength
        mrbumpResult.shelxeNumChains= shelxeP.numChains
    
    return

def analyseSolution( ampleResult=None,
                     nativePdbInfo=None,
                     nativePdbSingle=None,
                     refModelPdbInfo=None,
                     resSeqMap=None,
                     originInfo=None,
                     dsspLog=None,
                     workdir=None ):


    if ampleResult.mrProgram == "phaser":
        placedPdb = ampleResult.phaserPdb
    elif ampleResult.mrProgram == "molrep":
        placedPdb = ampleResult.molrepPdb
    else:
        assert False


    if placedPdb is None:
     print "NO PDB FOR ",ampleResult
     return
    #else:
    # print "GOT PDB FOR ",ampleResult

    # debug - copy into work directory as reforigin struggles with long pathnames
    shutil.copy(placedPdb, os.path.join( workdir, os.path.basename( placedPdb ) ) )
    
    placedPdbInfo = pdbedit.get_info( placedPdb )
    
    # Get reforigin info
    if True:
    #try:
        rmsder = reforigin.ReforiginRmsd()
        rmsder.getRmsd(  nativePdbInfo=nativePdbInfo,
                         placedPdbInfo=placedPdbInfo,
                         refModelPdbInfo=refModelPdbInfo,
                         cAlphaOnly=True )
        ampleResult.reforiginRMSD = rmsder.rmsd
    #except Exception, e:
    #    print "ERROR: ReforiginRmsd with: {0} {1}".format( nativePdbInfo.pdb, placedPdbInfo.pdb )
    #    print "{0}".format( e )
    #    ampleResult.reforiginRMSD = 9999
         
    #
    # SHELXE PROCESSING
    #
    if not ampleResult.shelxePdb is None and os.path.isfile( ampleResult.shelxePdb ):
        # Need to copy to avoid problems with long path names
        shelxePdb = os.path.join(workdir, os.path.basename( ampleResult.shelxePdb ) )
        shutil.copy( ampleResult.shelxePdb, shelxePdb )
        
        csym                           = csymmatch.Csymmatch()
        shelxeCsymmatchPdb             = ample_util.filename_append( 
                                                                    filename=shelxePdb, 
                                                                    astr="csymmatch", 
                                                                    directory=workdir )
        
        # Had problem with shelx losing origin information
        csym.run( refPdb=nativePdbInfo.pdb, inPdb=shelxePdb, outPdb=shelxeCsymmatchPdb )
        shelxeCsymmatchOrigin = csym.origin()
        
        # See if this origin is valid
        ampleResult.floatingOrigin = originInfo.isFloating()
        ampleResult.csymmatchOriginOk = True
        if not shelxeCsymmatchOrigin or \
        ( shelxeCsymmatchOrigin not in originInfo.redundantAlternateOrigins() and not ampleResult.floatingOrigin ):
            ampleResult.csymmatchOriginOk   = False
            shelxeCsymmatchOrigin  = None
        
        ampleResult.shelxeCsymmatchShelxeScore  = csym.averageScore()
        shelxeCsymmatchPdbSingle       = ample_util.filename_append( filename=shelxeCsymmatchPdb, 
                                                                     astr="1chain", 
                                                                     directory=workdir )
        pdbedit.to_single_chain(shelxeCsymmatchPdb, shelxeCsymmatchPdbSingle)
        
        # Compare the traced model to the native with maxcluster
        # We can only compare one chain so we extracted this earlier
        maxComp = maxcluster.Maxcluster()
        d = maxComp.compareSingle( nativePdb=nativePdbSingle,
                                   modelPdb=shelxeCsymmatchPdbSingle,
                                   sequenceIndependant=True,
                                   rmsd=False
                                 )
        ampleResult.shelxeTM = d.tm
        ampleResult.shelxeTMPairs = d.pairs
        
        d = maxComp.compareSingle( nativePdb=nativePdbSingle,
                                   modelPdb=shelxeCsymmatchPdbSingle,
                                   sequenceIndependant=True,
                                   rmsd=True )
        ampleResult.shelxeRMSD = d.rmsd

        # Now calculate contacts
    
        # Only bother when we have a floating origin if the csymmatch origin is ok
        if not ampleResult.floatingOrigin or ( ampleResult.floatingOrigin and ampleResult.csymmatchOriginOk ):
            ccalc = contacts.Contacts()
            #try:
            if True:
                ccalc.getContacts( placedPdbInfo=placedPdbInfo,
                                   nativePdbInfo=nativePdbInfo,
                                   resSeqMap=resSeqMap,
                                   originInfo=originInfo,
                                   shelxeCsymmatchOrigin=shelxeCsymmatchOrigin,
                                   workdir=workdir,
                                   dsspLog=dsspLog
                                )
            #except Exception, e:
            #    print "ERROR WITH CONTACTS: {0}".format( e )
       
            if ccalc.best:
                ampleResult.contactData        = ccalc.best
                ampleResult.numContacts        = ccalc.best.numContacts
                ampleResult.inregisterContacts = ccalc.best.inregister
                ampleResult.ooregisterContacts = ccalc.best.ooregister
                ampleResult.backwardsContacts  = ccalc.best.backwards
                ampleResult.contactOrigin      = ccalc.best.origin
                ampleResult.goodContacts       = ampleResult.inregisterContacts + ampleResult.ooregisterContacts
                ampleResult.nocatContacts      = ampleResult.numContacts - ampleResult.goodContacts
                ampleResult.helixSequence      = ccalc.best.helix
                if ccalc.best.helix:
                    ampleResult.lenHelix = len( ccalc.best.helix )
                
                gotHelix=False
                hfile = os.path.join( workdir, "{0}.helix".format( ampleResult.ensembleName ) )
                gotHelix =  ccalc.writeHelixFile( hfile )
                        
                # Just for debugging
                if ampleResult.shelxeCC >= 25 and ampleResult.shelxeAvgChainLength >= 10 and not gotHelix:
                    print "NO HELIX FILE"

    return

#
# MAIN STARTS HERE
#

rundir = "/home/jmht/Documents/test/CC/single_model"
singleModelDir="/media/data/shared/coiled-coils/single_model"
dataDir="/media/data/shared/coiled-coils/ensemble/"

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
    if pdbCode != "1M5I":
        continue
    
    workdir = os.path.join( rundir, pdbCode )
    if not os.path.isdir( workdir ):
        os.mkdir( workdir )
    os.chdir( workdir )
        
    print "\nResults for ",pdbCode
    
    # Get path to native Extract all the nativeInfo from it
    nativePdb = os.path.join( dataDir, pdbCode,  "{0}.pdb".format( pdbCode ) )
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

    # For maxcluster comparsion of shelxe model we need a single chain from the native so we get this here
    if len( nativePdbInfo.models[0].chains ) > 1:
#         chainID = nativePdbInfo.models[0].chains[0]
#         nativePdbSingle  = ample_util.filename_append( filename=nativePdbInfo.pdb,
#                                                        astr="chain{0}".format( chainID ), 
#                                                        directory=workdir )
#         pdbedit.extract_chain( nativePdbInfo.pdb, nativePdbSingle, chainID=chainID )
        chainID = nativePdbInfo.models[0].chains[0]
        nativePdbSingle  = ample_util.filename_append( filename=nativePdbInfo.pdb,
                                                       astr="1chain".format( chainID ), 
                                                       directory=workdir )
        pdbedit.to_single_chain( nativePdbInfo.pdb, nativePdbSingle )
    else:
        nativePdbSingle = nativePdbInfo.pdb
    
    # Get information on the origins for this spaceGroup
    originInfo = pdb_model.OriginInfo( spaceGroupLabel=nativePdbInfo.crystalInfo.spaceGroup )
    dsspLog = os.path.join( dataDir, pdbCode,  "{0}.dssp".format( pdbCode  )  )

    # Get hold of a full model so we can do the mapping of residues
    refModelPdb = os.path.join( dataDir, pdbCode,  "models/S_00000001.pdb".format( pdbCode ) )
    refModelPdbInfo = pdbedit.get_info( refModelPdb )
    
    resSeqMap = residue_map.residueSequenceMap()
    resSeqMap.fromInfo( refInfo=refModelPdbInfo,
                        refChainID=refModelPdbInfo.models[0].chains[0],
                        targetInfo=nativePdbInfo,
                        targetChainID=nativePdbInfo.models[0].chains[0]
                       )
    
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
        
        try:
            analyseSolution( ampleResult=ampleResult,
                             nativePdbInfo=nativePdbInfo,
                             nativePdbSingle=nativePdbSingle,
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
