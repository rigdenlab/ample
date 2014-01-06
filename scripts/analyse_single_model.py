import copy
import sys

from import analyse_run import AmpleResult
from mrbump_results import MrBumpResult

def processMrbump( mrbumpResult ):
    
    # Need to remove last component as we recored the refmac directory
    mrbumpResult.mrDir = os.sep.join( mrbumpResult.resultDir.split(os.sep)[:-1] )
    
    mrbumpResult.ensembleName = mrbumpResult.name[9:-6]
    
    mrbumpP = MrbumpLogParser( mrbumpResult.mrbumpLog )
    mrbumpResult.estChainsASU = mrbumpP.noChainsTarget
    
    if mrbumpResult.program == "phaser":
        
        phaserPdb = os.path.join( resultDir,"refine","{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(mrbumpResult.program, ensembleName) )
        if os.path.isfile( phaserPdb ):
            phaserP = phaser_parser.PhaserPdbParser( phaserPdb )
            mrbumpResult.phaserLLG = phaserP.LLG
            mrbumpResult.phaserTFZ = phaserP.TFZ
            mrbumpResult.phaserPdb = phaserPdb
            
            phaserLog = os.path.join( resultDir, "{0}_loc0_ALL_{1}_UNMOD.log".format(mrbumpResult.program, ensembleName) )
            mrbumpResult.phaserLog = phaserLog
            phaserP = phaser_parser.PhaserLogParser( phaserLog )
            mrbumpResult.phaserTime = phaserP.time
        
    else:
        molrepLog = os.path.join( resultDir, "molrep.log" )
        mrbumpResult.molrepLog = molrepLog
        molrepP = MolrepLogParser( molrepLog )
        mrbumpResult.molrepScore = molrepP.score
        mrbumpResult.molrepTime = molrepP.time
        
        molrepPdb = os.path.join( resultDir,"refine","{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(mrbumpResult.program, ensembleName) )
        if is os.path.isfile( molrepPdb ):
            mrbumpResult.molrepPdb = molrepPdb
        
    #
    # SHELXE PROCESSING
    #
    # Now read the shelxe log to see how we did
    shelxePdb = os.path.join( resultDir, "build/shelxe", "shelxe_{0}_loc0_ALL_{1}_UNMOD.pdb".format( mrbumpResult.program, ensembleName ) )
    if os.path.isfile( shelxePdb):
        mrbumpResult.shelxePdb
        
    shelxeLog = os.path.join( resultDir, "build/shelxe/shelxe_run.log" )
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
                     refModelPdb=None,
                     resSeqMap=None,
                     originInfo=None,
                     dsspLog=None,
                     workdir=None ):
    
    if ampleResult.program == "phaser":
        placedPdb = ampleResult.phaserPdb
    elif ampleResult.program == "molrep":
        placedPdb = ampleResult.molrepPdb
    
    # debug - copy into work directory as reforigin struggles with long pathnames
    shutil.copy(placedPdb, os.path.join( workdir, os.path.basename( placedPdb ) ) )
    
    placedInfo = pdbedit.get_info( placedPdb )
    
    # Get reforigin info
    try:
        rmsder = reforigin.ReforiginRmsd( nativePdbInfo=nativeInfo,
                                          placedPdbInfo=placedInfo,
                                          refModelPdb=refModelPdb,
                                          cAlphaOnly=True )
        ampleResult.reforiginRMSD = rmsder.rmsd
    except Exception, e:
        print "ERROR: ReforiginRmsd with: {0} {1}".format( nativePdbInfo.pdb, placedPdbInfo.pdb )
        print "{0}".format( e )
        ampleResult.reforiginRMSD = 9999
         
    #
    # SHELXE PROCESSING
    #
    if os.path.isfile( ampleResult.shelxePdb ):
        shelxePdb = os.path.join(workdir, os.path.basename( origShelxePdb ) )
        shutil.copy( ampleResult.shelxePdb, shelxePdb )
        
        csym                           = csymmatch.Csymmatch()
        shelxeCsymmatchPdb             = ample_util.filename_append( 
                                                                    filename=shelxePdb, 
                                                                    astr="csymmatch", 
                                                                    directory=workdir )
        csym.run( refPdb=nativePdb, inPdb=shelxePdb, outPdb=shelxeCsymmatchPdb )
        shelxeCsymmatchOrigin = csym.origin()
        
        # See if this origin is valid
        ampleResult.floatingOrigin = originInfo.isFloating()
        ampleResult.csymmatchOriginOk = True
        if not shelxeCsymmatchOrigin or \
        ( shelxeCsymmatchOrigin not in originInfo.redundantAlternateOrigins() and not ampleResult.floatingOrigin ):
            ampleResult.csymmatchOriginOk   = False
            shelxeCsymmatchOrigin  = None
        
        shelxeCsymmatchShelxeScore     = csym.averageScore()
        ampleResult.shelxeCsymmatchShelxeScore  = shelxeCsymmatchShelxeScore
        
        shelxeCsymmatchPdbSingle       = ample_util.filename_append( filename=shelxeCsymmatchPdb, 
                                                                     astr="1chain", 
                                                                     directory=workdir )
        pdbedit.to_single_chain(shelxeCsymmatchPdb, shelxeCsymmatchPdbSingle)
        
        # Compare the traced model to the native with maxcluster
        d = maxComp.compareSingle( nativePdb=nativePdbInfo.pdb,
                               modelPdb=shelxeCsymmatchPdbSingle,
                               sequenceIndependant=True,
                               rmsd=False )
        
        ampleResult.shelxeTM = d.tm
        ampleResult.shelxeTMPairs = d.pairs
        d = maxComp.compareSingle( nativePdb=nativePdbInfo.pdb,
                               modelPdb=shelxeCsymmatchPdbSingle,
                               sequenceIndependant=True,
                               rmsd=True )
        
        ampleResult.shelxeRMSD = d.rmsd
        

        # Now calculate contacts
    
        # Only bother when we have a floating origin if the csymmatch origin is ok
        if not ampleResult.floatingOrigin or ( ampleResult.floatingOrigin and ampleResult.csymmatchOriginOk ):
            ccalc = contacts.Contacts()
            try:
            #if True:
                ccalc.getContacts( placedPdb=placedPdb,
                                   resSeqMap=resSeqMap,
                                   nativeInfo=nativeInfo,
                                   originInfo=originInfo,
                                   shelxeCsymmatchOrigin=shelxeCsymmatchOrigin,
                                   workdir=workdir,
                                   dsspLog=dsspLog
                                )
            except Exception, e:
                print "ERROR WITH CONTACTS: {0}".format( e )
       
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
                hfile = os.path.join( workdir, "{0}.helix".format( ensembleName ) )
                gotHelix =  ccalc.writeHelixFile( hfile )
                        
                # Just for debugging
                if ampleResult.shelxeCC >= 25 and ampleResult.shelxeAvgChainLength >= 10 and not gotHelix:
                    print "NO HELIX FILE"

    return

#
# MAIN STARTS HERE
#

rundir = ""
singleModelDir=""
dataDir=""

# Unpickle the original results
pfile = os.path.join( dataRoot,"ar_results.pkl" )
with open( pfile ) as f:
    ensembleResults = cPickle.load( f )

# Unpickle the mrbump results for this job
pfile = os.path.join( singleModelDir,"results.pkl" )
with open( pfile ) as f:
    mrbumpResults = cPickle.load( f )

pdbedit = pdb_edit.PDBEdit()


for pdbcode in sorted( mrbumpResults.keys() ):
    
    workdir = os.path.join( rundir, pdbcode )
    if not os.path.isdir( workdir ):
        os.mkdir( workdir )
    os.chdir( workdir )
        
    print "\nResults for ",pdbcode
    
    # Get path to native Extract all the nativeInfo from it
    nativePdb = os.path.join( dataDir, "{0}.pdb".format( pdbcode ) )
    pdbedit = pdb_edit.PDBEdit()
    nativeInfo = pdbedit.get_info( nativePdb )
    
    # First check if the native has > 1 model and extract the first if so
    if len( nativeInfo.models ) > 1:
        print "nativePdb has > 1 model - using first"
        nativePdb1 = ample_util.filename_append( filename=nativePdb, astr="model1", directory=workdir )
        pdbedit.extract_model( nativePdb, nativePdb1, modelID=nativeInfo.models[0].serial )
        nativePdb = nativePdb1
        
    # Standardise the PDB to rename any non-standard AA, remove solvent etc
    nativePdbStd = ample_util.filename_append( filename=nativePdb, astr="std", directory=workdir )
    pdbedit.standardise( nativePdb, nativePdbStd )
    nativePdb = nativePdbStd
    
    # Get the new Info about the native
    nativeInfo = pdbedit.get_info( nativePdb )
    
    # Get information on the origins for this spaceGroup
    originInfo = pdb_model.OriginInfo( spaceGroupLabel=nativeInfo.crystalInfo.spaceGroup )
    dsspLog = os.path.join( dataDir, "{0}.dssp".format( pdbcode  )  )

    # Get hold of a full model so we can do the mapping of residues
    refModelPdb = os.path.join( dataDir, "models/S_00000001.pdb".format( pdbcode ) )
    resSeqMap = residue_map.residueSequenceMap()
    
    modelInfo = pdbedit.get_info( refModelPdb )
    
    resSeqMap.fromInfo( refInfo=modelInfo,
                            refChainID=modelInfo.models[0].chains[0],
                            targetInfo=nativeInfo,
                            targetChainID=nativeInfo.models[0].chains[0]
                            )
    
    # Loop through results for that job
    for mrbumpResult in mrbumpResults[ pdbcode ]:
        
        # Need to specify the log
        mrbumpLog = 
        mrbumpResult.mrbumpLog = os.path.join( dataDir, "ROSETTA_MR_0/MRBUMP/cluster_1/",
                                               "{0}.sub.log".format( mrbumpResult.name[9:-6] ) )
    
        # Update the Mrbump result object and set all values in the Ample Result
        processMrbump( mrbumpResult )
        
        ampleResult.solution =  mrbumpResult.solution
        ampleResult.resultDir = mrbumpResult.mrDir
        ampleResult.rfact =  mrbumpResult.rfact
        ampleResult.rfree =  mrbumpResult.rfree
        ampleResult.mrProgram =  mrbumpResult.program
    
        ampleResult.phaserLog = mrbumpResult.phaserLog
        ampleResult.phaserLLG = mrbumpResult.phaserLLG
        ampleResult.phaserTFZ = mrbumpResult.phaserTFZ
        ampleResult.phaserPdb = mrbumpResult.phaserPdb
        ampleResult.phaserTime = mrbumpResult.phaserTime
        
        ampleResult.molrepLog = mrbumpResult.molrepLog
        ampleResult.molrepScore = mrbumpResult.molrepScore
        ampleResult.molrepTime = mrbumpResult.molrepTime
        ampleResult.molrepPdb = mrbumpResult.molrepPdb
        
        ampleResult.shelxeLog = mrbumpResult.shelxeLog
        ampleResult.shelxeCC = mrbumpResult.shelxeCC
        ampleResult.shelxeAvgChainLength = mrbumpResult.shelxeAvgChainLength
        ampleResult.shelxeMaxChainLength = mrbumpResult.shelxeMaxChainLength
        ampleResult.shelxeNumChains= mrbumpResult.shelxeNumChains
        
        # Find the original result object
        ampleResult=None
        for e in ensembleResults:
            if e.pdbCode == pdbCode and e.ensembleName == mrbumpResult.ensembleName:
                ampleResult = copy.deepcopy( e )
        
        assert ampleResult is not None,"Could not find result: {0} {1}".format( pdbCode,ensembleName  )
        
        analyseSolution( ampleResult=ar,
                         nativePdbInfo=nativePdbInfo,
                         refModelPdb=refModelPdb,
                         resSeqMap=resSeqMap,
                         originInfo=originInfo,
                         dsspLog=dsspLog,
                         workdir=workdir )
