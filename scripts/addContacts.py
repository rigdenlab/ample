import sys
sys.path.append( "/opt/ample-dev1/python" )
sys.path.append( "/opt/ample-dev1/scripts" )

import cPickle
import csv
import os
import shutil

import ample_util
import contacts
from analyse_run import AmpleResult
import pdb_edit
import pdb_model
import residue_map


root = "/media/data/shared/coiled-coils/ensemble"
runDir = os.getcwd()
os.chdir( runDir )

pfile = os.path.join( runDir, "ar_results_bucc.pkl" )
pfile = "/home/jmht/Documents/test/CC/timings/ar_results_bucc_timings.pkl"
with open( pfile ) as f:
    ensembleResults = cPickle.load( f )

# Hack to add extra attributes
a = AmpleResult()

#for pdbCode in [ "3H7Z" ]:

NI = {}

for r in ensembleResults:
    
    
    print "processing ",r.pdbCode, r.ensembleName
    
    # Need to add the extra attributes
    r.orderedAttrs = a.orderedAttrs
    r.orderedTitles = a.orderedTitles
    r.numAllContacts = None
    r.allContactsOrigin = None
    
    if r.floatingOrigin == True:
        continue
    
    dataDir = os.path.join( root, r.pdbCode )
    workdir = os.path.join( runDir, r.pdbCode )
    if not os.path.isdir( workdir ):
        os.mkdir( workdir )
    
    os.chdir( workdir )
    
    
    if r.pdbCode not in NI:
    
        # Get path to native Extract all the nativePdbInfo from it
        nativePdb = os.path.join( dataDir, "{0}.pdb".format( r.pdbCode ) )
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
        
        # Get the new Info about the native
        nativePdbInfo = pdbedit.get_info( nativePdb )
        
        # Get information on the origins for this spaceGroup
        originInfo = pdb_model.OriginInfo( spaceGroupLabel=nativePdbInfo.crystalInfo.spaceGroup )
        origins = originInfo.nonRedundantAlternateOrigins()
    
        # Get hold of a full model so we can do the mapping of residues
        refModelPdb = os.path.join( dataDir, "models/S_00000001.pdb".format( r.pdbCode ) )
        resSeqMap = residue_map.residueSequenceMap()
        refModelPdbInfo = pdbedit.get_info( refModelPdb )
        resSeqMap.fromInfo( refInfo=refModelPdbInfo,
                            refChainID=refModelPdbInfo.models[0].chains[0], # Only 1 chain in model
                            targetInfo=nativePdbInfo,
                            targetChainID=nativePdbInfo.models[0].chains[0]
                          )
        
        NI[ r.pdbCode ] = (nativePdbInfo,
                           resSeqMap,
                           origins)


    # For all jobs add up phaser and shelxe times to get overall time
    mrbdir = os.path.join( root, r.pdbCode, "ROSETTA_MR_0/MRBUMP/cluster_1" )
    #ensembles = [ os.path.splitext( os.path.basename( l ) )[0] for l in glob.glob( mrbdir + "/*.sub") ]
    #for ensemble in ensembles:
    ensemble = r.ensembleName
        
    mrDir = os.path.join( mrbdir,
                          "search_{0}_mrbump".format( ensemble ),
                          "data",
                          "loc0_ALL_{0}".format( ensemble ),
                          "unmod/mr/phaser"
                            )
    
    phaserPdb = os.path.join( mrDir, "refine", "phaser_loc0_ALL_{0}_UNMOD.1.pdb".format( ensemble ) )
    if os.path.isfile( phaserPdb ):
        
        placedPdbInfo = pdbedit.get_info( phaserPdb )

        ccalc = contacts.Contacts()
        
        ccalc.allContacts( placedPdbInfo=placedPdbInfo,
                           nativePdbInfo=NI[ r.pdbCode ][0],
                           resSeqMap=NI[ r.pdbCode ][1],
                           origins=NI[ r.pdbCode ][2] ,
                           workdir=workdir
                           )
        
        r.numAllContacts = ccalc.numAllContacts
        r.allContactsOrigin = ccalc.allContactsOrigin
        
os.chdir( runDir )

pfile = "ar_results_contacts.pkl"
f = open( pfile, 'w' )
ampleDict = cPickle.dump( ensembleResults, f  )

cpath = os.path.join( runDir, 'results_contacts.csv' )
csvfile =  open( cpath, 'wb')
csvwriter = csv.writer(csvfile, delimiter=',',
                        quotechar='"', quoting=csv.QUOTE_MINIMAL)

header=False
for r in ensembleResults:
    if not header:
        #csvwriter.writerow( r.titlesAsList() )
        csvwriter.writerow( r.valueAttrAsList() )
        header=True
    csvwriter.writerow( r.valuesAsList() )
    
csvfile.close()
