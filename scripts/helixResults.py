#!/usr/bin/env ccp4-python

import csv
import cPickle
import glob
import os
import sys
sys.path.insert(0,"/opt/ample-dev1/python")
sys.path.insert(0,"/opt/ample-dev1/scripts")

from analyse_run import AmpleResult
import phaser_parser
import parse_arpwarp
import parse_buccaneer
import parse_refmac
import parse_shelxe


"""
list all search directories and extract the number of residues from the directory

get mr directory
get phaser log & get data & calculate time
move pdb onto  native origin
get refmac and find rfree
get shelxe and get data and time
move pdb onto  native origin
get buccaneer pdb and get data
move pdb onto  native origin
get arpwarp pdb and get data
move pdb onto  native origin
"""

nativeOrigin=False
if nativeOrigin:
    import matchToNative

dataRoot=""
nativeData=""

results=[]
for pdbCode in [ l.strip() for l in open( os.path.join(dataRoot,"dirs.list") ) if not l.startswith("#") ]:
    
        if pdbCode != "1BYZ":
            continue
    
        mrbumpDir=os.path.abspath(os.path.join(dataRoot,pdbCode))
        
        # Get a list of the jobs
        # For now we just use the submission scripts and assume all have .sh or .sub extension
        jobDirs = [ os.path.splitext( os.path.basename(e) )[0] for e in glob.glob( os.path.join( mrbumpDir, "*.sh") ) ]
        if not len(jobDirs):
            # legacy - try .sub
            jobDirs=[ os.path.splitext( os.path.basename(e) )[0] for e in glob.glob( os.path.join( mrbumpDir, "*.sub") ) ]

        if not len(jobDirs):
            print("Could not extract any results from directory: {0}".format( mrbumpDir ) )
            return False
        
        for jobDir in jobDirs:
            
            mrPdbs=[] # For moving onto native origin
            
            # Sometimes mrbump is appended - other times not
            if jobDir.endswith("_mrbump"):
                jobName=jobDir[7:-7]
            else:
                jobName=jobDir[7:]
                
            mrProgram="phaser"
            # Get the molecular replacement directory
            mrDir=os.path.join( jobDir,
                                  "data",
                                  "loc0_ALL_{0}".format(jobName),
                                  "unmod",
                                  "mr",
                                  mrProgram
                                 )
            
            if not os.path.isdir(mrDir):
                print("missing phaser directory: {0}".format(mrDir))
                continue
            
            result=AmpleResult()
            
            #
            # Phaser and Refmac processing
            #
            phaserLog=os.path.join(mrDir, "{0}_loc0_ALL_{1}_UNMOD.log".format(mrProgram,jobName))
            if os.path.isfile(phaserLog):
                phaserP = phaser_parser.PhaserLogParser(phaserLog, onlyTime=True)
                result.phaserLog    = phaserLog
                result.phaserTime   = phaserP.time
                result.phaserKilled = phaserP.killed

            phaserPdb=os.path.join( mrDir,"refine","{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(mrProgram,jobName))
            if os.path.isfile(phaserPdb):
                mrPdbs.append(phaserPdb)
                phaserP = phaser_parser.PhaserPdbParser(phaserPdb)
                result.phaserLLG = phaserP.LLG
                result.phaserTFZ = phaserP.TFZ
                result.phaserPdb = phaserPdb
                
            refmacLog=os.path.join(mrDir,"refine","refmac_{0}_loc0_ALL_{1}_UNMOD.1.log".format(mrProgram,jobName))
            if os.path.isfile(refmacLog):
                refmacP = parse_refmac.RefmacLogParser(refmacLog)
                result.rfact = refmacP.finalRfact
                result.rfree = refmacP.finalRfree
            
            #
            # Shelxe processing 
            #
            shelxeDir = os.path.join(mrDir,"build","shelxe")
            shelxePdb = os.path.join(shelxeDir,
                                    "shelxe_{0}_loc0_ALL_{1}_UNMOD.pdb".format( mrProgram,jobName))
            if os.path.isfile(shelxePdb):
                mrPdbs.append(shelxePdb)
                
            shelxeLog = os.path.join(shelxeDir,
                                    "shelxe_run.log" )
            
            if os.path.isfile(shelxeLog):
                shelxeP = parse_shelxe.ShelxeLogParser(shelxeLog)
                result.shelxeCC = shelxeP.CC
                result.shelxeAvgChainLength = shelxeP.avgChainLength
                result.shelxeMaxChainLength = shelxeP.maxChainLength
                result.shelxeNumChains= shelxeP.numChains
                
            #
            # Buccaneer Rebuild Processing
            #
            buccaneerPdb = os.path.join(shelxeDir,
                                        "rebuild",
                                        "buccaneer",
                                        "buccSX_output.pdb")
            if os.path.isfile(buccaneerPdb):
                mrPdbs.append(buccaneerPdb)
            
            buccaneerLog = os.path.join(shelxeDir,
                                        "rebuild",
                                        "buccaneer",
                                        "buccaneer.log" )
            
            if os.path.isfile(buccaneerLog):
                bp = parse_buccaneer.BuccaneerLogParser()
                bp.parse( buccaneerLog )
                result.buccFinalRfree = bp.finalRfree
                result.buccFinalRfact = bp.finalRfact
            
            #
            # Arpwarp Rebuild Processing
            #
            arpwarpPdb = os.path.join(shelxeDir,
                                      "rebuild",
                                      "arpwarp",
                                      "refmacSX_output_warpNtrace.pdb")
            if os.path.isfile(arpwarpPdb):
                mrPdbs.append(arpwarpPdb)
                
            arpwarpLog = os.path.join(shelxeDir,
                                      "rebuild",
                                      "arpwarp",
                                      "arpwarp.log")
            
            if os.path.isfile(arpwarpLog):
                ap = parse_arpwarp.ArpwarpLogParser()
                ap.parse(arpwarpLog)
                result.arpWarpFinalRfact=ap.finalRfree
                result.arpWarpFinalRfact=ap.finalRfact
                
            # Save the result object
            results.append(result)
                
            # See if we want to copy all the pdb's onto the same origin as the native
            if nativeOrigin:
                resultDir=os.path.join(nativeData,pdbCode)
                if not os.path.isdir(resultDir):
                    os.mkdir(resultDir)
                if len(mrPdbs):
                    nativeDir=os.path.join(nativeData,pdbCode)
                    nativePdb=os.path.join(nativeDir,pdbCode+".pdb")
                    nativeMtz=os.path.join(nativeDir,pdbCode+"-cad.pdb")
                    matchToNative.run(nativePdb, nativeMtz, mrPdbs, outDir=mrPdbs)
                

runDir=os.getcwd()
pfile = os.path.join( runDir, "results.pkl")
f = open( pfile, 'w' )
cPickle.dump( results, f  )

cpath = os.path.join( runDir, 'results.csv' )
csvfile =  open( cpath, 'wb')
csvwriter = csv.writer(csvfile, delimiter=',',
                        quotechar='"', quoting=csv.QUOTE_MINIMAL)

header=False
for r in results:
    if not header:
        #csvwriter.writerow( r.titlesAsList() )
        csvwriter.writerow( r.valueAttrAsList() )
        header=True
    csvwriter.writerow( r.valuesAsList() )
    
csvfile.close()

sys.exit()

wdir = os.getcwd()
hdir = "/media/data/shared/coiled-coils/ideal_helices/ideal_helices.run2"

helixResults = []

# # unpickle helix dict
# hpkl = os.path.join( hdir, "ensemble2Helix.pkl" )
# with open( hpkl ) as h:
#     e2h = cPickle.load(h)

os.chdir( hdir )
for d in glob.glob("[0-9]*"):
    os.chdir(d)
    for mrd in glob.glob("search_*"):
        rs = mrbump_results.ResultsSummary()
        resultsTable = os.path.abspath( os.path.join( mrd,"results", "resultsTable.dat" ) )
        if not os.path.exists(resultsTable):
            print(" -- Could not find file: {0}".format( resultsTable ) )
            continue
            
        # Extract the result
        result =  rs.parseTableDat(resultsTable)[0]
        
        # Hack to add in other attributes
        result.pdbCode        = d
        helixSequence         = mrd.split("_")[1]
        result.helixSideChain = "_".join( mrd.split("_")[2:] )
        result.helixSequence  = helixSequence
        #result.ensembles      = []
        #for seq in e2h[ result.pdbCode ]:
        #    if seq == helixSequence:
        #        result.ensembles  =  e2h[ result.pdbCode ][ seq ]
                
        helixResults.append( result )

    os.chdir(hdir)
    #break

pfile = os.path.join( wdir, "helix_results.pkl")
with open( pfile, 'w' ) as f:
    cPickle.dump( helixResults, f  )