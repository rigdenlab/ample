#!/usr/bin/env ccp4-python

import csv
import cPickle
import glob
import os
import sys
sys.path.insert(0,"/opt/ample-dev1/python")
sys.path.insert(0,"/opt/ample-dev1/scripts")

import phaser_parser
import parse_arpwarp
import parse_buccaneer
import parse_refmac
import parse_shelxe


def alignPdb(nativeMap,mrPdb,outPdb=None,outDir=None):
    origin =  phenixer.ccmtzOrigin( nativeMap, mrPdb )
    
    # offset.pdb is the mrPdb moved onto the new origin
    offsetPdb = "offset.pdb"
    print "Found origin: {0}\nOffset pdb is: {1}".format( origin, offsetPdb )
    
    # Run csymmatch to map offsetted pdb onto native
    if outPdb is None:
        outPdb = ample_util.filename_append( filename=mrPdb, astr="csymmatch", directory=outDir )
    print "Running csymmatch to wrap {0} onto native {1}".format( offsetPdb, nativePdb )
    csymmatch.Csymmatch().run( refPdb=nativePdb, inPdb=offsetPdb, outPdb=outPdb, originHand=False )
    
    return

def matchResult(result,nativeMap,outDir):

    assert os.path.isfile(nativeMap)
    assert os.path.isdir(outDir)
    
    if os.path.isfile(result.phaserPdb):
        alignPdb(nativeMap,result.phaserPdb,outPdb=None,outDir=outDir)
    if os.path.isfile(result.shelxePdb):
        alignPdb(nativeMap,result.shelxePdb,outPdb=None,outDir=outDir)
    if os.path.isfile(result.buccaneerPdb):
        outPdb=os.path.join(outDir,"buccaneer_{0}_csymmatch.pdb".format(result.jobName))
        alignPdb(nativeMap,result.buccaneerPdb,outPdb=outPdb,outDir=outDir)
    if os.path.isfile(result.arpwarpPdb):
        outPdb=os.path.join(outDir,"arpwarp_{0}_csymmatch.pdb".format(result.jobName))
        alignPdb(nativeMap,result.arpwarpPdb,outPdb=outPdb,outDir=outDir)
    
    return

def processMrDirectory(result):
    
    assert os.path.isdir(result.mrDir)
    assert result.mrProgram

    #
    # Phaser and Refmac processing
    #
    phaserLog=os.path.join(result.mrDir, "{0}_loc0_ALL_{1}_UNMOD.log".format(result.mrProgram,result.jobName))
    if os.path.isfile(phaserLog):
        phaserP = phaser_parser.PhaserLogParser(phaserLog, onlyTime=True)
        #result.phaserLog    = phaserLog
        result.phaserTime   = phaserP.time
        result.phaserKilled = phaserP.killed

    phaserPdb=os.path.join( result.mrDir,"refine","{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(result.mrProgram,result.jobName))
    if os.path.isfile(phaserPdb):
        phaserP = phaser_parser.PhaserPdbParser(phaserPdb)
        result.phaserLLG = phaserP.LLG
        result.phaserTFZ = phaserP.TFZ
        result.phaserPdb = phaserPdb
        
    refmacLog=os.path.join(result.mrDir,"refine","refmac_{0}_loc0_ALL_{1}_UNMOD.log".format(result.mrProgram,result.jobName))
    if os.path.isfile(refmacLog):
        refmacP = parse_refmac.RefmacLogParser(refmacLog)
        result.rfact = refmacP.finalRfact
        result.rfree = refmacP.finalRfree
    
    #
    # Shelxe processing 
    #
    shelxeDir = os.path.join(result.mrDir,"build","shelxe")
    shelxePdb = os.path.join(shelxeDir,
                            "shelxe_{0}_loc0_ALL_{1}_UNMOD.pdb".format(result.mrProgram,result.jobName))
    if os.path.isfile(shelxePdb):
        pass
        
    shelxeLog = os.path.join(shelxeDir,
                            "shelxe_run.log" )
    
    if os.path.isfile(shelxeLog):
        shelxeP = parse_shelxe.ShelxeLogParser(shelxeLog)
        result.shelxeCC = shelxeP.CC
        result.shelxeAvgChainLength = shelxeP.avgChainLength
        result.shelxeMaxChainLength = shelxeP.maxChainLength
        result.shelxeNumChains= shelxeP.numChains
        result.shelxeTime=shelxeP.cputime
        
    #
    # Buccaneer Rebuild Processing
    #
    buccaneerPdb = os.path.join(shelxeDir,
                                "rebuild",
                                "buccaneer",
                                "buccSX_output.pdb")
    if os.path.isfile(buccaneerPdb):
        pass
    
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
        pass
        
    arpwarpLog = os.path.join(shelxeDir,
                              "rebuild",
                              "arpwarp",
                              "arpwarp.log")
    
    if os.path.isfile(arpwarpLog):
        ap = parse_arpwarp.ArpwarpLogParser()
        ap.parse(arpwarpLog)
        result.arpWarpFinalRfact=ap.finalRfact
        result.arpWarpFinalRfree=ap.finalRfree
    
    return

class HelixResult(object):
    """Results for an ample solution"""
    
    def __init__(self):

        # The attributes we will be holding
        self.orderedAttrs = [ 
                              'pdbCode',
                              'polyaLength',
                              'phaserLLG',
                              'phaserTFZ',
                              'phaserTime',
                              'phaserKilled',
                              'rfact',
                              'rfree',
                              'buccFinalRfact',
                              'buccFinalRfree',
                              'arpWarpFinalRfact',
                              'arpWarpFinalRfree',
                              'shelxeCC',
                              'shelxeAvgChainLength',
                              'shelxeMaxChainLength',
                              'shelxeNumChains',
                              'shelxeTime'
                              ]
        
        # The matching titles
        self.orderedTitles = [  
                              'pdbCode',
                              'polyaLength',
                              'phaserLLG',
                              'phaserTFZ',
                              'phaserTime',
                              'phaserKilled',
                              'rfact',
                              'rfree',
                              'buccFinalRfact',
                              'buccFinalRfree',
                              'arpWarpFinalRfact',
                              'arpWarpFinalRfree',
                              'shelxeCC',
                              'shelxeAvgChainLength',
                              'shelxeMaxChainLength',
                              'shelxeNumChains',
                              'shelxeTime'
                                 ]

        # Things not to output
        self.skip = [ ]
        
        # Set initial values
        for a in self.orderedAttrs:
            setattr( self, a, None )
        
        return
    
    def valueAttrAsList(self):
        return [ a for a in self.orderedAttrs if a not in self.skip ]

    def valuesAsList(self):
        return [ getattr(self, a) for a in self.orderedAttrs if a not in self.skip ]
    
    def titlesAsList(self):
        return [ t for i, t in enumerate( self.orderedTitles ) if self.orderedAttrs[i] not in self.skip ]
    
    def asDict(self):
        """Return ourselves as a dict"""
        d = {}
        for a in self.orderedAttrs:
            d[ a ] = getattr(self, a)
        return d
    
    def __str__(self):
        
        s = ""
        for i, t in enumerate( self.orderedTitles ):
            if self.orderedAttrs[i] not in self.skip:
                s += "{0:<26} : {1}\n".format( t, getattr( self, self.orderedAttrs[i] ) )
        return s
    
# End HelixResult


matchNative=True
if matchNative:
    import ample_util
    import csymmatch
    import phenixer
    import MTZ_parse

dataRoot="/home/jmht/foo"
analysisDir="/home/jmht/foo/analysis"

results=[]
runDir=os.getcwd()
for pdbCode in [ l.strip() for l in open( os.path.join(dataRoot,"dirs.list") ) if not l.startswith("#") ]:
    
        mrbumpDir=os.path.abspath(os.path.join(dataRoot,pdbCode))
        outDir=os.path.join(analysisDir,pdbCode)
        if not os.path.isdir(outDir):
            os.mkdir(outDir)
        
        # Get a list of the jobs
        # For now we just use the submission scripts and assume all have .sh or .sub extension
        jobs = [ os.path.splitext( os.path.basename(e) )[0] for e in glob.glob( os.path.join( mrbumpDir, "*.sh") ) ]
        if not len(jobs):
            # legacy - try .sub
            jobs=[ os.path.splitext( os.path.basename(e) )[0] for e in glob.glob( os.path.join( mrbumpDir, "*.sub") ) ]

        if not len(jobs):
            print("Could not extract any results from directory: {0}".format( mrbumpDir ) )
            continue
        
        if matchNative:
            # Get paths to files
            nativePdb=os.path.join(analysisDir,pdbCode+".pdb")
            nativeMtz=os.path.join(analysisDir,pdbCode+"-cad.mtz")
            
            # Create the map
            mtzp = MTZ_parse.MTZ_parse()
            mtzp.run_mtzdmp( nativeMtz )
            nativeMap = phenixer.generateMap(nativeMtz,
                                             nativePdb,
                                             FP= mtzp.F,
                                             SIGFP=mtzp.SIGF,
                                             FREE=mtzp.FreeR_flag)
        
        for jobName in jobs:
            
            # Sometimes mrbump is appended - other times not
            jobDir=os.path.join(mrbumpDir,"search_{0}_mrbump".format(jobName))
            if not os.path.isdir(jobDir):
                jobDir=os.path.join(mrbumpDir,"search_{0}".format(jobName))
            if not os.path.isdir(jobDir):
                print("MISSING JOB DIRECTORY: {0}".format(jobDir))
                continue
                
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
            
            # Create empty result object
            result=HelixResult()
            
            result.pdbCode=pdbCode
            result.jobName=jobName
            result.mrProgram=mrProgram
            result.mrDir=mrDir
            
            # Number of residues
            result.polyaLength=int(jobName.split("_")[1])
            
            processMrDirectory(result)
            results.append(result)
            
            if matchNative:
                alignPdbs(result,nativeMap,outDir)
                
        # End of individual job loop

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