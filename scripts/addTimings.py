#!/usr/bin/env ccp4-python
import sys
sys.path.insert(0, "/opt/ample-dev1/python" )
sys.path.insert(0, "/opt/ample-dev1/scripts" )


import cPickle
import csv

sys.path.append("/usr/lib/python2.7/dist-packages")
from dateutil import parser
import glob
import os

import ample_util
from analyse_run import AmpleResult
#import phaser_parser
import parse_shelxe

def sgeLogTime( log ):
    with open( log, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("Started at "):
                start = parser.parse( line[10:] )
            if line.startswith("Results reported at "):
                end = parser.parse( line[19:] )
                
            if line.startswith("CPU time"):
                ctime = float( line.split()[3] )
                break
    
    
    # Calculate walltime
    wtime = end - start
    wtime = wtime.seconds
    ctime = int( ctime )
    return wtime, ctime

def ampleLogTime( log ):
    with open( log, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("ALL DONE"):
                thours =  float( line.split()[3] ) 
    
    return thours * 60 * 60

def logTime(log):
    with open( log, 'r') as f:
        lines = f.readlines()
        
        # times are first two chunks of file
        start = " ".join(lines[0].strip().split()[0:2])
        # strip milliseconds - why the FUCK!!! is the default log format a time that can't be parsed with
        # standard python tools?!?!?!?!?!
        start = start.split(",")[0]
        
        l=None
        for i in range(1,5):
            l = lines[-i].strip()
            if l:
                break
        assert l
        end = " ".join(l.split()[0:2])
        end = end.split(",")[0]
        
        tstart = parser.parse( start )
        tend = parser.parse( end )
        
    delta = tend-tstart
    return delta.seconds



#l = "/media/data/shared/coiled-coils/ensemble/ensemble_redo_failures1/1G1J/ROSETTA_MR_0/ensemble.log"
#t = logTime(l)
#print "TIME ",t
#sys.exit()


e1root = "/media/seagate/coiled-coils/ensemble/ensemble.run1"

runDir = os.getcwd()
os.chdir( runDir )

pfile = os.path.join( runDir, "final_results.pkl" )
with open( pfile ) as f:
    ensembleResults = cPickle.load( f )


# Map targets to directories
pdb2dir = {}
for jd in [ l.strip() for l in open( "/media/data/shared/coiled-coils/ensemble/final_results/dirs.list") if not l.startswith("#") ]:
    directory = "/".join(jd.split("/")[0:-1])
    pdbCode = jd.split("/")[-1]
    pdb2dir[pdbCode]=directory

# Hack to add extra attributes
a = AmpleResult()


#for pdbCode in [ "3H7Z" ]:
for r in ensembleResults:
    
    #if r.pdbCode not in pdb2dir:
    #    continue
    
    print "processing ",r.pdbCode, r.ensembleName
    # Need to add the extra attributes
    r.orderedAttrs = a.orderedAttrs
    r.orderedTitles = a.orderedTitles
    
    dataRoot=pdb2dir[r.pdbCode]
    # Always use the old models as we dont' have the times for the redo ones
    #if dataRoot == "/media/data/shared/coiled-coils/ensemble/ensemble.run2":
    #    modelsDir=os.path.join(e1root,r.pdbCode,"models")
    #else:
    #    modelsDir=os.path.join(dataRoot,r.pdbCode,"models")
    
    modelsDir=os.path.join(e1root,r.pdbCode,"models")
    
    # Find all fragment logs and add up times
    fdir = os.path.join( e1root, r.pdbCode, "fragments" )
    ftime = 0.0
    for flog in glob.glob( fdir + "/frags_*log*"):
        w, c = sgeLogTime( flog )
        ftime += c
    
    #print "GOT ftime ",ftime
    r.fragmentTime = ftime

    # Parse model time
    mlog = glob.glob( modelsDir + "/models_*")[0]
    w, mtime = sgeLogTime( mlog )
    
    #print "MTIME ",mtime
    r.modelTime = mtime
    
    # Get the ensembling time from the ensemble log
    if dataRoot == "/media/data/shared/coiled-coils/ensemble/ensemble.run2":
        elog=os.path.join(e1root,r.pdbCode,"ROSETTA_MR_0","ensemble.log")
        etime = logTime(elog)
        # First run ran 3 ensembles
        etime = etime/3
    else:
        elog=os.path.join(dataRoot,r.pdbCode,"ROSETTA_MR_0","ensemble.log")
        etime = logTime(elog)
    
    # Parse ample log to get ample time
    #alog = os.path.join( dataRoot, r.pdbCode, "run_ample.sh.out" )
    #atime = ampleLogTime( alog )
    
    #print "ETIME ",etime
    r.ensembleTime = etime

    # For all jobs add up phaser and shelxe times to get overall time
    if dataRoot == "/media/data/shared/coiled-coils/ensemble/ensemble.run2":
        mrbdir = os.path.join(dataRoot,r.pdbCode)
    else:
        mrbdir = os.path.join( dataRoot, r.pdbCode, "ROSETTA_MR_0/MRBUMP/cluster_1" )
    
    mrDir = os.path.join( mrbdir,
                          "search_{0}_mrbump".format( r.ensembleName ),
                          "data",
                          "loc0_ALL_{0}".format( r.ensembleName ),
                          "unmod/mr/phaser"
                            )
    #Already calculated the phaser log time    
#     phaserLog = os.path.join( mrDir, "phaser_loc0_ALL_{0}_UNMOD.log".format( r.ensembleName ) )
#     ptime = 0.0
#     if os.path.isfile( phaserLog ):
#         phaserP = phaser_parser.PhaserLogParser( phaserLog )
#         ptime = phaserP.time
    
    shelxeLog = os.path.join( mrDir, "build/shelxe/shelxe_run.log" )
    stime = 0.0
    if os.path.isfile( shelxeLog ):
        shelxeP = parse_shelxe.ShelxeLogParser( shelxeLog )
        stime = shelxeP.cputime
    
    r.shelxeTime = stime
    
    #print "PTIME ",r.phaserTime
    #print "STIME ",stime



pfile = ample_util.filename_append( pfile, astr="timings")
f = open( pfile, 'w' )
ampleDict = cPickle.dump( ensembleResults, f  )

cpath = os.path.join( runDir, 'final_results_timings.csv' )
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

