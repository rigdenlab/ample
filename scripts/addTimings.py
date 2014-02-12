import sys
sys.path.append( "/opt/ample-dev1/python" )
sys.path.append( "/opt/ample-dev1/scripts" )

import cPickle
import csv
from dateutil import parser
import datetime
import glob
import os

import ample_util
from analyse_run import AmpleResult
import phaser_parser
import shelxe_log

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
    
    return thours * 60

root = "/media/data/shared/coiled-coils/ensemble"
for pdbCode in [ "3H7Z" ]:
    
    # Find all fragment logs and add up times
    fdir = os.path.join( root, pdbCode, "fragments" )
    ftime = 0.0
    for flog in glob.glob( fdir + "/frags_*"):
        w, c = sgeLogTime( flog )
        ftime += c
    
    print "GOT ftime ",ftime

    # Parse model time
    mdir = os.path.join( root, pdbCode, "models" )
    mlog = glob.glob( mdir + "/models_*")[0]
    w, mtime = sgeLogTime( mlog )
    
    print "MTIME ",mtime
    

    # Parse ample log to get ample time
    alog = os.path.join( root, pdbCode, "run_ample.sh.out" )
    atime = ampleLogTime( alog )
    
    print "aTIME ",atime

    # For all jobs add up phaser and shelxe times to get overall time
    mrbdir = os.path.join( root, pdbCode, "ROSETTA_MR_0/MRBUMP/cluster_1" )
    ensembles = [ os.path.splitext( os.path.basename( l ) )[0] for l in glob.glob( mrbdir + "/*.sub") ]
    
    for ensemble in ensembles:
        
        mrDir = os.path.join( mrbdir,
                              "search_{0}_mrbump".format( ensemble ),
                              "data",
                              "loc0_ALL_{0}".format( ensemble ),
                              "unmod/mr/phaser"
                                )
        
        phaserLog = os.path.join( mrDir, "phaser_loc0_ALL_{0}_UNMOD.log".format( ensemble ) )
        ptime = 0.0
        if os.path.isfile( phaserLog ):
            phaserP = phaser_parser.PhaserLogParser( phaserLog )
            ptime = phaserP.time
        
        shelxeLog = os.path.join( mrDir, "build/shelxe/shelxe_run.log" )
        stime = 0.0
        if os.path.isfile( shelxeLog ):
            shelxeP = shelxe_log.ShelxeLogParser( shelxeLog )
            stime = shelxeP.cputime
        
        print "PTIME ",ptime
        print "STIME ",stime


