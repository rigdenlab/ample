#!/usr/bin/env python
'''
Created on May 17, 2013

@author: jmht

Simple taskfarming harness for LSB.

Uses LSB environment variables and a shared filesystem
to control a set of jobs.
'''

import glob
import os
import subprocess
import sys
import time


# Each node sets these
#hostfile = os.environ['LSB_DJOB_HOSTFILE']
hosts = os.environ['LSB_HOSTS'].split()
nhosts = len(hosts)
jobid = int(os.environ['LSB_JOBID'])
hostname = subprocess.check_output(["hostname"]).strip()
BARRIER_FILE="barrier.{}".format( jobid )
BARRIER_SLEEP=5

lockstem="ample_lockfile.{0}".format( jobid )
lockfile="{0}.{1}".format( lockstem, hostname )


def root():
    """Return true if we are the first node"""
    global hostname, hosts
    if hostname == hosts[0]:
        return True
    return False

def setup():
    """Set things up"""
    global hosts, nhosts, jobid, BARRIER_FILE
    print "Root setting things up ",hosts
    f = open( BARRIER_FILE, "w" )
    f.close()
    return


def getLock():
    """Get the lock"""
    
    global hostname, nhosts, jobid, lockstem, lockfile
    
    waitInterval=0.2
    # See if there are any lockfiles
    while len( glob.glob( lockstem+"*") ):
        time.sleep( waitInterval )
    
    # no lockfiles so we create ours
    f = open( lockfile, 'w' )
    f.close()
    
    # Check if any others were create while ours was being written
    lockfiles =  glob.glob( lockstem+"*" )
    if len( lockfiles ) > 1:
        print "GOT DEADLOCK!!"
        # The first host gets the lock
        lockhost = [ l.split(".")[-1] for l in lockfiles ][0]
        
        if not hostname == lockhost:
            releaseLock()
            return False
    
    return True

def releaseLock():
    """Delete the lockfile"""
    os.unlink(lockfile)
    return

def popJob():
    """Return and remove the next job in the job list"""
    global hosts, nhosts, jobid
    
    while not getLock():
        time.sleep(1)
    
    # Get a list of all jobs
    jobfile = "alljobs.list"
    alljobs = [ l.strip() for l in open( jobfile, 'r' ) ]
    njobs= len( alljobs )
    
    # We've finished
    if not njobs:
        releaseLock()
        return False
    
    # Write all bar the last back to file
    f = open( jobfile, 'w' )
    for j in alljobs[:-1]:
        f.write( j +"\n" )
    f.close()
    
    # Release the lock
    releaseLock()
    
    # Return the last job
    job = alljobs[ njobs-1 ]
    return job


def getJobsDict():
    """Create a dictionary of the jobs for each node.
    Will need if we create serial jobs
    """
    
    global hosts, nhosts, jobid
    
    # Get a list of all jobs
    jobfile="alljobs.list"
    alljobs = [ l.strip() for l in open( jobfile, 'r') ]
    njobs=len(alljobs)

    # split jobs between hosts
    jobsPerHost = njobs / nhosts  
    remainder = njobs % nhosts
    
    jobs = {} # maps host : list of jobs
    count=0
    for i in range( 0, njobs-remainder, jobsPerHost ):
        jobs[ hosts[count] ] = alljobs[ i : i+jobsPerHost ]
        count+=1
    
    # Now separate out remainder
    if remainder > 0:
        remainJobs = alljobs[ njobs-remainder: ]
        for i, h in enumerate( jobs.keys() ):
            if i < remainder:
                jobs[ h ].append( remainJobs[ i] )
    return jobs

def checkStart():
    """Check nothing from a previous run is left over"""
    global BARRIER_FILE, lockfile
    if os.path.exists( BARRIER_FILE ):
        print "BARRIER_FILE ALREADY EXISTS: {0}".format( BARRIER_FILE )
        os.unlink( BARRIER_FILE )
    if os.path.exists( lockfile ):
        print "lockfile ALREADY EXISTS: {0}".format( lockfile )
        os.unlink( lockfile )
    return


# 
# if root():
#     setup()
# else:
#     while not os.path.exists( BARRIER_FILE ):
#         time.sleep( BARRIER_SLEEP )
#         
# print "Host {0} carrying on".format(hostname)

checkStart()

job = popJob()
while job:
    print "Host {0} got job {1}".format( hostname, job)
    time.sleep(2)
    job = popJob()
    
print "Host {0} finished".format( hostname )

