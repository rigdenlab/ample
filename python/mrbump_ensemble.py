'''
Created on Feb 28, 2013

@author: jmht
'''

# python imports
import logging
import multiprocessing
import subprocess
import os
import time
import unittest

# our imports
import clusterize
import mrbump_cmd

def mrbump_ensemble_cluster( ensembles, mrBuildClusterDir, amoptd, clusterID="X" ):
    """
    Process the list of ensembles using MrBump on a cluster.
    
    Args:
    ensembles -- list of pdb files, each file containing a group of ensembles
    amoptd -- dictionary object containing options
    clusterId -- number/id of the cluster that is being processed
    """
    logger = logging.getLogger()
    logger.info("Running MR and model building on a cluster\n\n")

    mrBuild = clusterize.ClusterRun()
    mrBuild.QTYPE = amoptd['submit_qtype']

    # Reset the queue list
    mrBuild.qList = []
    jobID = 0
    for pdbfile in ensembles:
        mrBuild.mrBuildOnCluster(mrBuildClusterDir, pdbfile, jobID, amoptd )
        jobID = jobID + 1
    
    # Monitor the cluster queue to see when all jobs have finished
    mrBuild.monitorQueue()
    
    # Cleanup code
    #shutil.rmtree(work_dir + '/fine_cluster_' + str(clusterID))
    # shutil.rmtree(work_dir+'/pre_models')
    #for l in os.listdir(work_dir + '/spicker_run'):
    #    if os.path.splitext(l)[1] == 'pdb':
    #        os.remove(work_dir + '/spicker_run/' + l)
    #os.remove(work_dir + '/spicker_run/rep1.tra1')

    # for each_run in os.listdir(mrBuildClusterDir ):
        #   if os.path.isdir(  os.path.join(mrBuildClusterDir, each_run)):
        #      name=re.split('_', each_run)
        #      mrBuildOutputDir=os.path.join(bump_dir, "cluster_run"+str(cluster)+"result"+name[1])
        #      os.mkdir(mrBuildOutputDir)
        #      shutil.move (os.path.join(mrBuildClusterDir, each_run, "search_"+name[1]+"_mrbump","phaser_shelx" ),mrBuildOutputDir  )
        #      shutil.move (os.path.join(mrBuildClusterDir, each_run, "search_"+name[1]+"_mrbump","molrep_shelx" ),mrBuildOutputDir  )
        #      shutil.move (os.path.join(mrBuildClusterDir, each_run, "search_"+name[1]+"_mrbump","data" ),mrBuildOutputDir  )
        #      shutil.move (os.path.join(mrBuildClusterDir, each_run, "logs" ),mrBuildOutputDir  )
        # shutil.rmtree(mrBuildClusterDir)
    
##End mrbump_ensemble_cluster

def mrbump_ensemble_local( ensembles, amoptd, clusterID="X" ):
    """
    Process the list of ensembles using MrBump on a local machine
    
    Args:
    ensembles -- list of pdb files, each file containing a group of ensembles
    amoptd -- dictionary with job options
    clusterId -- number/id of the cluster that is being processed
    """

    logger = logging.getLogger()
    logger.info("Running MR and model building on a local machine")

    # Queue to hold the jobs we want to run
    queue = multiprocessing.Queue()
    
    # Create all the run scripts and add to the queue
    logger.info("Generating MRBUMP runscripts in: {0}".format( os.getcwd() ) )
    for pdbfile in ensembles:
        # Name is filename without extension
        name = os.path.splitext( os.path.split(pdbfile)[1] )[0]
        script_path = create_jobscript(name, pdbfile, amoptd)
        queue.put(script_path)

    # Now start the jobs
    processes = []
    for i in range( amoptd['nproc'] ):
        process = multiprocessing.Process(target=worker, args=( queue, amoptd['early_terminate'] ) )
        process.start()
        processes.append(process)

    
    # Loop through the processes checking if any are done
    timeout=10*60
    timeout=1*60
    while len(processes):
        
        for i, process in enumerate(processes):
            
            # Join process for timeout seconds and if we haven't finished by then
            # move onto the next process
            process.join(timeout)
            
            if not process.is_alive():
                #print "CHECKING COMPLETED PROCESS {0} WITH EXITCODE {1}".format(process,process.exitcode)
                # Remove from processes to check
                del processes[i]
                
                # Finished so see what happened
                if process.exitcode == 0 and amoptd['early_terminate']:
                    if not queue.empty():
                        logger.info( "Process {0} was successful so removing remaining jobs from queue".format(process.name) )
                        # Remove all remaining processes from the queue. We do this rather than terminate the processes
                        # as terminating leaves the MRBUMP processes running. This way we hang around until all our
                        # running processes have finished
                        while not queue.empty():
                            job = queue.get()
                            logger.debug( "Removed job [{0}] from queue".format(job) )
                        
    # need to wait here as sometimes it takes a while for the results files to get written
    time.sleep(3)
##End mrbump_ensemble_local

def create_jobscript( name, pdb, amoptd, directory=None ):
    """
    Create the script to run MrBump for this PDB.
    
    Args:
    name -- used to identify job and name the run script
    pdb -- the path to the pdb file for this job
    amoptd -- dictionary with job options
    directory -- directory to write script to - defaults to cwd
    
    Returns:
    path to the script
    """
    
    # Path
    if not directory:
        directory = os.getcwd()
    script_path = directory + os.sep + name + '.sub'
    
    # Write script
    job_script = open(script_path, "w")
    job_script.write('#!/bin/sh\n')
    jobstr = mrbump_cmd.mrbump_cmd( amoptd, jobid=name, ensemble_pdb=pdb )
    job_script.write(jobstr)
    job_script.close()
    
    # Make executable
    os.chmod(script_path, 0o777)
    
    return script_path
##End create_jobscript

def worker( queue, early_terminate=False ):
    """
    Worker process to run MrBump jobs until no more left.
    
    Args:
    queue -- a python Queue object
    early_terminate -- bool - terminate on first success or continue running
    
    Returns:
    0 if molecular replacement worked
    1 if nothing found
    
    We keep looping, removing jobs from the queue until there are no more left.
    
    REM: This needs to import the main module that it lives in so maybe this should
    live in a separate module?
    """
    
    while True:
        if queue.empty():
            #print "worker {0} got empty queue {1}".format(os.getpid(),e)
            break

        mrb_script = queue.get()
        
        # Got a script so run
        # Get name from script
        name = os.path.splitext( os.path.split(mrb_script)[1] )[0]
        logfile = name + ".log"
        f = open( logfile, "w")
        print "Worker {0} running job {1}".format(multiprocessing.current_process().name, name)
        retcode = subprocess.call( [ mrb_script ], stdout=f, stderr=subprocess.STDOUT, cwd=None )
        
        # Can we use the retcode to check?
        # REM - is retcode object
        if retcode != 0:
            print "WARNING! Worker {0} got retcode {1}".format(multiprocessing.current_process().name, retcode )
        
        # Run directory is name of script with search_ prepended and _mrbump appended
        directory = os.path.join( os.getcwd(), "search_" + name + "_mrbump" )
        
        # Now check the result if early terminate
        if early_terminate:
            if check_success( directory ):
                print "Worker {0} job succeeded".format(multiprocessing.current_process().name)
                return 0
        
    #print "worker {0} FAILED!".format(multiprocessing.current_process().name)
    return 1
##End worker

def check_success( directory ):
    """
    Check if a job ran successfully.
    
    Args:
    directory -- directory mr bump ran the job
    
    Returns:
    True if success
    
    Success is assumed as a SHELX CC score of >= SHELXSUCCESS
    """
    
    SHELXSUCCESS = 25.0
    
    rfile = os.path.join(directory, 'results/resultsTable.dat')
    #print "{0} checking for file: {1}".format(multiprocessing.current_process().name,rfile)
    if not os.path.isfile(rfile):
        print "{0} cannot find results file: {1}".format(multiprocessing.current_process().name,rfile)
        return False
        
    f = open(rfile, 'r')
    
    # First line is
    hline =  f.readline().strip()
    headers = hline.split()
    
    # For now we assume we are using SHELXE to check results
    # otherwise we would check 'final_Rfree'
    scol = headers.index('SHELXE_CC')
    
    for line in f:
        fields = line.strip().split()
        score = float(fields[scol])
        if score >= SHELXSUCCESS:
            return True
    
    # Nothing good enough
    return False

class Test(unittest.TestCase):


    def XtestLocal(self):
        
        import glob
        
        d = {
             'mtz' : '/opt/ample-dev1/examples/toxd-example/1dtx.mtz',
             'fasta' : '/opt/ample-dev1/examples/toxd-example/toxd_.fasta',
             'mrbump_programs' : ' molrep ',
             'use_buccaneer' : False ,
             'buccaneer_cycles' : 5,
             'use_arpwarp' : False,
             'arpwarp_cycles': 10,
             'use_shelxe' : True,
             'shelx_cycles' : 10,
             'FREE' : 'FreeR_flag',
             'F' : 'FP',
             'SIGF' : 'SIGFP',
             'nproc' : 3,
             'domain_all_chains_pdb' : None,
             'ASU': None,
             'mr_keys': [],
             'early_terminate' : False
             }
        
        ensemble_dir = "/home/Shared/ample-dev1/examples/toxd-example/ROSETTA_MR_0/ensembles_1"
        ensembles = []
        for infile in glob.glob( os.path.join( ensemble_dir, '*.pdb' ) ):
            ensembles.append(infile)
            
        mrbump_ensemble_local( ensembles, d )
            
    
    def XtestSuccess(self):
        
        dir = "/home/jmht/t/ample-dev1/examples/toxd-example/ROSETTA_MR_3/MRBUMP/cluster_2/search_All_atom_trunc_0.551637_rad_3_mrbump"
        
        self.assertTrue( check_success( dir ) )
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    #unittest.main()
    # Nothing here - see notes in worker
    pass
    
