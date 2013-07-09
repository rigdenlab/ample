'''
Created on Feb 28, 2013

@author: jmht
'''

# python imports
import logging
import multiprocessing
import subprocess
import os
import sys
import time
import unittest

# our imports
import ample_util
import clusterize
import mrbump_cmd


def generate_jobscripts( ensemble_pdbs, amoptd ):
    """Write the MRBUMP shell scripts for all the ensembles.

    Args:
    ensemble_pdbs -- list of the ensembles, each a single pdb file
    amoptd -- dictionary with job options
    
    The split_mr option used here was added for running on the hartree
    wonder machine where job times were limited to 12 hours, but is left
    in in case it's of use elsewhere.
    """
    
    # Remember programs = also used for looping
    if amoptd['split_mr']:
        mrbump_programs = amoptd['mrbump_programs']
        nproc = amoptd['nproc']
    
    job_scripts = []
    for ensemble_pdb in ensemble_pdbs:
        
        # Get name from pdb path
        name = os.path.splitext( os.path.basename( ensemble_pdb ) )[0]
        
        # May need to run MR separately
        if amoptd['split_mr']:
            # create multiple jobs
            for program in mrbump_programs:
                jname = "{0}_{1}".format( name, program )
                amoptd['mrbump_programs'] = [ program ]
                # HACK - molrep only runs on a single processor
                # Can't do this any more as any job < 16 can't run on the 12 hour queue
                #if program == "molrep":
                #    amoptd['nproc'] = 1
                script = write_jobscript( name=jname, pdb=ensemble_pdb, amoptd=amoptd )
                #amoptd['nproc'] = nproc
                job_scripts.append( script )
        else:
            # Just run as usual
            script = write_jobscript( name=name, pdb=ensemble_pdb, amoptd=amoptd )
            job_scripts.append( script )
            
    # Reset amoptd
    if amoptd['split_mr']:
       amoptd['mrbump_programs'] = mrbump_programs
            
    if not len( job_scripts ):
        msg = "No job scripts created!"
        logging.critical( msg )
        raise RuntimeError, msg
    
    return job_scripts
        
def write_jobscript( name, pdb, amoptd, directory=None ):
    """
    Create the script to run MrBump for this PDB.
    
    Args:
    name -- used to identify job and name the run script
    pdb -- the path to the pdb file for this job
    amoptd -- dictionary with job options
    directory -- directory to write script to - defaults to cwd
    
    Returns:
    path to the script
    
    There is an issue here as the code to add the parallel job submission
    script header is required here, so we create a ClusterRun object.
    Should think about a neater way to do this rather then split the parallel stuff
    across two modules.
    """
    
    # Path
    if not directory:
        directory = os.getcwd()
    script_path = directory + os.sep + name + '.sub'
    
    # Write script
    job_script = open(script_path, "w")
    
    # If on cluster, insert queue header 
    if amoptd['submit_cluster']:
        # Messy - create an instance to get the script header. Could pass one in to save creating one each time
        mrBuild = clusterize.ClusterRun()
        mrBuild.QTYPE = amoptd['submit_qtype']
        logFile=os.path.join( directory, name + ".log" )
        script_header = mrBuild.subScriptHeader( nProc=amoptd['nproc'], logFile=logFile, jobName=name)
        job_script.write( script_header )
        job_script.write("pushd " + directory + "\n\n")
        # Required on the RAL cluster as the default tmp can be deleted on the nodes
        if amoptd['cluster_qtype'] == "SGE":
            job_script.write("setenv CCP4_SCR $TMPDIR\n\n")

    else:
        job_script.write('#!/bin/sh\n')
    
    jobstr = mrbump_cmd.mrbump_cmd( amoptd, jobid=name, ensemble_pdb=pdb )
    job_script.write(jobstr)
    
    if amoptd['submit_cluster']:
        job_script.write('\n\npopd\n\n')
        
    job_script.close()
    
    # Make executable
    os.chmod(script_path, 0o777)
    
    return script_path
##End create_jobscript

def mrbump_ensemble_cluster( job_scripts, amoptd ):
    """
    Process the list of ensembles using MrBump on a cluster.
    
    Args:
    job_scripts -- list of scripts to run mrbump
    amoptd -- dictionary object containing options
    """
    logger = logging.getLogger()
    logger.info("Running MR and model building on a cluster\n\n")

    mrBuild = clusterize.ClusterRun()
    mrBuild.QTYPE = amoptd['submit_qtype']
    for script in job_scripts:
        job_number = mrBuild.submitJob( subScript=script )

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

def mrbump_ensemble_local( job_scripts, amoptd ):
    """
    Process the list of ensembles using MrBump on a local machine
    
    Args:
    job_scripts -- list of scripts to run mrbump
    amoptd -- dictionary with job options
    """

    logger = logging.getLogger()
    logger.info("Running MR and model building on a local machine")

    # Queue to hold the jobs we want to run
    queue = multiprocessing.Queue()
    
    # Add jobs to the queue
    #logger.info("Generating MRBUMP runscripts in: {0}".format( os.getcwd() ) )
    for script in job_scripts:
        queue.put(script)

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
    return

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
        name = os.path.splitext( os.path.basename(mrb_script) )[0]
        #logfile = name + ".log"
        #f = open( logfile, "w")
        print "Worker {0} running job {1}".format(multiprocessing.current_process().name, name)
        
        retcode = ample_util.run_command( [ mrb_script ], logfile=name + ".log", dolog=False)
        #retcode = subprocess.call( [ mrb_script ], stdout=f, stderr=subprocess.STDOUT, cwd=None )
        
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
                #return 0
                sys.exit(0)
        
    #print "worker {0} FAILED!".format(multiprocessing.current_process().name)
    #return 1
    sys.exit(1)
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
    
    SHELX_SUCCESS = 25.0
    
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
    if 'SHELXE_CC' in headers:
        scol = headers.index('SHELXE_CC')
    else:
        return False
    
    for line in f:
        fields = line.strip().split()
        score = float(fields[scol])
        if score >= SHELX_SUCCESS:
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
        
        tdir = "/home/jmht/t/ample-dev1/examples/toxd-example/ROSETTA_MR_3/MRBUMP/cluster_2/search_All_atom_trunc_0.551637_rad_3_mrbump"
        
        self.assertTrue( check_success( tdir ) )
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    #unittest.main()
    # Nothing here - see notes in worker
    pass
    
