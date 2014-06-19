'''
Created on Feb 28, 2013

@author: jmht
'''

# python imports
import logging
import multiprocessing
import os
import unittest

# our imports
import clusterize
import mrbump_cmd
import mrbump_results
import workers


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
        if amoptd['submit_qtype'] == "SGE":
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
    
    if amoptd['submit_array']:
        mrBuild.submitArrayJob(job_scripts)
    else:
        for script in job_scripts:
            job_number = mrBuild.submitJob( subScript=script )

    # Monitor the cluster queue to see when all jobs have finished
    mrBuild.monitorQueue()
    
    if amoptd['submit_array']:
        mrBuild.cleanUpArrayJob()
    
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
    """Run ensembling locally"""
    js = workers.JobServer()
    js.setJobs( job_scripts )
    js.start( nproc=amoptd['nproc'],
              early_terminate=amoptd['early_terminate'],
              check_success=check_success )
    return

def check_success( job ):
    """
    Check if a job ran successfully.
    
    Args:
    directory -- directory mr bump ran the job
    
    Returns:
    True if success
    
    Success is assumed as a SHELX CC score of >= SHELXSUCCESS
    """
    
    
    directory, script = os.path.split( job )
    scriptname = os.path.splitext( script )[0]
    rfile = os.path.join(directory, 'search_'+scriptname+'_mrbump','results/resultsTable.dat')
    #print "{0} checking for file: {1}".format(multiprocessing.current_process().name,rfile)
    if not os.path.isfile(rfile):
        print "{0} cannot find results file: {1}".format(multiprocessing.current_process().name,rfile)
        return False
    
    # Results summary object to parse table file
    mrbR = mrbump_results.ResultsSummary()
    
    # Put into order and take top one
    results = mrbR.parseTableDat(rfile)
    mrbR.sortResults(results)
    r = results[0]

    success=False
    rFreeSuccess=0.4
    if r.shelxCC and r.shelxCC != "--" and float(r.shelxCC) >= 25.0:
        success=True
    elif r.buccRfree and r.buccRfree != "--" and float(r.buccRfree) >=rFreeSuccess:
        success=True
    elif r.arpWarpRfree and r.arpWarpRfree != "--" and float(r.arpWarpRfree) >=rFreeSuccess:
        success=True
    elif r.rfree and float(r.rfree) >= rFreeSuccess:
        success=True
        
    return success

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
        
        pdbid = "/home/jmht/t/ample-dev1/examples/toxd-example/ROSETTA_MR_3/MRBUMP/cluster_2/search_All_atom_trunc_0.551637_rad_3_mrbump"
        
        self.assertTrue( check_success( pdbid ) )
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    #unittest.main()
    # Nothing here - see notes in worker
    pass
    
