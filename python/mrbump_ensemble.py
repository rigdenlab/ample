'''
Created on Feb 28, 2013

@author: jmht
'''

# python imports
import logging
import multiprocessing
import subprocess
import os
import Queue
import re
import shutil
import sys
import unittest

# our imports
#import clusterize
import Final_display_results
import printTable
import run_mr_bump_shelx_parallel
import mrbump_cmd


def mrbump_ensemble_cluster( ensembles, amoptd, clusterID="X" ):
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

    mrBuildClusterDir = os.path.join(amopt.d['mrbump_dir'], "cluster_run" + str(clusterID))
    os.mkdir(mrBuildClusterDir)
    #mrBuild.getMTZInfo(amopt.d['mtz'], mrBuildClusterDir)

    # Reset the queue list
    mrBuild.qList = []
    jobID = 0
    for pdbfile in ensembles:
        mrBuild.mrBuildOnCluster(mrBuildClusterDir, pdbfile, jobID, amoptd )
        jobID = jobID + 1
    mrBuild.monitorQueue()
    
    work_dir = amopt.d['work_dir']
    shutil.rmtree(work_dir + '/fine_cluster_' + str(clusterID))
    # shutil.rmtree(work_dir+'/pre_models')
    for l in os.listdir(work_dir + '/spicker_run'):
        if os.path.splitext(l)[1] == 'pdb':
            os.remove(work_dir + '/spicker_run/' + l)
    os.remove(work_dir + '/spicker_run/rep1.tra1')
    T = printTable.Table()

    T.bumppath = mrBuildClusterDir
    T.cluster = True
    table = T.maketable()
    out = sys.stdout
    T.pprint_table(out, table)

        # cleanup
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
    
    # Monitor the cluster queue to see when all jobs have finished
##End mrbump_ensemble_cluster

def mrbump_ensemble_local( ensembles, amoptd, clusterID="X" ):
    """
    Process the list of ensembles using MrBump on a local machine
    
    Args:
    ensembles -- list of pdb files, each file containing a group of ensembles
    amoptd -- dictionary with job options
    clusterId -- number/id of the cluster that is being processed
    """

    # Queue to hold the jobs we want to run
    queue = multiprocessing.Queue()
    
    # Create all the run scripts and add to the queue
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
    done=0
    timeout=10*60
    timeout=1*60
    killall=False # if we early terminate we check this to see if we kill any remaining jobs
    while done != len(processes):
        
        for process in processes:
            
            if killall:
                process.terminate()
            else:
                # Join process for timeout seconds and if we haven't finished by then
                # move onto the next process
                process.join(timeout)
                
            if not process.is_alive():
                # Finished so see what happened
                if process.exitcode == 0 and amoptd['early_terminate']:
                    # Got a successful completion                    
                    #print "EARLY TERMINATE"
                    killall=True
                done+=1
                
    print results_summary( os.getcwd() )
    
    
#    # At this point all jobs have finished
#    bump_dir = amopt.d['mrbump_dir']
#    # Process results
#    if amopt.d['use_shelxe']:
#        Final_display_results.make_log(bump_dir, os.path.join(amopt.d['work_dir'], 'Final_results.log'))
#        # print '\n\nFinal Results:\n\n'
#        # T=printTable.Table()
#        # T.bumppath = work_dir +'/MRBUMP_cluster'+str(cluster)
#        # T.cluster = False
#        # table = T.maketable()
#        # out = sys.stdout
#        # T.pprint_table(out, table)
#    else:
#        resultslog = open(amopt.d['work_dir'] + os.sep + 'Results.log', "w")
#        print 'getting results from: {}',format( bump_dir )
#        for mrbumplog in os.listdir(bump_dir):
#            if re.search('.log', mrbumplog):
#                # print mrbumplog
#                for line in open(mrbumplog):
#                    if re.search('^(\d)\s*loc0_', line):
#                        if not re.search('method', line):
#                            print line
#                            resultslog.write(line)
#        resultslog.close()
        
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

def  worker( queue, early_terminate=False ):
    """
    Worker process to run MrBump jobs until no more left.
    
    Args:
    queue -- a python Queue object
    early_terminate -- bool - terminate on first success or continue running
    
    Returns:
    0 if molecular replacement worked
    1 if nothing found
    
    We keep looping, removing jobs from the queue until there are no more left.
    """
    
    while True:
        try:
            # We use false to generate an exception when the Queue is empty
            mrb_script = queue.get(False)
            #print "worker {} got script {}".format(os.getpid(), mrb_script )
        except Queue.Empty:
            #print "worker {} got empty queue {}".format(os.getpid(),e)
            break
        
        # Got a script so run
        # Get name from script
        name = os.path.splitext( os.path.split(mrb_script)[1] )[0]
        logfile = name + ".log"
        f = open( logfile, "w")
        retcode = subprocess.call( [ mrb_script ], stdout=f, stderr=subprocess.STDOUT, cwd=None )
        
        # Can we use the retcode to check?
        # REM - is retcode object
        #print "worker {} GOT retcode!".format(os.getpid(), retcode )
        
        # Run directory is name of script with search_ prepended and _mrbump appended
        directory = "search_" + name + "_mrbump"
        
        # Now check the result if early terminate
        if early_terminate:
            if check_success( directory ):
                #print "worker {} GOT SUCCESS!".format(os.getpid() )
                return 0
        
    #print "worker {} FAILED!".format(os.getpid() )  
    return
    
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
    
    rfile = directory + os.sep + 'results/resultsTable.dat' 
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

def results_summary( run_dir ):
    """
    Generate a summary of the MR BUMP results
    
    Args:
    run_dir -- the directory the jobs were run from
    
    Returns:
    results results_file as a string
    """
    
    # Criteria for success with shelx
    SHELXSUCCESS = 25.0

    # Set up the results results_file headers
    resultsTable = []
    
    hline = None

    # Check all results files and add their contents to the resultsTable
    for search in os.listdir(run_dir):
        if os.path.isdir(run_dir+os.sep+search):
            res = run_dir + os.sep + search + '/results/resultsTable.dat'
            if os.path.exists(res):
                firstline=True
                for line in open(res):
                    if firstline:
                        if not hline:
                            hline = line.strip().split()
                        firstline=False
                        continue
                    if  re.search('PHASER', line)  or re.search('MOLREP', line):
                        resultsTable.append(line.split())

    # sort results
    use_shelx=True # For time being assume we always use shelx
    if use_shelx:
        y = hline.index('SHELXE_CC')
    # print resultsTable
    # print y
        resultsTable.sort(key=lambda x: float(x[y]))
        resultsTable.reverse()
        best = resultsTable[0][0]
        prog = resultsTable[0][1]
        if float( resultsTable[0][y]) >=25:
            DIE = True

    else:
        y = hline.index('final_Rfree')

        resultsTable.sort(key=lambda x: float(x[y]))
        best = resultsTable[0][0]
        prog = resultsTable[0][1]


    # Rebuild the path that generated the result
    n= re.sub('loc0_ALL','search', best)
    n = re.sub('UNMOD','mrbump', n)
    Best = run_dir + os.sep + n + '/data/'+re.sub('_UNMOD','', best)+'/unmod/mr/'+prog.lower()+'/refine'

    resultsTable.insert(0, hline)

    # Currently a hack - need a list with a write method
    class WriteList(list):
        def write(self,line):
            self.append(line)

    # Output the results results_file
    T=printTable.Table()
    #out = sys.stdout
    out = WriteList()
    T.pprint_table(out, resultsTable)
    
    
    header = """###########################################################################################
###########################################################################################
##                                                                                       ##
##                                                                                       ##

Overall Summary

"""
    
    r = header
    # now add list to string
    r += "".join(out)
    
    r += '\nBest results so far are in :\n'
    r +=  Best

    footer = """

##                                                                                       ##
##                                                                                       ##
###########################################################################################
###########################################################################################
"""
    
    r+= footer
        
    return r


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
             'domain_all_chains_fasta' : None,
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
        
        dir = "/home/Shared/ample-dev1/examples/toxd-example/ROSETTA_MR_0/MRBUMP_cluster1/search_All_atom_trunc_21.848434_rad_2_mrbump"
        
        self.assertTrue( check_success( dir ) )
        
    def testResultsSummary(self):
        
        dir = "/opt/ample-dev1/python"
        dir = "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_1/MRBUMP_cluster1"
        
        
#        T = printTable.Table()
#        T.bumppath = dir
#        T.cluster = False
#        table = T.maketable()
#        out = sys.stdout
#        T.pprint_table(out, table)

        print results_summary( dir )


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
    