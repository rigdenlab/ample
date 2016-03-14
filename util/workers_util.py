'''
Created on Feb 28, 2013

@author: jmht
'''

# python imports
import glob
import logging
import multiprocessing
import os
import time

# our imports
from ample.python import clusterize
from ample.util import ample_util
from ample.util import worker

class JobServer(object):
    def __init__(self):
        self.inqueue = None
        self.outqueue = None
        self.logger = logging.getLogger()
        self.logger.info("Running jobs on a local machine")
    
    def setJobs(self, jobs):
        """Add the list of jobs we are to run"""
        if self.inqueue: raise RuntimeError,"NOT THOUGHT ABOUT MULTIPLE INVOCATIONS!"
        
        # Queue to hold the jobs we want to run
        queue = multiprocessing.Queue()
        
        # Add jobs to the inqueue
        #logger.info("Generating MRBUMP runscripts in: {0}".format( os.getcwd() ) )
        for job in jobs:
            if not os.path.isfile(job): raise RuntimeError,"JobServer cannot find job: {0}".format(job)
            queue.put(job)
            
        self.inqueue = queue
        #self.inqueue.close()
        # We can't call inqueue.close() even though we've finished as it all goes horribly wrong
        # We sleep to allow enough time for the objects to be picked and put on the queue
        time.sleep(2)
        
        return
    
    def start(self, nproc=None, early_terminate=False, check_success=None, monitor=None, chdir=False):
        
        assert nproc != None

        # Now start the jobs
        processes = []
        for i in range(nproc):
            process = multiprocessing.Process(target=worker.worker, args=(self.inqueue,
                                                                          early_terminate,
                                                                          check_success,
                                                                          chdir))
            process.start()
            processes.append(process)
        
        # Loop through the processes checking if any are done
        timeout=1*60
        
        # Broken on OSX 
        #qsize=self.outqueue.qsize()
        
        if monitor: monitor()
        
        success=True
        while len(processes):
            for i, process in enumerate(processes):
                # Join process for timeout seconds and if we haven't finished by then move onto the next process
                process.join(timeout)
                if not process.is_alive():
                    self.logger.debug("Checking completed process {0} with exitcode {1}".format(process,process.exitcode))
                    
                    # Set failed if any job failed
                    if process.exitcode != 0:
                        self.logger.critical("Process {0} failed with exitcode {1}".format(process,process.exitcode))
                        success=False
                    
                    # Finished so see what happened
                    if process.exitcode == 0 and early_terminate:
                        if not self.inqueue.empty():
                            print "Process {0} was successful so removing remaining jobs from inqueue".format(process.name) 
                            self.logger.info( "Process {0} was successful so removing remaining jobs from inqueue".format(process.name) )
                            # Remove all remaining processes from the inqueue. We do this rather than terminate the processes
                            # as terminating leaves the MRBUMP processes running. This way we hang around until all our
                            # running processes have finished
                            while not self.inqueue.empty():
                                job = self.inqueue.get()
                                self.logger.debug( "Removed job [{0}] from inqueue".format(job))
                        else:
                            print "Got empty queue - all jobs done"
                            
                    # Remove from processes to check
                    del processes[i]
                    
                # Run the monitor function
                if monitor: monitor()
                            
        # need to wait here as sometimes it takes a while for the results files to get written
        time.sleep(3)        
        return success

def run_scripts(job_scripts,
                monitor=None,
                check_success=None,
                early_terminate=None,
                chdir=False,
                nproc=None,
                job_time=None,
                job_name=None,
                submit_cluster=None,
                submit_qtype=None,
                submit_queue=None,
                submit_array=None,
                submit_max_array=None):
    if submit_cluster:
        return run_scripts_cluster(job_scripts,
                                   monitor=monitor,
                                   job_time=job_time,
                                   job_name=job_name,
                                   submit_cluster=submit_cluster,
                                   submit_qtype=submit_qtype,
                                   submit_queue=submit_queue,
                                   submit_array=submit_array,
                                   submit_max_array=submit_max_array,
                                   )
    else:
        return run_scripts_serial(job_scripts,
                                  nproc=nproc,
                                  monitor=monitor,
                                  early_terminate=early_terminate,
                                  check_success=check_success,
                                  chdir=chdir
                                  )

def run_scripts_cluster(job_scripts,
                        monitor=None,
                        job_time=None,
                        job_name=None,
                        submit_cluster=None,
                        submit_qtype=None,
                        submit_queue=None,
                        submit_array=None,
                        submit_max_array=None,
                        nproc=None):
    logger = logging.getLogger()
    logger.info("Running jobs on a cluster")
    cluster_run = clusterize.ClusterRun()
    qtype = submit_qtype
    cluster_run.QTYPE = submit_qtype
    if submit_array and len(job_scripts) > 1:
        cluster_run.submitArrayJob(job_scripts,
                                   job_time=job_time,
                                   job_name=job_name,
                                   qtype=qtype,
                                   queue=submit_queue,
                                   max_array_jobs=submit_max_array
                                   )
    else:
        for script in job_scripts:
            name=os.path.splitext(os.path.basename(script))[0]
            logfile="{0}.log".format(name)
            if not job_name: job_name = name
            if nproc is None: nproc = 1
            with open(script) as f: lines = f.readlines()
            slines = clusterize.ClusterRun().queueDirectives(nproc=nproc,
                                                             job_name=job_name,
                                                             job_time=job_time,
                                                             log_file=logfile,
                                                             queue=submit_queue,
                                                             qtype=submit_qtype)
            # We add the queue directives after the first line of the script
            with open(script,'w') as f: f.writelines("".join([lines[0]] + slines + lines[1:]))
            os.chmod(script, 0o777)
            cluster_run.submitJob(subScript=script)

    # Monitor the cluster queue to see when all jobs have finished
    cluster_run.monitorQueue(monitor=monitor)
    
    # Rename scripts for array jobs
    if submit_array and len(job_scripts) > 1: cluster_run.cleanUpArrayJob()
    return True

def run_scripts_serial(job_scripts,
                       nproc=None,
                       monitor=None,
                       early_terminate=None,
                       check_success=None,
                       chdir=False):
    success=False
    if len(job_scripts) > 1:
        # Don't need early terminate - check_success if it exists states what's happening
        js = JobServer()
        js.setJobs(job_scripts)
        success = js.start(nproc=nproc,
                           early_terminate=bool(early_terminate),
                           check_success=check_success,
                           monitor=monitor,
                           chdir=chdir)
    else:
        script=job_scripts[0]
        name=os.path.splitext(os.path.basename(script))[0]
        logfile="{0}.log".format(name)
        dir=os.path.dirname(script)
        if chdir: os.chdir(dir)
        rtn = ample_util.run_command([script], logfile=logfile)
        if rtn == 0: success = True
    return success

# Need this defined outside of the test or it can't be pickled on Windoze
def _check_success_test( job ):
    jobname = os.path.splitext(os.path.basename(job))[0]
    if jobname == "job_2": return True
    return False
 
