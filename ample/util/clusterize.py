#!/usr/bin/env python

# Class for submitting ample modelling and MR/Refine/Build components to a cluster
# queuing system.
#
# Ronan Keegan 25/10/2011
#

import logging
import os
import subprocess
import shlex
import time

if not "CCP4" in sorted(os.environ.keys()):
    raise RuntimeError('CCP4 not found')

logger = logging.getLogger(__name__)

class ClusterRun:

    def __init__(self):

        self.qList=[]
        self.runningQueueList=[]
        self.QTYPE=""

        self.modeller = None

        self.runDir=None
        self.logDir=None
        self.scriptDir=None
        self._scriptFile  = None
        self.debug=True
        
        return

    def cleanUpArrayJob(self,scriptFile=None,logDir=None):
        """Rename all the log files
        Args:
        logDir: directory that the logfiles should end up in
        """
        
        if not scriptFile:
            scriptFile=self._scriptFile
        assert os.path.isfile(scriptFile),"Cannot find scriptFile {0}".format(scriptFile)
        
        scriptFiles = []
        with open(scriptFile) as f:
            for line in f:
                scriptFiles.append(line.strip())
        
        for i, line in enumerate(scriptFiles):
            jobDir, script = os.path.split(line)
            jobName = os.path.splitext(script)[0]
            oldLog = "arrayJob_{0}.log".format(i+1)
            if logDir is None:
                # Put log in script directory
                newLog = os.path.join(jobDir, "{0}.log".format(jobName))
            else:
                newLog = os.path.join(logDir, "{0}.log".format(jobName))
                
            if os.path.isfile(oldLog):
                logger.debug("Moving {0} to {1}".format(oldLog, newLog))
                os.rename(oldLog, newLog ) 
            else:
                logger.critical("Cannot find logfile {0} to copy to {1}".format(oldLog, newLog))
        
        return

    # Currently unused
    def XgetJobStatus(self, qNumber):
        """ Check a job status int the cluster queue """

        status=1

        command_line='qstat -j %d' % qNumber

        process_args = shlex.split(command_line)
        p = subprocess.Popen(process_args, stdin = subprocess.PIPE,
                                    stdout = subprocess.PIPE, stderr=subprocess.PIPE)

        (child_stdout, child_stderr, child_stdin) = (p.stdout, p.stderr, p.stdin)

        # Write the keyword input
        child_stdin.close()

        child_stdout.close()

        # Watch the output for successful termination
        err=child_stderr.readline()

        while err:
            #sys.stdout.write(out)
            if self.QTYPE=="SGE":
                if "Following jobs do not exist" in err:
                    status=0
            err=child_stderr.readline()

        child_stderr.close()

        return status
    
    def getRunningJobList(self, user=""):
        """ Check a job status int the cluster queue 

            For LSF output is of form:
JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
35340   jxt15-d RUN   q1h32      ida7c42     ida2a40     *ep 5;done Mar 25 13:23
                                             ida2a40
                                             ida2a40
"""

        if self.QTYPE=="SGE":
            if user == "":
                command_line='qstat'
            else:
                command_line='qstat -u ' + user
        elif self.QTYPE=="LSF":
            if user == "":
                command_line='bjobs'
            else:
                command_line='bjobs -u ' + user

        log_lines=[]
        self.runningQueueList=[]

        process_args = shlex.split(command_line)
        p = subprocess.Popen(process_args, stdout = subprocess.PIPE)

        child_stdout = p.stdout

        # Read the output
        out=child_stdout.readline()

        while out:
            #sys.stdout.write(out)
            log_lines.append( out.strip() )
            out=child_stdout.readline()

        child_stdout.close()

        if log_lines != []:
            log_lines.pop(0)
            # SGE has extra header
            if self.QTYPE=="SGE":
                log_lines.pop(0)
            for i in log_lines:
                self.runningQueueList.append(i.split()[0])
        return

    def monitorQueue(self, user="",monitor=None):
        """ Monitor the Cluster queue to see when all jobs are completed """

        if not len(self.qList):
            raise RuntimeError,"No jobs found in self.qList!"

        logger.info("Jobs submitted to cluster queue, awaiting their completion...")

        # set a holder for the qlist
        runningList=self.qList
        newRunningList=[]

        while runningList!=[]:
            #print "runningList is ",runningList
            time.sleep(60)
            self.getRunningJobList(user)
            for job in runningList:
                if str(job) in self.runningQueueList:
                    newRunningList.append(job)
            if len(runningList) > len(newRunningList):
                logger.info("Queue Monitor: %d out of %d jobs remaining in cluster queue..." %  (len(newRunningList),len(self.qList)))
            if len(newRunningList) == 0:
                logger.info("Queue Monitor: All jobs complete!")
            runningList=newRunningList
            newRunningList=[]
            if monitor: monitor()
            
        return

    def queueDirectives(self,
                        nproc=None,
                        log_file=None,
                        job_name=None,
                        job_time=None,
                        submit_max_array=None,
                        submit_num_array_jobs=None,
                        submit_pe_sge='mpi',
                        submit_pe_lsf='#BSUB -R "span[ptile={0}]"',
                        submit_qtype=None,
                        submit_queue=None,
                        ):
        """
        Create a string suitable for writing out as the header of the submission script
        for submitting to a particular queueing system
        """
        sh = []
        if submit_qtype=="SGE":
            sh += ['#$ -j y\n',
                   '#$ -cwd\n',
                   '#$ -w e\n',
                   '#$ -V\n',
                   '#$ -S /bin/bash\n']
            if job_time: sh += ['#$ -l h_rt={0}\n'.format(job_time)]
            if submit_queue: sh += ['#$ -q {0}\n'.format(submit_queue)]
            if job_name: sh += ['#$ -N {0}\n'.format(job_name)]
            if submit_num_array_jobs:
                sh += ['#$ -o arrayJob_$TASK_ID.log\n']
                sh += ['#$ -t 1-{0}\n'.format(submit_num_array_jobs)]
                if submit_max_array: sh += ['#$ -tc {0}\n'.format(submit_max_array)]
            else:
                if log_file: sh += ['#$ -o {0}\n'.format(log_file)]
            # jmht hack for morrigan
            #if nproc and nproc > 1: sh += ['#$ -pe threaded {0}\n'.format(nproc)]
            if nproc: sh += ['#$ -pe {0} {1}\n'.format(submit_pe_sge, nproc)]
            sh += ['\n']
        elif submit_qtype=="LSF":
            assert not submit_num_array_jobs,"Array jobs not supported yet for LSF"
            if nproc: sh += [submit_pe_lsf.format(nproc) + os.linesep]
            if job_time:
                sh += ['#BSUB -W {0}\n'.format(job_time)]
            else:
                sh += ['#BSUB -W 4:00\n']
            if nproc: sh += ['#BSUB -n {0}\n'.format(nproc)]
            if submit_queue: sh += ['#BSUB -q {0}\n'.format(submit_queue)]
            if log_file: sh += ['#BSUB -o {0}\n'.format(log_file)]
            if job_name: sh += ['#BSUB -J {0}\n'.format(job_name)]       
            sh += ['\n']
        else:
            raise RuntimeError,"Unrecognised QTYPE: {0}".format(submit_queue)
        sh += ['\n']
        return sh
    
    def submitJob(self, subScript=None, jobDir=None):
        """
        Submit the job to the queue and return the job number.
        
        Args:
        subScript -- the path to the submission script
        jobDir -- the directory the job is submitted from - will run in
        
        Returns:
        job number as a string
        
         We cd to the job directory, submit and then cd back to where we came from
        """
        
        curDir=os.getcwd()
        if jobDir: os.chdir(jobDir)
        
        command_line=None
        stdin = None
        if self.QTYPE=="SGE":
            command_line='qsub -V %s' % subScript
        elif self.QTYPE=="LSF":
            command_line='bsub'
            stdin = open( subScript, "r")
        else:
            raise RuntimeError,"Unrecognised QTYPE: ".format(self.QTYPE)            

        logger.debug("Submitting job with command: {0}".format(command_line))
        process_args = shlex.split(command_line)
        try:
            p = subprocess.Popen(process_args,
                                 stdin = stdin,
                                 stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE)
        except Exception,e:
            raise RuntimeError("Error submitting job to queue with commmand: {0}\n{1}".format(command_line,e))

        child_stdout = p.stdout
        child_stderr = p.stderr

        # Check there were no errors
        stderr_str = child_stderr.readline()
        if self.QTYPE=="SGE" and "Unable to run job" in stderr_str:
            raise RuntimeError("Error submitting job to cluster queueing system: {0}".format(stderr_str))

        # Watch the output for successful termination
        out = child_stdout.readline()

        qNumber=0
        err_str = None
        while out:
            qNumber = None
            if self.QTYPE=="SGE":
                if "Your job-array" in out:
                    # Array jobs have different form
                    #Your job-array 19094.1-10:1 ("array.script") has been submitted
                    qNumber=int(out.split()[2].split(".")[0])
                    self.qList.append(qNumber)
                elif "Your job" in out:
                    qNumber=int(out.split()[2])
                    self.qList.append(qNumber)
            elif self.QTYPE=="LSF":
                # Job <35339> is submitted to queue <q1h32>.
                if "is submitted to queue" in out:
                    qStr=out.split()[1]
                    qNumber=int(qStr.strip("<>"))
                    self.qList.append(qNumber)                

            if qNumber:
                logger.debug("Submission script {0} submitted to queue as job {1}".format( subScript, qNumber ) )
            out=child_stdout.readline()
        child_stdout.close()
        os.chdir(curDir)
        return str(qNumber)
    
    def submitArrayJob(self,job_scripts,job_name=None,job_dir=None,job_time=None,queue=None,qtype=None,max_array_jobs=None):
        """Submit a list of jobs as an SGE array job"""
        
        if qtype != "SGE": raise RuntimeError,"Need to add code for non-SGE array jobs"
        if job_dir is None:
            if self.scriptDir and os.path.isdir(self.scriptDir):
                job_dir=self.scriptDir
            else:
                job_dir=os.getcwd()
        os.chdir(job_dir)
        
        # Create the list of scripts
        self._scriptFile = os.path.abspath(os.path.join(job_dir,"array.jobs"))
        nJobs=len(job_scripts)
        with open(self._scriptFile,'w') as f:
            for s in job_scripts:
                # Check the scripts are of the correct format - abspath and .sh extension
                if not s.startswith("/") or not s.endswith(".sh"):
                    raise RuntimeError,"Scripts for array jobs must be absolute paths with a .sh extension: {0}".format(s)
                f.write(s+"\n")
                
        # Generate the qsub array script
        arrayScript = os.path.abspath(os.path.join(job_dir,"array.script"))
 
        # Write head of script
        s = "#!/bin/sh\n"
        s += "".join(self.queueDirectives(nproc=None,
                                  log_file=None,
                                  job_name=job_name,
                                  job_time=job_time,
                                  queue=queue,
                                  qtype=qtype,
                                  num_array_jobs=nJobs,
                                  max_array_jobs=max_array_jobs
                                  ))
        # Command to run 
        s += "scriptlist={0}\n".format(self._scriptFile)

        # Add on the rest of the script - need to do in two bits or the stuff in here gets interpreted by format
        s += """
# Extract info on what we need to run
script=`sed -n "${SGE_TASK_ID}p" $scriptlist`

jobdir=`dirname $script`
jobname=`basename $script .sh`

# cd to jobdir and runit
cd $jobdir

# Run the script
$script
"""
        with open(arrayScript,'w') as f: f.write(s)
        self.submitJob(subScript=arrayScript, jobDir=job_dir)
        return

