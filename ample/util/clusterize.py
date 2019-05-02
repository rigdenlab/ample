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
import shutil
import time

logger = logging.getLogger(__name__)

class ClusterRun:

    def __init__(self):

        self.qList=[]
        self.runningQueueList=[]
        self.QTYPE=""

        self.modeller = None

        self.runDir=None
        self.logDir=None
        self._scriptFile  = None
        self.debug=True
        
        return

    def cleanUpArrayJob(self, scriptFile=None, logDir=None):
        """Rename all the log files
        Args:
        logDir: directory that the logfiles should end up in
        """
        
        if not scriptFile:
            scriptFile = self._scriptFile
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
                # WARNING: problems with os.rename() call on BADB
                # os.rename docs -->> The operation may fail on some Unix
                # flavors if src and dst are on different filesystems.
                shutil.move(oldLog, newLog)
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
            raise RuntimeError("No jobs found in self.qList!")

        logger.info("Jobs submitted to cluster queue, awaiting their completion...")

        # set a holder for the qlist
        runningList=self.qList
        newRunningList=[]

        while runningList!=[]:
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
        for submitting to a particular queueing system.
    
        Args:
        job_scripts -- the list of scripts to run as the array
        job_name -- the name of the job in the queue (required for LSF array jobs)
        job_dir -- the directory the job will run in
        job_time -- maximum job runtime in minutes (CHECK)
        submit_max_array -- maximum number of array jobs to run concurrently
        submit_queue -- the name of the queue to submit the job to
        submit_qtype -- the type of the queueing system (e.g. SGE)
        
        Returns:
        queue directives as a list of EOL-terminated strings
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
            if nproc and nproc > 1: sh += ['#$ -pe {0} {1}\n'.format(submit_pe_sge, nproc)]
            sh += ['\n']
        elif submit_qtype=="LSF":
            if nproc and submit_pe_lsf: sh += [submit_pe_lsf.format(nproc) + os.linesep]
            if job_time:
                sh += ['#BSUB -W {0}\n'.format(job_time/60)]
            else:
                sh += ['#BSUB -W 4:00\n']
            if nproc and nproc > 1: sh += ['#BSUB -n {0}\n'.format(nproc)]
            if submit_queue: sh += ['#BSUB -q {0}\n'.format(submit_queue)]
            if log_file: sh += ['#BSUB -o {0}\n'.format(log_file)]
            if submit_num_array_jobs:
                assert job_name,"LSF array job requires a job name"
                sh += ['#BSUB -o arrayJob_%I.log\n']
                if submit_max_array:
                    sh += ['#BSUB -J {0}[1-{1}]%{2}\n'.format(job_name, submit_num_array_jobs, submit_max_array)]
                else:
                    sh += ['#BSUB -J {0}[1-{1}]\n'.format(job_name, submit_num_array_jobs)]
            elif job_name: sh += ['#BSUB -J {0}\n'.format(job_name)]       
            sh += ['\n']
        else:
            raise RuntimeError("Unrecognised QTYPE: {0}".format(submit_queue))
        sh += ['\n']
        return sh
    
    def submitJob(self, subScript):
        """
        Submit the job to the queue and return the job number.
        
        Args:
        subScript -- the path to the submission script
        
        Returns:
        job number as a string
        
         We cd to the job directory, submit and then cd back to where we came from
        """
        
        command_line = None
        stdin = None
        if self.QTYPE == "SGE":
            command_line = 'qsub -V %s' % subScript
        elif self.QTYPE == "LSF":
            command_line='bsub'
            stdin = open( subScript, "r")
        else:
            msg = "Unrecognised QTYPE: {0}".format(self.QTYPE)
            raise RuntimeError(msg)

        logger.debug("Submitting job with command: {0}".format(command_line))
        process_args = shlex.split(command_line)
        try:
            p = subprocess.Popen(process_args, stdin = stdin, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        except Exception as e:
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
        while out:
            qNumber = None
            if self.QTYPE == "SGE":
                if "Your job-array" in out:
                    # Array jobs have different form
                    #Your job-array 19094.1-10:1 ("array.script") has been submitted
                    qNumber=int(out.split()[2].split(".")[0])
                    self.qList.append(qNumber)
                elif "Your job" in out:
                    qNumber=int(out.split()[2])
                    self.qList.append(qNumber)
            elif self.QTYPE == "LSF":
                # Job <35339> is submitted to queue <q1h32>.
                if "is submitted to queue" in out:
                    qStr=out.split()[1]
                    qNumber=int(qStr.strip("<>"))
                    self.qList.append(qNumber)                
            if qNumber:
                logger.debug("Submission script {0} submitted to queue as job {1}".format( subScript, qNumber ) )
            out = child_stdout.readline()
        child_stdout.close()
        return str(qNumber)
    
    def submitArrayJob(self,
                       job_scripts,
                       job_name=None,
                       job_time=None,
                       submit_max_array=None,
                       submit_queue=None,
                       submit_qtype=None
                       ):
        """Submit a list of jobs as an array job
        
        Args:
        job_scripts -- the list of scripts to run as the array
        job_name -- the name of the job in the queue (required for LSF array jobs)
        job_dir -- the directory the job will run in
        job_time -- maximum job runtime in minutes (CHECK)
        submit_max_array -- maximum number of array jobs to run concurrently
        submit_queue -- the name of the queue to submit the job to
        submit_qtype -- the type of the queueing system (e.g. SGE)
        
        """
        
        job_dir = os.getcwd()
        
        # Create the list of scripts
        self._scriptFile = os.path.abspath(os.path.join(job_dir,"array.jobs"))
        nJobs = len(job_scripts)
        with open(self._scriptFile,'w') as f:
            for s in job_scripts:
                # Check the scripts are of the correct format - abspath and .sh extension
                if not s.startswith("/") or not s.endswith(".sh"):
                    raise RuntimeError("Scripts for array jobs must be absolute paths with a .sh extension: {0}".format(s))
                f.write(s+"\n")
                
        # Generate the qsub array script
        arrayScript = os.path.abspath(os.path.join(job_dir,"array.script"))
 
        if submit_qtype == "SGE":
            task_env = 'SGE_TASK_ID'
        elif submit_qtype == "LSF":
            task_env = 'LSB_JOBINDEX'
        else:
            raise RuntimeError("Unsupported submission type: {0}".format(submit_qtype))
        
        # Write head of script
        s = "#!/bin/sh\n"
        # Queue directives
        s += "".join(self.queueDirectives(nproc=None,
                                          log_file=None,
                                          job_name=job_name,
                                          job_time=job_time,
                                          submit_max_array=submit_max_array,
                                          submit_num_array_jobs=nJobs,
                                          submit_queue=submit_queue,
                                          submit_qtype=submit_qtype
                                          ))
        # body
        s += """scriptlist={0}

# Extract info on what we need to run
script=`sed -n "${{{1}}}p" $scriptlist`

jobdir=`dirname $script`
jobname=`basename $script .sh`

# cd to jobdir and runit
cd $jobdir

# Run the script
$script
""".format(self._scriptFile, task_env)
        with open(arrayScript,'w') as f: f.write(s)
        self.submitJob(arrayScript)
        return

