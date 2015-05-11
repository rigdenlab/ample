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
import unittest
import time

if not "CCP4" in sorted(os.environ.keys()):
    raise RuntimeError('CCP4 not found')

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
        self.shelxClusterScript="python " + os.path.join(os.environ["CCP4"], "share", "ample", "python", "shelxe_trace.py")

        # Required when a specific python interpreter needs to be invoked in the nodes
        # See ensembleOnCluster
        self.pythonPath = "/home/rmk65/opt/python/python-2.7.2/bin/python"
        self.pythonPath = "ccp4-python"

        self.debug=True

        if os.name == "nt":
            self.pdbsetEXE=os.path.join(os.environ["CCP4"], "bin", "pdbset.exe")
        else:
            self.pdbsetEXE=os.path.join(os.environ["CCP4"], "bin", "pdbset")

        self.logger =  logging.getLogger()
        
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
                self.logger.debug("Moving {0} to {1}".format(oldLog, newLog))
                os.rename(oldLog, newLog ) 
            else:
                self.logger.critical("Cannot find logfile {0} to copy to {1}".format(oldLog, newLog))
        
        return

    def ensembleOnCluster(self, amoptd):
        """ Run the modelling step on a cluster """

        # write out script
        work_dir = amoptd['work_dir']
        script_path = os.path.join( work_dir, "submit_ensemble.sh" )
        qtype=amoptd['submit_qtype']
        self.QTYPE = qtype
        logFile= script_path+".log"
        with open(script_path, "w") as job_script:
            script_header = "#!/bin/sh\n"
            script_header += self.queueDirectives(nProc=1,
                                                 logFile=logFile,
                                                 jobName="ensemble",
                                                 jobTime="3600",
                                                 qtype=qtype,
                                                 queue=amoptd['submit_queue']
                                                 )
            job_script.write(script_header)
    
            # Find path to this directory to get path to python ensemble.py script
            pydir=os.path.abspath( os.path.dirname( __file__ ) )
            ensemble_script = os.path.join(pydir, "ensemble.py")
            job_script.write("{0} {1} {2} {3}\n".format(self.pythonPath, "-u", ensemble_script, amoptd['results_path']))

        # Make executable
        os.chmod(script_path, 0o777)

        # submit
        self.logger.info("Running ensembling on a cluster. Submitting to queue of type: {0}".format(self.QTYPE))
        self.submitJob( subScript=script_path, jobDir=amoptd['work_dir'] )

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
    
    def writeFragmentsSubscript(self, cmd=None, script_path=None, nProc=None, logFile=None, queue=None, qtype=None):
        """ Run the modelling step on a cluster """
        with open(script_path, "w") as job_script:
            script_header="#!/bin/sh\n"
            script_header+=self.queueDirectives(nProc=nProc,
                                               logFile=logFile,
                                               jobName="genFrags",
                                               qtype=qtype,
                                               queue=queue
                                               )
            job_script.write(script_header)
            job_script.write("\n{0}\n".format(cmd))

        # Make executable
        os.chmod(script_path, 0o777)
        return

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

        self.logger.info("Jobs submitted to cluster queue, awaiting their completion...")

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
                self.logger.info("Queue Monitor: %d out of %d jobs remaining in cluster queue..." %  (len(newRunningList),len(self.qList)))
            if len(newRunningList) == 0:
                self.logger.info("Queue Monitor: All jobs complete!")
            runningList=newRunningList
            newRunningList=[]
            if monitor: monitor()
            
        return

    def queueDirectives(self,
                        nProc=None,
                        logFile=None,
                        jobName=None,
                        jobTime=None,
                        queue=None,
                        qtype=None,
                        numArrayJobs=None,
                        maxArrayJobs=None
                        ):
        """
        Create a string suitable for writing out as the header of the submission script
        for submitting to a particular queueing system
        """
        sh = ""
        if qtype=="SGE":
            sh += '#$ -j y\n'
            sh += '#$ -cwd\n'
            sh += '#$ -w e\n'
            sh += '#$ -V\n'
            sh += '#$ -S /bin/bash\n'
            if jobTime: sh += '#$ -l h_rt={0}\n'.format(jobTime)
            if queue: sh += '#$ -q {0}\n'.format(queue)
            if numArrayJobs:
                sh += '#$ -o arrayJob_$TASK_ID.log\n'
                sh += '#$ -t 1-{0}\n'.format(numArrayJobs)
                if maxArrayJobs: sh += '#$ -tc {0}\n'.format(maxArrayJobs)
            else:
                if logFile: sh += '#$ -o {0}\n'.format(logFile)
                if jobName: sh += '#$ -N {0}\n'.format(jobName)
            # jmht hack for morrigan
            if nProc and nProc > 1: sh += '#$ -pe threaded {0}\n'.format( nProc )
            sh += '\n'
        elif qtype=="LSF":
            assert not numArrayJobs,"Array jobs not supported yet for LSF"
            # jmht - hard-wired for hartree wonder
            if nProc and nProc < 16:
                sh += '#BSUB -R "span[ptile={0}]"\n'.format(nProc)
            else:
                sh += '#BSUB -R "span[ptile=16]"\n'
            if jobTime:
                sh += '#BSUB -W {0}\n'.format(jobTime)
            else:
                sh += '#BSUB -W 4:00\n'
            if nProc: sh += '#BSUB -n {0}\n'.format(nProc) 
            if queue: sh += '#BSUB -q {0}\n'.format(queue)
            if logFile: sh += '#BSUB -o {0}\n'.format(logFile)
            if jobName: sh += '#BSUB -J {0}\n'.format(jobName)         
            sh += '\n'
        else:
            raise RuntimeError,"Unrecognised QTYPE: {0}".format(qtype)
        sh += '\n'
        
        return sh
    
    def setupModellingDir(self, RunDir):
        """ A function to create the necessary directories for the modelling step """

        self.runDir=RunDir

        # Create the directories for the submission scripts
        #if not os.path.isdir(os.path.join(RunDir, "models")):
        #    os.mkdir(os.path.join(RunDir, "models"))
        if not os.path.isdir(os.path.join(RunDir, "pre_models")):
            os.mkdir(os.path.join(RunDir, "pre_models"))
        
        self.scriptDir=os.path.join(RunDir, "pre_models", "submit_scripts")
        if not os.path.isdir(self.scriptDir):
            os.mkdir(self.scriptDir)
        
        self.logDir=os.path.join(RunDir, "pre_models", "logs")
        if not os.path.isdir(self.logDir):
            os.mkdir(self.logDir)
        
        return
    
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
        if jobDir:
            os.chdir(jobDir)
        
        command_line=None
        stdin = None
        if self.QTYPE=="SGE":
            command_line='qsub -V %s' % subScript
        elif self.QTYPE=="LSF":
            command_line='bsub'
            stdin = open( subScript, "r")
        else:
            raise RuntimeError,"Unrecognised QTYPE: ".format(self.QTYPE)            

        self.logger.debug("Submitting job with command: {0}".format(command_line))
        process_args = shlex.split(command_line)
        try:
            p = subprocess.Popen(process_args,
                                 stdin = stdin,
                                 stdout = subprocess.PIPE)
        except Exception,e:
            raise RuntimeError,"Error submitting job to queue with commmand: {0}\n{1}".format(command_line,e)

        child_stdout = p.stdout
        # Watch the output for successful termination
        out=child_stdout.readline()

        qNumber=0
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
                self.logger.debug("Submission script {0} submitted to queue as job {1}".format( subScript, qNumber ) )
            out=child_stdout.readline()
        child_stdout.close()
        os.chdir(curDir)
        return str(qNumber)
    
    def submitArrayJob(self,jobScripts,jobDir=None,jobTime=None,queue=None,qtype=None,maxArrayJobs=None):
        """Submit a list of jobs as an SGE array job"""
        
        if qtype != "SGE": raise RuntimeError,"Need to add code for non-SGE array jobs"
        if jobDir is None:
            if self.scriptDir and os.path.isdir(self.scriptDir):
                jobDir=self.scriptDir
            else:
                jobDir=os.getcwd()
        os.chdir(jobDir)
        
        # Create the list of scripts
        self._scriptFile = os.path.abspath(os.path.join(jobDir,"array.jobs"))
        nJobs=len(jobScripts)
        with open(self._scriptFile,'w') as f:
            for s in jobScripts:
                # Check the scripts are of the correct format - abspath and .sh extension
                if not s.startswith("/") or not s.endswith(".sh"):
                    raise RuntimeError,"Scripts for array jobs must be absolute paths with a .sh extension: {0}".format(s)
                f.write(s+"\n")
                
        # Generate the qsub array script
        arrayScript = os.path.abspath(os.path.join(jobDir,"array.script"))
 
        # Write head of script
        s = "#!/bin/sh\n"
        s += self.queueDirectives(nProc=None,
                                  logFile=None,
                                  jobName=None,
                                  jobTime=jobTime,
                                  queue=queue,
                                  qtype=qtype,
                                  numArrayJobs=nJobs,
                                  maxArrayJobs=maxArrayJobs
                                  )
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
        self.submitJob(subScript=arrayScript, jobDir=jobDir)
        return

class Test(unittest.TestCase):
    
    def testSubmit(self):    
        # Test array jobs
        
        # Create run scripts
        jobScripts=[]
        for i in range(10):
            s = """#!/bin/bash
echo "I am script {0}"
""".format(i)
            script=os.path.abspath(os.path.join(os.getcwd(),"script_{0}.sh".format(i)))
            with open(script,'w') as f:
                f.write(s)
            os.chmod(script, 0o777)
            jobScripts.append(script)
        
        c=ClusterRun()
        qtype="SGE"
        c.QTYPE=qtype
        
        c.submitArrayJob(jobScripts,qtype=qtype)
        c.monitorQueue()
        c.cleanUpArrayJob()
        
        return
