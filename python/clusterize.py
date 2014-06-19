#!/usr/bin/env python

# Class for submitting ample modelling and MR/Refine/Build components to a cluster
# queuing system.
#
# Ronan Keegan 25/10/2011
#

import logging
import os
import sys
import subprocess
import shlex
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
        self.pythonPath = "python"
        self.pythonPath = "/home/rmk65/opt/python/python-2.7.2/bin/python"
        

        self.debug=True

        if os.name == "nt":
            self.pdbsetEXE=os.path.join(os.environ["CCP4"], "bin", "pdbset.exe")
        else:
            self.pdbsetEXE=os.path.join(os.environ["CCP4"], "bin", "pdbset")

        self.logger =  logging.getLogger()
        
        return

    def cleanUpArrayJob(self,logDir=None):
        """Rename all the log files
        Args:
        logDir: directory that the logfiles should end up in
        """
        
        assert os.path.isfile(self._scriptFile)
        
        scriptFiles = []
        with open( self._scriptFile ) as f:
            for line in f:
                scriptFiles.append(line.strip())
        
        for i, line in enumerate( scriptFiles ):
            jobDir, script = os.path.split( line )
            jobName = os.path.splitext(script)[0]
            oldLog = "arrayJob_{0}.log".format( i+1 )
            if logDir is None:
                # Put log in script directory
                newLog = os.path.join( jobDir, "{0}.log".format( jobName ) )
            else:
                newLog = os.path.join( logDir, "{0}.log".format( jobName ) )
                
            self.logger.debug( "Moving {0} to {1}".format( oldLog, newLog ) )
            os.rename( oldLog, newLog ) 
        
        return

    def ensembleOnCluster(self, amoptd):
        """ Run the modelling step on a cluster """

        # write out script
        work_dir = amoptd['work_dir']
        script_path = os.path.join( work_dir, "submit_ensemble.sh" )
        job_script = open(script_path, "w")

        self.QTYPE = amoptd['submit_qtype']
        logFile= script_path+".log"
        script_header = self.subScriptHeader( nProc=1, logFile=logFile, jobName="ensemble", jobTime="1:00")
        job_script.write( script_header )

        # Find path to this directory to get path to python ensemble.py script
        pydir=os.path.abspath( os.path.dirname( __file__ ) )
        ensemble_script = os.path.join( pydir, "ensemble.py" )

        job_script.write("{0} {1} {2}\n".format( self.pythonPath, ensemble_script, amoptd['results_path'] ) )
        job_script.close()

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
    
    def generateFragmentsOnCluster(self, cmd=None, fragmentsDir=None, nProc=None, logFile=None ):
        """ Run the modelling step on a cluster """

        # write out script
        script_path = os.path.join( fragmentsDir, "submit_fragments.sh" )
        job_script = open(script_path, "w")

        script_header = self.subScriptHeader( nProc=nProc, logFile=logFile, jobName="genFrags")
        job_script.write( script_header )
        job_script.write("\n{0}\n".format( cmd ) )
        job_script.close()

        # Make executable
        os.chmod(script_path, 0o777)
  
        self.submitJob( subScript=script_path, jobDir=fragmentsDir )

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

    def NMRmodelOnCluster(self, RunDir, proc, jobNumber, ROSETTA_PATH, ROSETTA_DB, FASTA, frags_3_mers, frags_9_mers, ideal_homolog,  ALI, seed, MR_ROSETTA ):
        """ Farm out the modelling step on a cluster (SGE) """
        # Set the scriptFile number according to the job number
        if jobNumber<10:
            fileNumber="0000000" + str(jobNumber)
        elif jobNumber>=10 and jobNumber<100:
            fileNumber="000000" + str(jobNumber)
        elif jobNumber>=100 and jobNumber<1000:
            fileNumber="00000" + str(jobNumber)
        elif jobNumber>=1000 and jobNumber<10000:
            fileNumber="0000" + str(jobNumber)
        elif jobNumber>=10000 and jobNumber<100000:
            fileNumber="000" + str(jobNumber)
        elif jobNumber>=100000 and jobNumber<1000000:
            fileNumber="00" + str(jobNumber)
        else:
            sys.stdout.write("No. of Models exceeds program limits (Max=999999)\n")
            sys.exit()
        preModelDir=os.path.join(RunDir, "pre_models", "model_" + str(jobNumber))

        if os.path.isdir(preModelDir) == False:
            os.mkdir(preModelDir)

        n = os.path.split(ideal_homolog)
        n ='S_' +n[-1].rstrip('.pdb')+'_0001.pdb'


        PDBInFile     = os.path.join(preModelDir, n)
        PDBSetOutFile = os.path.join(preModelDir, "pdbsetOut_" + str(jobNumber) + ".pdb")
        PDBScwrlFile  = os.path.join(preModelDir, "scwrlOut_" + str(jobNumber) + ".pdb")
        SEQFile       = os.path.join(preModelDir, "S_" + fileNumber + ".seq")
        PDBOutFile    = os.path.join(RunDir, "models", "1_S_" + fileNumber + ".pdb")

        # Create a cluster submission script for this modelling job
        jobName="model_" + str(proc) + "_" + str(seed)
        sub_script=os.path.join(RunDir, "pre_models", "submit_scripts", "job_" + jobName + ".sub")

        #self.jobLogsList.append(os.path.join(runDir, "pre_models", "logs", jobName + '.log'))

        logFile = os.path.join(RunDir, "pre_models", "logs", jobName + '.log')
        scriptFile=open(sub_script, "w")
        
        script_header = self.subScriptHeader(logFile=logFile, jobName=jobName)
        scriptFile.write(script_header)

        #scriptFile.write("export CCP4_SCR=$TMPDIR\n\n")

        # jmht - this needs to go in the rosetta object
        scriptFile.write('cd '+ os.path.join(RunDir, "pre_models", "model_" + str(jobNumber)) +'\n\n'+
             MR_ROSETTA +' \\\n'+
             '-database '+ROSETTA_DB+' \\\n'+
             '-MR:mode cm \\\n'+
             '-in:scriptFile:extended_pose 1 \\\n'+
             '-in:scriptFile:fasta '+FASTA+' \\\n'+
             '-in:scriptFile:alignment '+ALI+' \\\n'+
             '-in:scriptFile:template_pdb '+ideal_homolog+' \\\n'+
             '-loops:frag_sizes 9 3 1 \\\n'+
             '-loops:frag_files '+frags_9_mers+' '+frags_3_mers+' none \\\n'+
             '-loops:random_order \\\n'+
             '-loops:random_grow_loops_by 5 \\\n'+
             '-loops:extended \\\n'+
             '-loops:remodel quick_ccd \\\n'+
             '-loops:relax relax \\\n'+
             '-relax:default_repeats 4 \\\n'+
             '-relax:jump_move true    \\\n'+
             '-cm:aln_format grishin \\\n'+
             '-MR:max_gaplength_to_model 8 \\\n'+
             '-nstruct 1  \\\n'+
             '-ignore_unrecognized_res \\\n'+
             '-overwrite \n\n')

        scriptFile.write("pushd " + os.path.join(preModelDir) + "\n\n" +

        self.pdbsetEXE + " xyzin " + PDBInFile + " xyzout " + PDBSetOutFile + "<<eof\n" +
        "sequence single\n" +
        "eof\n\n" +

        "tail -n +2 SEQUENCE | sed s'/ //g' >> " + SEQFile + "\n" +
        "popd\n\n"  )
        if self.modeller.use_scwrl:
            scriptFile.write( self.modeller.scwrl_exe + " -i " + PDBInFile + " -o " + PDBScwrlFile + " -s " + SEQFile + "\n\n" +
            "head -n -1 " + PDBScwrlFile + " >> " + PDBOutFile + "\n" +
             "\n")
        else:
            scriptFile.write('cp ' + PDBInFile + ' ' +  PDBOutFile + "\n" )

        # Clean up non-essential files unless we are debugging
        if self.debug == False:
            scriptFile.write("rm " + PDBSetOutFile + "\n" +
            "rm " + os.path.join(preModelDir, "SEQUENCE") + "\n" +
            "rm " + PDBScwrlFile + "\n\n")

        scriptFile.close()

        jobDir = os.path.join(self.modeller.work_dir, "pre_models", "submit_scripts")
        
        self.submitJob(subScript=sub_script, jobDir=jobDir)
        
        return
    
    def modelOnCluster(self,modeller,amoptd):

        self.logger.info('Submitting {0} jobs to a queueing system of type: {1}\n'.format( amoptd['nmodels'],amoptd['submit_qtype'] ) )
        if amoptd['submit_array']:
            self.logger.info('Submitting SGE array job')

        self.QTYPE = amoptd['submit_qtype']
        
        self.modeller=modeller
        self.modeller.generate_seeds( amoptd['nmodels'] )
        if self.modeller.transmembrane:
            self.modeller.generate_tm_predict()    
        
        self.setupModellingDir( self.modeller.models_dir )
        
        jobScripts = []
        # loop over the number of models and submit a job to the cluster
        for i in range( amoptd['nmodels'] ):
            jobNumber=i+1
            jobScript, jobDir = self.writeModelScript(jobNumber)
            if amoptd['submit_array']:
                jobScripts.append(jobScript)
            else:
                self.submitJob(subScript=jobScript, jobDir=jobDir)
        
        # For array jobs we submit as one
        if amoptd['submit_array']:
            self.submitArrayJob(jobScripts, jobDir=amoptd['work_dir'])
            
        # Monitor the cluster queue to see when all jobs have finished
        self.monitorQueue()
        
        if amoptd['submit_array']:
            self.cleanUpArrayJob(logDir=self.logDir)
        
        return

    def monitorQueue(self, user=""):
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
            
        return
    
    def setupModellingDir(self, RunDir):
        """ A function to create the necessary directories for the modelling step """

        self.runDir=RunDir

        # Create the directories for the submission scripts
        if not os.path.isdir(os.path.join(RunDir, "models")):
            os.mkdir(os.path.join(RunDir, "models"))
        if not os.path.isdir(os.path.join(RunDir, "pre_models")):
            os.mkdir(os.path.join(RunDir, "pre_models"))
        
        self.scriptDir=os.path.join(RunDir, "pre_models", "submit_scripts")
        if not os.path.isdir(self.scriptDir):
            os.mkdir(self.scriptDir)
        
        self.logDir=os.path.join(RunDir, "pre_models", "logs")
        if not os.path.isdir(self.logDir):
            os.mkdir(self.logDir)
        
        return

    def subScriptHeader(self, nProc=None, logFile=None, jobName=None, jobTime=None):
        """
        Create a string suitable for writing out as the header of the submission script
        for submitting to a particular queueing system
        """
        
        sh = ""
        if self.QTYPE=="SGE":
            sh += '#!/bin/sh\n'
            sh += '#$ -j y\n'
            sh += '#$ -cwd\n'
            sh += '#$ -w e\n'
            sh += '#$ -V\n'
            sh += '#$ -S /bin/bash\n'
            sh += '#$ -o {0}\n'.format(logFile) 
            sh += '#$ -N {0}\n\n'.format(jobName)
            # jmht hack for morrigan
            if nProc and nProc > 1:
                sh += '#$ -pe threaded {0}\n\n'.format( nProc )
        elif self.QTYPE=="LSF":
            sh += '#!/bin/sh\n'
            # jmht - hard-wired for hartree wonder
            if nProc and nProc < 16:
                sh += '#BSUB -R "span[ptile={0}]"\n'.format(nProc)
            else:
                sh += '#BSUB -R "span[ptile=16]"\n'
            if jobTime:
                sh += '#BSUB -W {0}\n'.format(jobTime)
            else:
                sh += '#BSUB -W 4:00\n'
            if nProc:
                sh += '#BSUB -n {0}\n'.format(nProc) 
            sh += '#BSUB -o {0}\n'.format(logFile) 
            sh += '#BSUB -J {0}\n\n'.format(jobName)         
        else:
            raise RuntimeError,"Unrecognised QTYPE: {0}".format(self.QTYPE)
        
        # Make sure the CCP4 scratch directory is available
        sh += '[[ ! -d $CCP4_SCR ]] && mkdir $CCP4_SCR\n\n'
        
        return sh+'\n\n'
    
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

        process_args = shlex.split(command_line)
        p = subprocess.Popen(process_args, stdin = stdin,
                                      stdout = subprocess.PIPE)

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
    
    def submitArrayJob(self,jobScripts,jobDir=None):
        """Submit a list of jobs as an SGE array job"""
        
        if self.QTYPE != "SGE":
            raise RuntimeError,"Need to add code for non-SGE array jobs"

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
        s = """#!/bin/bash

# Set up SGE variables
#$ -j y
#$ -cwd
#$ -w e
#$ -V
#$ -o arrayJob_$TASK_ID.log
#$ -t 1-{0}
#
# Ignore for now as we always run single processor jobs
##$ -pe smp 16

scriptlist={1}
""".format(nJobs,self._scriptFile)

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
        
        with open(arrayScript,'w') as f:
            f.write(s)
        
        # submit the array script
        self.submitJob( subScript=arrayScript, jobDir=jobDir )
        
        return

    def writeModelScript(self, jobNumber):
        """ Farm out the modelling step on a cluster (SGE) """
        
        nProc=1

        # Set the scriptFile number according to the job number
        if jobNumber<10:
            fileNumber="0000000" + str(jobNumber)
        elif jobNumber>=10 and jobNumber<100:
            fileNumber="000000" + str(jobNumber)
        elif jobNumber>=100 and jobNumber<1000:
            fileNumber="00000" + str(jobNumber)
        
        elif jobNumber>=1000 and jobNumber<10000:
            fileNumber="0000" + str(jobNumber)
        elif jobNumber>=10000 and jobNumber<100000:
            fileNumber="000" + str(jobNumber)
        elif jobNumber>=100000 and jobNumber<1000000:
            fileNumber="00" + str(jobNumber)
        else:
            sys.stdout.write("No. of Models exceeds program limits (Max=999999)\n")
            sys.exit()

        preModelDir=os.path.join(self.runDir, "pre_models", "model_" + str(jobNumber))

        if not os.path.isdir(preModelDir):
            os.mkdir(preModelDir)

        PDBInFile     = os.path.join(preModelDir, "S_00000001.pdb")
        PDBSetOutFile = os.path.join(preModelDir, "pdbsetOut_" + str(jobNumber) + ".pdb")
        PDBScwrlFile  = os.path.join(preModelDir, "scwrlOut_" + str(jobNumber) + ".pdb")
        SEQFile       = os.path.join(preModelDir, "S_" + fileNumber + ".seq")
        PDBOutFile    = os.path.join(self.modeller.models_dir, "1_S_" + fileNumber + ".pdb")

        # Get the seed for this job
        seed = self.modeller.seeds[jobNumber-1]

        # Create a cluster submission script for this modelling job
        jobName="model_" + str(nProc) + "_" + str(seed)
        scriptPath=os.path.join(self.runDir, "pre_models", "submit_scripts", "job_" + jobName + ".sh")
        logFile = os.path.join(self.runDir, "pre_models", "logs", jobName + '.log')
        
        with open(scriptPath, "w") as scriptFile:
            # Modelling always run on single processor
            script_header = self.subScriptHeader( nProc=1, logFile=logFile, jobName=jobName)
            scriptFile.write(script_header+"\n\n")
            #scriptFile.write("export CCP4_SCR=$TMPDIR\n\n")
    
            # Build up the rosetta command
            nstruct=1 # 1 structure
            rcmd = self.modeller.modelling_cmd( preModelDir, nstruct, seed )
            cmdstr = " ".join(rcmd) + "\n\n"
            scriptFile.write( cmdstr )
    
            scriptFile.write("pushd " + os.path.join(preModelDir) + "\n\n" +
    
            self.pdbsetEXE + " xyzin " + PDBInFile + " xyzout " + PDBSetOutFile + "<<eof\n" +
            "sequence single\n" +
            "eof\n\n" +
    
            "tail -n +2 SEQUENCE | sed s'/ //g' >> " + SEQFile + "\n" +
            "popd\n\n"  )
            if self.modeller.use_scwrl:
                scriptFile.write( self.modeller.scwrl_exe + " -i " + PDBInFile + " -o " + PDBScwrlFile + " -s " + SEQFile + "\n\n" +
                "head -n -1 " + PDBScwrlFile + " >> " + PDBOutFile + "\n" +
                 "\n" )
            else:
                scriptFile.write('cp ' + PDBInFile + ' ' +  PDBOutFile + "\n" )
    
            # Clean up non-essential files unless we are debugging
            if not self.debug:
                scriptFile.write("rm " + PDBSetOutFile + "\n" +
                "rm " + os.path.join(preModelDir, "SEQUENCE") + "\n" +
                "rm " + PDBScwrlFile + "\n\n")

        # Make executble
        os.chmod( scriptPath, 0o777)
        jobDir = os.path.join(self.runDir, "pre_models", "submit_scripts")
        
        return scriptPath,jobDir

if __name__ == "__main__":
    
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
    c.QTYPE="SGE"
    
    c.submitArrayJob(jobScripts)
    c.monitorQueue()
    c.cleanUpArrayJob()
