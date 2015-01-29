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
import unittest

# our imports
import worker

class JobServer(object):
    
    def __init__(self):
        
        
        self.inqueue = None
        self.outqueue = None
        self.logger = logging.getLogger()
        self.logger.info("Running jobs on a local machine")
        
    
    def setJobs(self, jobs):
        """Add the list of jobs we are to run"""

        if self.inqueue:
            raise RuntimeError,"NOT THOUGHT ABOUT MULTIPLE INVOCATIONS!"
        
        # Queue to hold the jobs we want to run
        queue = multiprocessing.Queue()
        
        # Add jobs to the inqueue
        #logger.info("Generating MRBUMP runscripts in: {0}".format( os.getcwd() ) )
        for job in jobs:
            queue.put( job )
            
        self.inqueue = queue
        #self.inqueue.close()
        # We can't call inqueue.close() even though we've finished as it all goes horribly wrong
        # We sleep to allow enough time for the objects to be picked and put on the queue
        time.sleep( 2 )
        
        return
    
    def start( self, nproc=None, early_terminate=False, check_success=None ):
        
        assert nproc != None

        # Now start the jobs
        processes = []
        for i in range( nproc ):
            process = multiprocessing.Process(target=worker.worker, args=( self.inqueue,
                                                                           early_terminate,
                                                                           check_success ) )
            process.start()
            processes.append(process)
        
        # Loop through the processes checking if any are done
        timeout=10*60
        timeout=1*60
        
        # Broken on OSX 
        #qsize=self.outqueue.qsize()
        while len(processes):
            
            for i, process in enumerate(processes):
                
                # Join process for timeout seconds and if we haven't finished by then
                # move onto the next process
                process.join(timeout)
                
                if not process.is_alive():
                    print "Checking completed process {0} with exitcode {1}".format(process,process.exitcode)
                    # Remove from processes to check
                    del processes[i]
                    
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
                                self.logger.debug( "Removed job [{0}] from inqueue".format(job) )
                                print "Removed job [{0}] from inqueue".format(job)
                        else:
                            print "Got empty queue - all jobs done"
                            
        # need to wait here as sometimes it takes a while for the results files to get written
        time.sleep( 3 )
        return

# Need this defined outside of the test or it can't be pickled on Windoze
def _check_success_test( job ):
    import os
    jobname = os.path.splitext( os.path.basename( job ) )[0]
    if jobname == "job_2":
        return True
    return False
 
class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -1 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        return
     
    def makeJob(self, name):
         
        script = """#!/usr/bin/python
import sys,time
print "I am job: {0}"
time.sleep( 3 )
sys.exit(0)
""".format( name )
 
        f = open( name, 'w' )
        f.write( script )
        f.close()
        os.chmod( name, 0o777)
         
        return name
 
    def testJobServer(self):
         
        #print "running in ",os.getcwd()
        jobs = []
        for j in range( 15 ):
            j = os.path.abspath("job_{0}.py".format( j ))
            jobs.append( self.makeJob(j))
             
        js = JobServer()
        js.setJobs( jobs )
        js.start( nproc=2, early_terminate=True, check_success=_check_success_test )
         
        # Cleanup
        for j in jobs:
            os.unlink(j)
        for l in glob.glob("job_*.log"):
            os.unlink(l)
         
        pass
 
def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testJobServer'))
    return suite
     
#
# Run unit tests
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(testSuite())
 
       
