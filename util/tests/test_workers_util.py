"""Test functions for util.workers_util"""

import glob
import os
import tempfile
import unittest
from ample.util import workers_util

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
        cls.ample_dir = os.sep.join( paths[ : -2 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"testing")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
     
    def makeJob(self, name):
        
        script = """#!/usr/bin/python
import sys,time
print "I am job: {0}"
time.sleep( 3 )
sys.exit(0)
""".format( name )
    
        f = tempfile.NamedTemporaryFile("w+b", prefix=name, suffix="py", delete=False)
        f.write(script)
        f.close()
        os.chmod(f.name, 0o77)
         
        return f.name
 
    def test_jobServer(self):
        jobs = []
        for j in range(15):
            j = os.path.abspath("job_{0}".format(j))
            jobs.append(self.makeJob(j))
             
        js = workers_util.JobServer()
        js.setJobs(jobs)
        js.start(nproc=2, early_terminate=True, check_success=workers_util._check_success_test)
         
        # Cleanup
        for j in jobs: os.unlink(j)
        for l in glob.glob("job_*.log"): os.unlink(l)
        pass
    
if __name__ == "__main__":
    unittest.main()
