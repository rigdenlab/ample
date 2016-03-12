
import os
import unittest
from ample.python import clusterize

class Test(unittest.TestCase):
    
    def test_submit(self):    
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
        
        c=clusterize.ClusterRun()
        qtype="SGE"
        c.QTYPE=qtype
        
        c.submitArrayJob(jobScripts,qtype=qtype)
        c.monitorQueue()
        c.cleanUpArrayJob()
        
if __name__ == "__main__":
    unittest.main()