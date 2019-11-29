import os
import unittest

from ample.util import ample_util, clusterize


def on_cluster():
    try:
        retcode = ample_util.run_command(["qstat"])
    except:
        retcode = -1
    return True if retcode == 0 else False


class Test(unittest.TestCase):
    @unittest.skipUnless(on_cluster(), "not on SGE cluster")
    def test_submit_SGE(self):
        jobScripts = []
        for i in range(10):
            s = """#!/bin/bash
echo "I am script {0}"
""".format(
                i
            )
            script = os.path.abspath(os.path.join(os.getcwd(), "script_{0}.sh".format(i)))
            with open(script, 'w') as f:
                f.write(s)
            os.chmod(script, 0o777)
            jobScripts.append(script)

        c = clusterize.ClusterRun()
        qtype = "SGE"
        c.QTYPE = qtype

        c.submitArrayJob(jobScripts, submit_qtype=qtype)
        c.monitorQueue()
        c.cleanUpArrayJob()


if __name__ == "__main__":
    unittest.main()
