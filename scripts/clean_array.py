import glob
import os
import shutil
import sys

ample_root="/home/jmht/ample-dev1"
sys.path.insert(0,os.path.join(ample_root,"python"))

import clusterize

pdb_codes=["1MIX", "1P9G", "1UCS", "1XKR", "2BL2", "2EFR", "2FM9", "2JKU", "2QIH", "2QSK", "2UUI", "2XFD", "2YKT", "3CI9", "3CVF", "3GD8", "3GHF", "3HAP", "3HFE"]
#pdb_codes=["1MIX"]
root="/data2/jmht/testset/thresh"

for pdb in pdb_codes:
    mrbd=os.path.join(root,pdb,"AMPLE_0","MRBUMP")
    os.chdir(mrbd)
    jobsFile = os.path.abspath(os.path.join(mrbd,"array.jobs"))
    #arrayScript = os.path.abspath(os.path.join(mrbd,"array.script"))
    if os.path.isfile(jobsFile): clusterize.ClusterRun().cleanUpArrayJob(jobsFile) 

