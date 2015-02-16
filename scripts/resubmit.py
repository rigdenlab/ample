import os
import shutil
import sys

ample_root="/opt/ample-dev1"
sys.path.insert(0,os.path.join(ample_root,"python"))
import mrbump_results

#with open("")
pdb_codes=["1UCS"]
root="/volatile/jmht42/testset/percent"

for pdb in pdb_codes:
    mrbd=os.path.join(root,pdb,"ROSETTA_MR_0","MRBUMP")
    os.chdir(mrbd)
    
    for res in ResultsSummary().extractResults(mrbd).results:
        print "GOT ",r['JobDirectory'],r['ensemble_name']
