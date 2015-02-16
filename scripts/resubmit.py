import os
import shutil
import sys

ample_root="/home/jmht42/ample-dev1"
sys.path.insert(0,os.path.join(ample_root,"python"))
import mrbump_results

#with open("")
pdb_codes=["1UCS"]
root="/volatile/jmht42/testset/percent"

for pdb in pdb_codes:
    mrbd=os.path.join(root,pdb,"ROSETTA_MR_0","MRBUMP")
    os.chdir(mrbd)
    
    res_sum = mrbump_results.ResultsSummary()
    res_sum.extractResults(mrbd)
    ensembles=[]
    for r in res_sum.results:
        search_dir=r['Search_directory']
        ensemble=r['ensemble_name']
        x=os.path.join(search_dir,"results","finished.txt")
        if r['Solution_Type']=='unfinished':
           ensembles.append(ensemble) 
