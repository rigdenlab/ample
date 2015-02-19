import glob
import os
import shutil
import sys

ample_root="/home/jmht42/ample-dev1"
sys.path.insert(0,os.path.join(ample_root,"python"))

import ample_util
import mrbump_results

def prep_array(jobScripts,jobDir):
    os.chdir(jobDir)
    
    # Create the list of scripts
    scriptFile = os.path.abspath(os.path.join(jobDir,"resub.jobs"))
    nJobs=len(jobScripts)
    with open(scriptFile,'w') as f:
        for s in jobScripts:
            # Check the scripts are of the correct format - abspath and .sh extension
            if not s.startswith("/") or not s.endswith(".sh"):
                raise RuntimeError,"Scripts for array jobs must be absolute paths with a .sh extension: {0}".format(s)
            f.write(s+"\n")
            
    # Generate the qsub array script
    arrayScript = os.path.abspath(os.path.join(jobDir,"resub.script"))
    
    # Write head of script
    s = """#!/bin/bash
# Set up SGE variables
#$ -j y
#$ -cwd
#$ -w e
#$ -V
#$ -l h_rt=172800
#$ -o arrayJob_$TASK_ID.log
#$ -t 1-{0}
#$ -S /bin/bash
#
# Ignore for now as we always run single processor jobs
##$ -pe smp 16

scriptlist={1}
""".format(nJobs,scriptFile)

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
    
    return arrayScript
            
            
pdb_codes=["1MIX", "1P9G", "1UCS", "1XKR", "2BL2", "2EFR", "2FM9", "2JKU", "2QIH", "2QSK", "2UUI", "2XFD", "2YKT", "3CI9", "3CVF", "3GD8", "3GHF", "3HAP", "3HFE"]
#pdb_codes=["2BL2", "2FM9", "2JKU", "2QSK", "2UUI", "2XFD", "2YKT", "3CVF", "3GD8", "3GHF", "3HAP", "3HFE"]
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
        if r['Solution_Type']=='unfinished' or r['Solution_Type']=='no_job_directory' or (r['Solution_Type']=='PHASER_FAIL' and r['SHELXE_CC'] is not None):
           ensembles.append(ensemble)
           #if r['Solution_Type']=='unfinished' and os.path.isdir(search_dir): shutil.rmtree(search_dir)
           if os.path.isdir(search_dir): shutil.rmtree(search_dir)
           log=ensemble+".log"
           if os.path.isfile(log): os.unlink(log)

    # Remove old array scripts
    #old=glob.glob("arrayJob*")
    #for s in old:
    #    os.unlink(s)
    
    if len(ensembles):
        scripts=[os.path.abspath(e+".sh") for e in ensembles]
        jscript=prep_array(scripts,mrbd)
    
        # Submit job
        ample_util.run_command(["qsub",jscript])

