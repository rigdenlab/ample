import cPickle
import glob
import os
import shutil
import sys

ample_root="/home/jmht42/ample-dev1"
sys.path.insert(0,os.path.join(ample_root,"python"))

import ample_util
import mrbump_results

            
#pdb_codes=["1MIX", "1P9G", "1UCS", "1XKR", "2BL2", "2EFR", "2FM9", "2JKU", "2QIH", "2QSK", "2UUI", "2XFD", "2YKT", "3CI9", "3CVF", "3GD8", "3GHF", "3HAP", "3HFE"]
pdb_codes=["1MIX"]
root="/volatile/jmht42/testset/percent"

for pdb in pdb_codes:
    mrbd=os.path.join(root,pdb,"ROSETTA_MR_0","MRBUMP")
    os.chdir(mrbd)
    
    res_sum = mrbump_results.ResultsSummary()
    res_sum.extractResults(mrbd)
    for r in res_sum.results:
        if r['Solution_Type']=='PHASER_FAIL' or r['Solution_Type']=='REFMAC_FAIL':
            search_dir=r['Search_directory']
            pkl=os.path.join(mrbd,search_dir,'results','resultsTable.pkl')
            with open(pkl) as f: rd=cPickle.load(f)
            name=dr.keys()[0]
            pd=dr[name]['PHASER']
            jd=pd['Job_directory']
            if pd['PHASER_error']:
                print "FOUND FAIL FOR {0}: {1}".format(name,pd['PHASER_error'])
            
            
#         if r['Solution_Type']=='unfinished' or r['Solution_Type']=='no_job_directory' or (r['Solution_Type']=='PHASER_FAIL' and r['SHELXE_CC'] is not None):


