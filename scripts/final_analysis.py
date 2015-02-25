import cPickle
import os
import shutil
import sys

ample_root="/home/jmht42/ample-dev1"
sys.path.insert(0,os.path.join(ample_root,"python"))
mrbumpp="/home/jmht42/mrbump-trunk/include/parsers"
sys.path.insert(0,mrbumpp)

import ample_util
import benchmark
import mrbump_results

root="/volatile/jmht42/testset/percent"
#pdb_codes=["1MIX", "1P9G", "1UCS", "1XKR", "2BL2", "2EFR", "2FM9", "2JKU", "2QIH", "2QSK", "2UUI", "2XFD", "2YKT", "3CI9", "3CVF", "3GD8", "3GHF", "3HAP", "3HFE"]
pdb_codes=["1MIX", "1P9G", "1UCS",  "2BL2", "2EFR", "2FM9", "2JKU", "2QIH", "2QSK", "2UUI", "2XFD", "2YKT", "3CI9", "3CVF", "3GD8", "3GHF", "3HAP", "3HFE"]
for pdb in pdb_codes:
    print "PDB ",pdb
    ample_dir=os.path.join(root,pdb,"ROSETTA_MR_0")
    pkl=os.path.join(ample_dir,'resultsd.pkl')
    mrbumpd=os.path.join(ample_dir,'MRBUMP')
    with open(pkl) as f: amoptd=cPickle.load(f)
    res_sum = mrbump_results.ResultsSummary()
    res_sum.extractResults(mrbumpd)
    amoptd['mrbump_results'] = res_sum.results
    shutil.rmtree(amoptd['benchmark_dir'])
    os.mkdir(amoptd['benchmark_dir'])
    benchmark.analyse(amoptd)
    ample_util.saveAmoptd(amoptd)

