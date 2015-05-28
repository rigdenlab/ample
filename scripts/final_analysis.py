import cPickle
import os
import shutil
import sys

ample_root="//opt/ample-dev1"
sys.path.insert(0,os.path.join(ample_root,"python"))
mrbumpp="/opt/mrbump-trunk/include/parsers"
sys.path.insert(0,mrbumpp)

import ample_util
import benchmark
import mrbump_results

root="/media/data/shared/testset/truncate_single"
#pdb_codes=["1MIX", "1P9G", "1UCS", "1XKR", "2BL2", "2EFR", "2FM9", "2JKU", "2QIH", "2QSK", "2UUI", "2XFD", "2YKT", "3CI9", "3CVF", "3GD8", "3GHF", "3HAP", "3HFE"]
pdb_codes=["2UUI", "2YKT", "3CI9", "3CVF", "3GD8", "3GHF", "3HAP"]
for pdb in pdb_codes:
    print "PDB ",pdb
    ample_dir=os.path.join(root,pdb,"ROSETTA_MR_0")
    if not os.path.isdir(ample_dir): ample_dir=os.path.join(root,pdb,"AMPLE_0")
    pkl=os.path.join(ample_dir,'resultsd.pkl')
    mrbumpd=os.path.join(ample_dir,'MRBUMP')
    with open(pkl) as f: amoptd=cPickle.load(f)
    res_sum = mrbump_results.ResultsSummary()
    res_sum.extractResults(mrbumpd)
    amoptd['mrbump_results'] = res_sum.results
    bd=os.path.join(ample_dir,"benchmark")
    shutil.rmtree(bd)
    os.mkdir(bd)
    benchmark.analyse(amoptd,newroot=ample_dir)
    amoptd['results_path']=os.path.join(ample_dir,"resultsd.pkl")
    ample_util.saveAmoptd(amoptd)

