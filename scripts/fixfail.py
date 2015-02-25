import cPickle
import glob
import os
import shutil
import sys

ample_root="/home/jmht42/ample-dev1"
sys.path.insert(0,os.path.join(ample_root,"python"))
mrbumpp="/home/jmht42/mrbump-trunk/include/parsers"
sys.path.insert(0,mrbumpp)

import ample_util
import mrbump_results
import parse_phaser
import parse_refmac
            
pdb_codes=["1MIX", "1P9G", "1UCS", "1XKR", "2BL2", "2EFR", "2FM9", "2JKU", "2QIH", "2QSK", "2UUI", "2XFD", "2YKT", "3CI9", "3CVF", "3GD8", "3GHF", "3HAP", "3HFE"]
#pdb_codes=["1UCS"]
root="/volatile/jmht42/testset/percent"

for pdb in pdb_codes:
    print "PDB ",pdb
    mrbd=os.path.join(root,pdb,"ROSETTA_MR_0","MRBUMP")
    os.chdir(mrbd)
    
    res_sum = mrbump_results.ResultsSummary()
    res_sum.extractResults(mrbd)
    for r in res_sum.results:
        if r['Solution_Type']=='PHASER_FAIL' or r['Solution_Type']=='REFMAC_FAIL':
            search_dir=r['Search_directory']
            pkl=os.path.join(mrbd,search_dir,'results','resultsTable.pkl')
            with open(pkl) as f: dr=cPickle.load(f)
            name=dr.keys()[0]
            pd=dr[name]['PHASER']
            #for k in sorted(pd.keys()): print k,pd[k]
            #jd=pd['JobDirectory']
            pl=pd['PHASER_logfile']
            pp=parse_phaser.PhaserLogParser(pl)
            pd["PHASER_version"] = pp.version
            pd["PHASER_LLG"] = pp.LLG
            pd["PHASER_TFZ"] = pp.TFZ
            pd["PHASER_RFZ"] = pp.RFZ
            pd["PHASER_time"] = pp.time
            pd["PHASER_killed"] = pp.killed
            pd["PHASER_error"] = pp.errorStr
            #print "GOT TFZ for {0}: {1}".format(name,pd["PHASER_TFZ"])
            print "GOT OLD SOLUTION ",pd['Solution_Type']
            if pd['PHASER_error']:
                print "GOT ERROR for {0}: {1}".format(name,pd["PHASER_error"])
                pd['Solution_Type']='PHASER_ERROR'
            elif not pd["PHASER_LLG"] and not pd['SHELXE_CC']:
                pd['Solution_Type']='PHASER_FAIL'
            else:
                rl=pd['REFMAC_logfile']
                rp=parse_refmac.RefmacLogParser(rl)
                pd['final_Rfact']=rp.finalRfact
                pd['final_Rfree']=rp.finalRfree
                diff=rp.finalRfree-rp.initRfree
                if pd['final_Rfree'] <= 0.35:
                    pd['Solution_Type']='GOOD'
                elif (pd['final_Rfree'] <= 0.5 and diff  <= 0.8) or pd['final_Rfree'] <= 0.48 or (pd['final_Rfree'] <= 0.55 and diff  <= 0.95):
                    pd['Solution_Type']='MARGINAL'
                else:
                    pd['Solution_Type']='POOR'
            print "GOT NEW SOLUTION ",pd['Solution_Type'],pd['final_Rfree']
            #with open(pkl,'w') as f: cPickle.dump(dr,f)
            
            

