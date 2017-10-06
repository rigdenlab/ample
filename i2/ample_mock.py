#!/usr/bin/env ccp4-python

import argparse
import copy
import cPickle
import logging
import os
import pyrvapi
import shutil
import sys
import time

#from ample.util import ample_util
from ample.util.ample_util import I2DIR, amoptd_fix_path
from ample.util.pyrvapi_results import AmpleOutput
from ample.constants import AMPLE_PKL

logging.basicConfig()

# Get the path to program.xml
parser = argparse.ArgumentParser()
parser.add_argument('-ample_pkl')
parser.add_argument('-ccp4i2_xml')
parser.add_argument('-rvapi_document', default=None)
parser.add_argument('-own_gui', default=False, action='store_true')
opt, _ = parser.parse_known_args()

# Create working directory
work_dir = os.path.abspath(I2DIR)
if os.path.isdir(work_dir): shutil.rmtree(work_dir)
os.mkdir(work_dir)

# Copy in amopt pkl
with open(opt.ample_pkl) as f: od = cPickle.load(f)

# update paths and copy across old files
amoptd_fix_path(od, newroot='/home/ccp4/tmp/ample',i2mock=False)
#amoptd_fix_path(od, newroot=work_dir, i2mock=True)

# Need to add these
od['no_gui'] = False
if opt.rvapi_document: od['rvapi_document'] = opt.rvapi_document
od['work_dir'] = work_dir
od['ccp4i2_xml'] = opt.ccp4i2_xml

with open(os.path.join(work_dir,AMPLE_PKL), 'w') as w: cPickle.dump(od, w)

# Run gui and create jsrview files from dict
AR = AmpleOutput(od, own_gui=opt.own_gui)
if True:    
    AR.display_results(od)
else:
    SLEEP = 2
    
    newd = copy.copy(od)
    newd['ensembles_data'] = None
    newd['mrbump_results'] = None
       
    AR.display_results(newd)
    time.sleep(SLEEP)
       
    #for i in range(10):
    newd['ensembles_data'] = od['ensembles_data']
    AR.display_results(newd)
    time.sleep(SLEEP)
    sys.exit()
       
    mrbump_results = []
    for r in ample_dict['mrbump_results'][0:3]:
        r['SHELXE_CC'] = None
        r['SHELXE_ACL'] = None
        mrbump_results.append(r)
    view1_dict['mrbump_results'] = mrbump_results
    AR.display_results(view1_dict)
    time.sleep(SLEEP)
       
    view1_dict['mrbump_results'] = ample_dict['mrbump_results'][0:5]
    AR.display_results(view1_dict)  
    time.sleep(SLEEP)
       
    view1_dict['mrbump_results'] = ample_dict['mrbump_results']
    AR.display_results(view1_dict) 

#pyrvapi.rvapi_store_document2('jens.rvapi')
