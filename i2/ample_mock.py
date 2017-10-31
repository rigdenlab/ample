#!/usr/bin/env ccp4-python

import argparse
import copy
import cPickle
import logging
import os
import shutil
import sys
import time

#from ample.util import ample_util
from ample.util.ample_util import I2DIR, amoptd_fix_path
from ample.util.pyrvapi_results import AmpleOutput
from ample.constants import AMPLE_PKL

logging.basicConfig(level=logging.DEBUG)

# Get the path to program.xml
parser = argparse.ArgumentParser()
parser.add_argument('-pkl')
parser.add_argument('-ccp4i2_xml')
parser.add_argument('-rvapi_document', default=None)
parser.add_argument('-own_gui', default=False, action='store_true')
opt, other = parser.parse_known_args()
logging.debug("Script {0} got known arguments: {1} unknown {2}".format(sys.argv[0], opt, other))

# Copy in amopt pkl
ample_pkl = opt.pkl
if opt.rvapi_document:
    mroot='/home/jscofe/tmp/ample'
    #mroot = '/opt/ample.git/ample_testing'
    ample_pkl =  os.path.join(mroot,'from_existing_models','resultsd.pkl')

if not os.path.isfile(ample_pkl):
    sys.stderr.write("Cannot find AMPLE pkl file: {0}\n".format(ample_pkl))
    sys.exit(1)

# Load AMPLE dictionary
with open(ample_pkl) as f: od = cPickle.load(f)

if opt.rvapi_document:
    amoptd_fix_path(od, newroot=mroot, i2mock=False)
else:
    # Create working directory
    work_dir = os.path.abspath(I2DIR)
    if os.path.isdir(work_dir): shutil.rmtree(work_dir)
    logging.info("Making work directory: {0}".format(work_dir))
    os.mkdir(work_dir)
    
    # update paths and copy files into run directory
    amoptd_fix_path(od, newroot=work_dir, i2mock=True)

# Need to add these
od['no_gui'] = False
if opt.rvapi_document:
    # JSCOFE HACK
    od['rvapi_document'] = opt.rvapi_document
    work_dir = os.getcwd()
od['work_dir'] = work_dir
if opt.ccp4i2_xml: logging.info("Setting ccp4i2_xml: {0}".format(opt.ccp4i2_xml))
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
