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
parser.add_argument('-ccp4i2_xml')
parser.add_argument('-rvapi_document')
opt, other = parser.parse_known_args()
logging.debug("Script {0} got known arguments: {1} unknown {2}".format(sys.argv[0], opt, other))

ample_pkl = '/opt/ample.git/from_existing_models/resultsd.pkl'
show_gui = True
work_dir = os.getcwd()
new_root = os.sep.join(os.path.dirname(ample_pkl).split(os.sep)[:-1])
i2mock = False

if not (ample_pkl and os.path.isfile(ample_pkl)):
    sys.stderr.write("Cannot find AMPLE pkl file: {0}\n".format(ample_pkl))
    sys.exit(1)

# Load AMPLE dictionary
with open(ample_pkl) as f: old_dict = cPickle.load(f)

if opt.rvapi_document:
    old_dict['rvapi_document'] = opt.rvapi_document
elif opt.ccp4i2_xml:
    # Create working directory
    work_dir = os.path.abspath(I2DIR)
    if os.path.isdir(work_dir): shutil.rmtree(work_dir)
    logging.info("Making work directory: {0}".format(work_dir))
    os.mkdir(work_dir)
    new_root = work_dir
    i2mock = True
    logging.info("Setting ccp4i2_xml: {0}".format(opt.ccp4i2_xml))
    old_dict['ccp4i2_xml'] = opt.ccp4i2_xml

amoptd_fix_path(old_dict, newroot=new_root, i2mock=i2mock)
# Need to add these
old_dict['show_gui'] = show_gui
old_dict['work_dir'] = work_dir
with open(os.path.join(work_dir, AMPLE_PKL), 'w') as w: cPickle.dump(old_dict, w)

# Run gui and create jsrview files from dict
AR = AmpleOutput(old_dict)
if True:
    AR.display_results(old_dict)
else:
    SLEEP = 2

    newd = copy.copy(old_dict)
    newd['ensembles_data'] = None
    newd['mrbump_results'] = None

    AR.display_results(newd)
    time.sleep(SLEEP)

    #for i in range(10):
    newd['ensembles_data'] = old_dict['ensembles_data']
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
