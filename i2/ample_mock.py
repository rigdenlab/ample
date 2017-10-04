#!/usr/bin/env ccp4-python

import argparse
import cPickle
import logging
import os
import shutil
import sys

#from ample.util import ample_util
from ample.util.ample_util import I2DIR, amoptd_fix_path
from ample.util.pyrvapi_results import AmpleOutput

logging.basicConfig()

# Get the path to program.xml
parser = argparse.ArgumentParser()
parser.add_argument('-ccp4i2_xml')
opt, _ = parser.parse_known_args()

opkl = '/opt/ample.git/ample_testing/from_existing_models/resultsd.pkl'
work_dir = os.path.abspath(I2DIR)
# Create working directory
os.mkdir(work_dir)

# Copy in amopt pkl
with open(opkl) as f: od = cPickle.load(f)

# update paths and copy across old files
amoptd_fix_path(od, newroot=work_dir, i2mock=True)

# Need to create
od['work_dir'] = work_dir
od['ccp4i2_xml'] = opt.ccp4i2_xml

# Run gui and create jsrview files from dict
AR = AmpleOutput(od)
AR.display_results(od)
