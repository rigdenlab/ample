#!/usr/bin/env ccp4-python

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "22 Mar 2019"
__version__ = "1.0"

import argparse
import os
import logging

from ample.util import argparse_util
from ample.util import logging_util
from ample.modelling.rosetta_model import RosettaModel

parser = argparse.ArgumentParser(description="AMPLE Modelling Module")
argparse_util.add_core_options(parser)
argparse_util.add_rosetta_options(parser)
argparse_util.add_cluster_submit_options(parser)

# Start logging to the console
logging_util.setup_console_logging()

logger = logging.getLogger()
logger.info("*** AMPLE ROSETTA modelling package ***")

# Get cmdline options or set defaults
args = parser.parse_args()

if args.rosetta_flagsfile is None:
    raise RuntimeError("Need to supply a ROSETTA flagsfile with the -rosetta_flagsfile argument")
flagsfile = os.path.abspath(args.rosetta_flagsfile)
rosetta_binary = args.rosetta_AbinitioRelax


nproc = args.nproc
if nproc is None:
    nproc = 1
submit_cluster = args.submit_cluster
if submit_cluster is None:
    submit_cluster = False
submit_qtype = args.submit_qtype
if submit_qtype is None:
    submit_qtype = 'SGE'
submit_queue = args.submit_queue
submit_array = args.submit_array
if submit_array is None:
    submit_array = True
submit_max_array = args.submit_max_array
rosetta_dir = args.rosetta_dir

nmodels = args.nmodels
if nmodels is None:
    nmodels = 1000
work_dir = args.work_dir
if work_dir is None:
    work_dir = os.path.abspath('rosetta_modelling')
if not os.path.isdir(work_dir):
    os.mkdir(work_dir)
models_dir = os.path.join(work_dir, "models")

rm = RosettaModel()
if rosetta_dir and os.path.isdir(rosetta_dir):
    rm.set_paths(rosetta_dir=rosetta_dir)

rm.nmodels = nmodels
rm.work_dir = work_dir
rm.models_dir = models_dir

rm.nproc = nproc
rm.submit_cluster = submit_cluster
rm.submit_qtype = submit_qtype
rm.submit_queue = submit_queue
rm.submit_array = submit_array
rm.submit_max_array = submit_max_array

logger.info("Running binary {} with flagsfile: {}".format(rosetta_binary, flagsfile))
rm.model_from_flagsfile(flagsfile, rosetta_binary)
