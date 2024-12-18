#!/usr/bin/env ccp4-python

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "01 Oct 2016"
__version__ = "1.0"

import argparse
import os
import sys

from ample.constants import AMPLE_PKL
from ample import ensembler
from ample.util import ample_util, config_util, exit_util, logging_util, process_models
from ample.util import argparse_util
from ample.util.options_processor import process_ensemble_options

ENSEMBLE_DIRNAME = 'ample_ensemble'

parser = argparse.ArgumentParser(description="AMPLE Ensembling Module")
argparse_util.add_general_options(parser)
argparse_util.add_cluster_submit_options(parser)
argparse_util.add_ensembler_options(parser)

# Get command-line arguments and see if we have a restart_pkl option as this
# is how we pass in an existing ample dictionary when we are running the ensembling
# as a standalone job on a cluster
optd = vars(parser.parse_args())

# Track restart as it determines if we need to unpack models
restart = False
if 'restart_pkl' in optd and optd['restart_pkl']:
    if not os.path.isfile(optd['restart_pkl']):
        msg = 'Cannot find ensemble pkl file: {0}'.format(optd['restart_pkl'])
        exit_util.exit_error(msg)
    try:
        optd = ample_util.read_amoptd(optd['restart_pkl'])
    except Exception as e:
        msg = "Error unpickling ensemble pkl: {0}".format(e.message)
        exit_util.exit_error(msg, sys.exc_info()[2])
    restart = True
else:
    # We're running purely from command-line arguments
    amopt = config_util.AMPLEConfigOptions()
    amopt.populate(optd)
    optd = amopt.d

# Start logging to the console
logger = logging_util.setup_console_logging()

# Make sure we have models if in standalone mode
if not restart and not ('models' in optd and optd['models'] and os.path.exists(optd['models'])):
    msg = 'AMPLE ensembler requires a -models argument with a file/directory of pdbs'
    exit_util.exit_error(msg, sys.exc_info()[2])

# Set up the working directory if one doesn't already exist
if not ('work_dir' in optd and optd['work_dir']):
    optd['work_dir'] = os.path.join(os.path.abspath(os.path.curdir), ENSEMBLE_DIRNAME)
if not os.path.isdir(optd['work_dir']):
    try:
        os.mkdir(optd['work_dir'])
    except OSError as e:
        msg = 'Error making ensemble workdir {0} : {1}'.format(optd['work_dir'], e)
        exit_util.exit_error(msg, sys.exc_info()[2])

assert os.path.isdir(optd['work_dir'])

# Start logging to a file
logging_util.setup_file_logging(os.path.join(optd['work_dir'], "ensemble.log"))
try:
    if not restart:
        results = process_models.extract_and_validate_models(optd)
        process_models.handle_model_import(optd, results)
        process_ensemble_options(optd)
        optd['ensemble_ok'] = os.path.join(optd['work_dir'], 'ensemble.ok')
        optd['results_path'] = os.path.join(optd['work_dir'], AMPLE_PKL)
    ensembler.create_ensembles(optd)
    ample_util.save_amoptd(optd)
except Exception as e:
    msg = "Error running ensembling: {0}".format(e.message)
    exit_util.exit_error(msg, sys.exc_info()[2])
