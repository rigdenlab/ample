#!/usr/bin/env ccp4-python

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "01 Oct 2016"
__version__ = "1.0"

import argparse
import os
import sys

from ample import ensembler
from ample.util import ample_util, config_util, exit_util, logging_util
from ample.util import argparse_util

ENSEMBLE_DIRNAME = 'ample_ensemble'
    
# Set up the command-line parsing
parser = argparse.ArgumentParser(description="AMPLE Ensembling Module")
# Add options for running as a standalone module
argparse_util.add_general_options(parser)
argparse_util.add_cluster_submit_options(parser)
# Ensemble options
ensembler.add_argparse_options(parser)

# Get command-line arguments and see if we have a restart_pkl option as this
# is how we pass in an existing ample dictionary when we are running the ensembling
# as a standalone job on a cluster
args = parser.parse_args()
optd = vars(args)

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
    amopt.populate(args)  
    optd = amopt.d 

# Start logging to the console
logging_util.setup_console_logging()

# Make sure we have models if in standalone mode
if not restart and not ('models' in optd and optd['models'] and os.path.exists(optd['models'])):
    msg = 'AMPLE ensembler requires a -models argument with a file/directory of pdbs'
    exit_util.exit_error(msg, sys.exc_info()[2])

# Set up the working directory if one doesn't already exist
if not ('work_dir' in optd and optd['work_dir']):
    optd['work_dir'] = os.path.join(os.path.abspath(os.path.curdir),ENSEMBLE_DIRNAME)
if not os.path.isdir(optd['work_dir']):
    try:
        os.mkdir(optd['work_dir'])
    except OSError as e:
        msg = 'Error making ensemble workdir {0} : {1}'.format(optd['work_dir'],e)
        exit_util.exit_error(msg, sys.exc_info()[2])

assert os.path.isdir(optd['work_dir'])

# Start logging to a file
logging_util.setup_file_logging(os.path.join(optd['work_dir'],"ensemble.log"))
try:
    if not restart:
        ample_util.extract_models(optd)
        if optd['subcluster_program'] == 'gesamt':
            optd['gesamt_exe'] = ample_util.find_exe('gesamt')
        elif optd['subcluster_program'] == 'gesamt':
            optd['maxcluster_exe'] = ample_util.find_exe('maxcluster')
        else:
            raise RuntimeError("Unknown subcluster_program: {0}".format(optd['subcluster_program']))
        optd['theseus_exe'] = ample_util.find_exe('theseus')
        optd['ensemble_ok'] = os.path.join(optd['work_dir'],'ensemble.ok')
        optd['results_path'] = os.path.join(optd['work_dir'], "resultsd.pkl")
    ensembler.create_ensembles(optd)
    ample_util.save_amoptd(optd)
except Exception as e:
    msg = "Error running ensembling: {0}".format(e.message)
    exit_util.exit_error(msg, sys.exc_info()[2])

