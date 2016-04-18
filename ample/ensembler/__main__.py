import argparse
import cPickle
import logging
import os
import sys

from ample.util import ample_util, config_util, exit_util
from ample.ensembler import ensembler_argparse
from ample.ensembler.ensembler_util import create_ensembles

ENSEMBLE_DIRNAME = 'ample_ensemble'

def setup_logging():
    """Set up file and console logging"""
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fl = logging.FileHandler(os.path.join(optd['work_dir'],"ensemble.log"))
    fl.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fl.setFormatter(formatter)
    logger.addHandler(fl)
    try:
        cl = logging.StreamHandler(stream=sys.stdout)
    except TypeError:
        cl = logging.StreamHandler(strm=sys.stdout)
    cl.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s\n') # Always add a blank line after every print
    cl.setFormatter(formatter)
    logger.addHandler(cl)
    
# Set up the command-line parsing
parser = argparse.ArgumentParser(description="AMPLE Ensembling Module")
ensembler_argparse.add_ensembler_options(parser, standalone=True)

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
        with open(optd['restart_pkl'], "r") as f: optd = cPickle.load(f)
    except Exception as e:
        msg = "Error unpickling ensemble pkl: {0}".format(e.message)
        exit_util.exit_error(msg, sys.exc_info()[2])
    restart = True
else:
    # We're running purely from command-line arguments
    amopt = config_util.AMPLEConfigOptions()
    amopt.populate(args)  
    optd = amopt.d 
        
# Set up the working directory if one doesn't already exist
if not ('work_dir' in optd and optd['work_dir'] and os.path.isdir(optd['work_dir'])):
    optd['work_dir'] = os.path.join(os.path.abspath(os.path.curdir),ENSEMBLE_DIRNAME)
    try:
        os.mkdir(optd['work_dir'])
    except OSError as e:
        msg = 'Error making ensemble workdir {0} : {1}'.format(optd['work_dir'],e)
        exit_util.exit_error(msg, sys.exc_info()[2])

setup_logging()
    
# Make sure we have models if in standalone mode
if not restart and not ('models' in optd and optd['models'] and os.path.exists(optd['models'])):
    msg = 'AMPLE ensembler requires a -models argument with a file/directory of pdbs'
    exit_util.exit_error(msg, sys.exc_info()[2])
    
try:
    if not restart:
        ample_util.extract_models(optd)
        # Hack - we'll be using gesamt to subcluster so it's not worth wiring in the stuff to
        # find maxcluster so we just assume it's there
        optd['maxcluster_exe'] = ample_util.find_exe('maxcluster')
        optd['ensemble_ok'] = os.path.join(optd['work_dir'],'ensemble.ok')
        optd['results_path'] = os.path.join(optd['work_dir'], "resultsd.pkl")
    create_ensembles(optd)
    ample_util.saveAmoptd(optd)
except Exception as e:
    msg = "Error running ensembling: {0}".format(e.message)
    exit_util.exit_error(msg, sys.exc_info()[2])
    