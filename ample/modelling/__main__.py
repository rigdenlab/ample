#!/usr/bin/env ccp4-python

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "22 Mar 2019"
__version__ = "1.0"

import argparse
import os
import logging
import multiprocessing

from ample.util import argparse_util
from ample.util.sequence_util import Sequence
from ample.util import logging_util
from ample.modelling.rosetta_model import RosettaModel


def process_args(args):
    if args.rosetta_flagsfile:
        args.rosetta_flagsfile = os.path.abspath(args.rosetta_flagsfile)
    if args.nproc is None:
        if args.submit_cluster:
            args.nproc = 1
        else:
            args.nproc = multiprocessing.cpu_count()


# Handle the command-line
parser = argparse.ArgumentParser(description="AMPLE Modelling Module")
# Need to add seperately here as is usually part of MR options
parser.add_argument('-fasta', help='protein fasta file.')
parser.add_argument(
    '-nchains', help='The number of chains to select from the multimer for the final single-chain models'
)
argparse_util.add_core_options(parser)
argparse_util.add_rosetta_options(parser)
argparse_util.add_cluster_submit_options(parser)

work_dir = os.path.abspath('rosetta_modelling')
parser.set_defaults(submit_cluster=False, submit_qtype='SGE', submit_array=True, nmodels=1000, work_dir=work_dir)
args = parser.parse_args()
process_args(args)

# Start logging to the console
logger = logging_util.setup_console_logging()
logger.info("*** AMPLE ROSETTA modelling package ***")

if not os.path.isdir(args.work_dir):
    os.mkdir(args.work_dir)

rm = RosettaModel()
if args.rosetta_dir and os.path.isdir(args.rosetta_dir):
    rm.set_paths(rosetta_dir=args.rosetta_dir)

if args.fasta:
    rm.fasta = args.fasta
    fp = Sequence(fasta=args.fasta, canonicalise=False)
    rm.sequence_length = fp.length()

rm.nmodels = args.nmodels
rm.work_dir = args.work_dir
rm.models_dir = os.path.join(args.work_dir, "models")
rm.frags_3mers = args.frags_3mers
rm.frags_9mers = args.frags_9mers
rm.multimer_modelling = args.multimer_modelling
rm.nchains = args.nchains

rm.nproc = args.nproc
rm.submit_cluster = args.submit_cluster
rm.submit_qtype = args.submit_qtype
rm.submit_queue = args.submit_queue
rm.submit_array = args.submit_array
rm.submit_max_array = args.submit_max_array

logger.info("Running binary {} with flagsfile: {}".format(args.rosetta_executable, args.rosetta_flagsfile))
if args.multimer_modelling:
    rm.do_multimer_modelling()
else:
    rm.model_from_flagsfile(args.rosetta_flagsfile, args.rosetta_executable)
