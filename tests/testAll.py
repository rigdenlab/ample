#!/usr/bin/env ccp4-python
"""Script to run all the AMPLE tests"""
import os
import sys
import unittest

###############################################################
#
# Add the python directory to the path
#
###############################################################
thisd =  os.path.abspath( os.path.dirname( __file__ ) )
paths = thisd.split( os.sep )
ampleDir = os.sep.join( paths[ : -1 ] )
sys.path.insert(0,os.path.join( ampleDir,'python'))

###############################################################
#
# The testsuite to hold the series of tests that we will run
#
###############################################################
testsuite = unittest.TestSuite()


###############################################################
#
# Add all the tests
#
###############################################################
import ample_ensemble
testsuite.addTests(ample_ensemble.testSuite())

import csymmatch
testsuite.addTests(csymmatch.testSuite())

import dssp
testsuite.addTests(dssp.testSuite())

import ensemble
testsuite.addTests(ensemble.testSuite())

import fasta_parser
testsuite.addTests(fasta_parser.testSuite())

import mtz_util
testsuite.addTests(mtz_util.testSuite())

import octopus_predict
testsuite.addTests(octopus_predict.testSuite())

import parse_arpwarp
testsuite.addTests(parse_arpwarp.testSuite())

import parse_buccaneer
testsuite.addTests(parse_buccaneer.testSuite())

import parse_phaser
testsuite.addTests(parse_phaser.testSuite())

import parse_refmac
testsuite.addTests(parse_refmac.testSuite())

import pdb_edit
testsuite.addTests(pdb_edit.testSuite())

import pdb_model
testsuite.addTests(pdb_model.testSuite())

import residue_map
testsuite.addTests(residue_map.testSuite())

import rio
testsuite.addTests(rio.testSuite())

import rosetta_model
testsuite.addTests(rosetta_model.testSuite())

import subcluster
testsuite.addTests(subcluster.testSuite())

import workers
testsuite.addTests(workers.testSuite())

###############################################################
#
# Now run 'em all
#
###############################################################
unittest.TextTestRunner(verbosity=2).run(testsuite) 
