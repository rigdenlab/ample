"""Module containing a framework for integration testing of AMPLE modules"""

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "25 Mar 2016"
__version__ = "1.0"

from unittest import TestCase, TestLoader, TextTestRunner, TestSuite
import glob
import imp
import logging
import os
import shutil
import sys

from ample.constants import SHARE_DIR, AMPLE_PKL
from ample.testing.constants import CLUSTER_ARGS, EXTRA_ARGS
from ample.util import ample_util
from pyjob.factory import TaskFactory
from pyjob.script import ScriptCollector, Script


ENSEMBLER = 'ensembler'
MODELLING = 'modelling'

# Any modules required when ample results dictionaries are unpickled should be added here
import ample.ensembler

logger = logging.getLogger(__name__)

# Available packages. Hard-coded for now to show visually what we have in
# argparse module. Not needed otherwise
PACKAGES = [
    'from_existing_models',
    'from_quark_models',
    'from_single_model',
    'homologs',
    'ideal_helices',
    'import_cluster',
    'import_ensembles',
    'import_models',
    'nmr_truncate',
]
if not sys.platform.startswith("win"):
    PACKAGES += [
        'nmr_remodel',
        'rosetta_contacts',
        'rosetta_contacts_subselect',
        'rosetta_modelling',
        'rosetta_restraints',
    ]


def add_cmd_options(parser):
    parser.add_argument('-clean', action='store_true', default=False, help="Clean up all test files/directories")
    parser.add_argument('-nproc', type=int, default=1, help="Number of processors to run on (1 per job)")
    parser.add_argument('-dry_run', action='store_true', default=False, help="Don\'t actually run the jobs")
    parser.add_argument('-rosetta_dir', help="Location of rosetta installation directory")
    parser.add_argument('test_cases', nargs='*', help="[ {0} ]".format(" | ".join(PACKAGES)))
    parser.add_argument('-run_dir', type=str, default=None, help="directory to run jobs in")


class AMPLEBaseTest(TestCase):
    RESULTS_PKL = None
    AMPLE_DICT = None

    def setUp(self):
        self.assertTrue(os.path.isfile(self.RESULTS_PKL), "Missing pkl file: {0}".format(self.RESULTS_PKL))
        try:
            self.AMPLE_DICT = ample_util.read_amoptd(self.RESULTS_PKL)
        except ImportError as e:
            logger.exception(
                "Error importing module while unpickling ample results dictionary: '{}'"
                "Add any imports required to the module: {}".format(e, os.path.abspath(__file__))
            )
            raise (e)


class AMPLEIntegrationFramework(object):
    """Framework to run Ample integration testing"""

    def __init__(self, test_cases=None, run_dir=None):
        examples_dir = os.path.join(SHARE_DIR, "examples")
        self.test_dict = SuiteLoader().load_cases(examples_dir, test_cases=test_cases)
        if not len(self.test_dict):
            if len(test_cases):
                msg = 'Could not find test cases {0} in directory {1}'.format(test_cases, examples_dir)
            else:
                msg = "Could not find any test cases in directory: {0}".format(examples_dir)
            raise RuntimeError(msg)

        # Make a directory to keep all files together
        _root = os.path.abspath(run_dir) if run_dir else self.get_run_dir()
        self.run_dir = os.path.join(_root, "ample_testing")
        if not os.path.isdir(self.run_dir):
            os.mkdir(self.run_dir)

    def get_run_dir(self):
        return os.getcwd()

    def clean(self, clean_all=True, clean_dir=False):
        for name in self.test_dict.keys():
            os.chdir(self.run_dir)
            logger.info("Cleaning {0} in directory {1}".format(name, self.run_dir))
            work_dir = os.path.join(self.run_dir, name)
            if os.path.isdir(work_dir):
                shutil.rmtree(work_dir)
            logfile = work_dir + '.log'
            if os.path.isfile(logfile):
                os.unlink(logfile)
            if clean_all:
                script = work_dir + ample_util.SCRIPT_EXT
                if os.path.isfile(script):
                    os.unlink(script)
        if clean_dir and os.path.isdir(self.run_dir):
            shutil.rmtree(self.run_dir)

    def run(self, nproc=1, dry_run=False, clean_up=True, rosetta_dir=None, **kwargs):
        """Run the integration testing jobs and then the unittests to test them.
        
        In all cases jobs are run on a single processor. For running on a cluster, the
        ample job scripts have the queue directives added to them, and  each ample job
        is launched on the head node. The individual jobs then submit their various job
        stages to the queue and the integration test job just manages running all the 
        individual ample jobs until they have finished. Although this means lots of jobs
        running on the head node, the actual computation done on the head node should be minimal
        as all processing is submitted to the queue.
        
        Previously when running on a cluster we created a single single-processor serial ample script 
        for each job and then submitted a single array job to run all the jobs on the cluster. This
        approach had to be abandoned as (I think) the individual jobs timed out. 
        """
        logger.info("Writing files to: {0}".format(self.run_dir))

        if dry_run:
            clean_up = False

        if rosetta_dir and not os.path.isdir(rosetta_dir):
            logger.debug("Cannot find rosetta_dir: {0}".format(rosetta_dir))
            sys.exit(1)

        if clean_up:
            self.clean()

        collector = self._create_scripts(rosetta_dir, **kwargs)
        if not len(collector.scripts):
            raise RuntimeError("Could not find any test cases to run!")

        logger.info("The following test cases will be run:")
        for name in self.test_dict.keys():
            logger.info("{0}: {1}".format(name, self.run_dir))

        nprocesses = nproc
        ## Run all the jobs
        # If we're running on a cluster, we run on as many processors as there are jobs,
        # as the jobs are just sitting and monitoring the queue
        if kwargs['submit_qtype'] != 'local':
            logger.info("Jobs will be submitted to a cluster queueing system")
            nprocesses = len(collector.scripts)

        if not dry_run:
            with TaskFactory(
                    kwargs['submit_qtype'],
                    collector,
                    environment=kwargs['submit_pe'],
                    name='test',
                    processes=nprocesses,
                    max_array_size=kwargs['submit_max_array'],
                    queue=kwargs['submit_queue'],
                    shell="/bin/bash",
            ) as task:
                task.run()
                task.wait(interval=5, monitor_f=None)

        # Now check the results using the unittesting framework
        self.run_unittest_suite()
        return

    def _create_scripts(self, rosetta_dir, **kwargs):
        """Create scripts and set path to working directory"""
        collector = ScriptCollector(None)
        owd = os.getcwd()
        for name in self.test_dict.keys():
            os.chdir(self.run_dir)
            work_dir = os.path.join(self.run_dir, name)
            args = self.test_dict[name]['args']

            # Rosetta is the only think likely to change between platforms so we update the entry
            if rosetta_dir and self._is_in_args('-rosetta_dir', args):
                args = self._update_args(args, [['-rosetta_dir', rosetta_dir]])
            # Additional argumenst for submitting to a cluster
            args = self._update_cluster_args(args, **kwargs)
            if EXTRA_ARGS:
                args = self._update_args(args, EXTRA_ARGS)

            # We track different modules using the name of the test case
            if name.startswith(ENSEMBLER):
                testcase_type = ENSEMBLER
            elif name.startswith(MODELLING):
                testcase_type = MODELLING
            else:
                testcase_type = 'ample'
            if testcase_type != 'ample' and sys.platform.startswith('win'):
                logger.critical("Cannot run module testcases on windows due to multiprocessing bug")
                continue
            script = self.write_script(self.run_dir, name, args + [['-work_dir', work_dir]], testcase_type)
            collector.add(script)
            # Set path to the directory the case is run so we can pass it to the unittest
            self.test_dict[name]['work_dir'] = work_dir

            # Run the setup function if one is provided
            if 'setup' in self.test_dict[name] and callable(self.test_dict[name]['setup']):
                self.test_dict[name]['setup'](self.run_dir)

            os.chdir(owd)  # Back to where we started
        return collector

    def _is_in_args(self, argt, args):
        if type(argt) is str:
            key = argt
        else:
            key = argt[0]
        return key in [a[0] for a in args]

    def _replace_arg(self, new_arg, args):
        for i, a in enumerate(args):
            if a[0] == new_arg[0]:
                args[i] = new_arg
                return args
        assert False

    def _update_args(self, args, new_args):
        """Add/update any args"""
        for argt in new_args:
            if not self._is_in_args(argt, args):
                args.append(argt)
            else:
                self._replace_arg(argt, args)
        return args

    def _update_cluster_args(self, args, **kwargs):
        """Add the cluster submission arguments
        
        See if any of the clustering submission arguments are in **kwarg and append
        any non-None ones to args. Otherwise we use the non-None arguments from CLUSTER_ARGS"""
        if not 'submit_cluster' in kwargs or not kwargs['submit_cluster']:
            return args
        for k, v in kwargs.iteritems():
            value = None
            if k in CLUSTER_ARGS.keys():
                if v is not None:
                    value = v
                elif CLUSTER_ARGS[k] is not None:
                    value = CLUSTER_ARGS[k]
            if value:
                # Need to add the hypen on to the key so it can be used as a command-line arg
                args.append(["-" + k, value])
        return args

    def run_unittest_suite(self):
        suite = TestSuite()
        for name in self.test_dict.keys():
            testClass = self.test_dict[name]['test']
            testClass.WORK_DIR = self.test_dict[name]['work_dir']
            testClass.RESULTS_PKL = os.path.join(testClass.WORK_DIR, AMPLE_PKL)
            _suite = TestLoader().loadTestsFromTestCase(testClass)
            suite.addTests(_suite)
        TextTestRunner(verbosity=2).run(suite)

    def write_script(self, work_dir, name, args, testcase_type):
        """Write script"""
        linechar = "^" if sys.platform.startswith('win') else "\\"

        script = Script(directory=work_dir, stem=name)
        test_exe = os.path.join(os.environ["CCP4"], "bin", "ample")
        test_exe = test_exe + ample_util.SCRIPT_EXT if sys.platform.startswith("win") else test_exe
        if testcase_type == ENSEMBLER:
            test_exe = '{0} -m ample.ensembler'.format(os.path.join(os.environ["CCP4"], "bin", "ccp4-python"))
        elif testcase_type == MODELLING:
            test_exe = '{0} -m ample.modelling'.format(os.path.join(os.environ["CCP4"], "bin", "ccp4-python"))

        # All arguments need to be strings
        args = [map(str, a) for a in args]
        script.append("{0} {1}".format(test_exe, linechar))
        for argt in args:
            script.append(" ".join(argt) + " " + linechar)

        return script


class SuiteLoader(object):
    """Loader designed to obtain all test cases in a package"""

    def load_cases(self, directory, test_cases=None, pattern="test_cases"):
        """function to load a integration test suite"""
        search_pattern = os.path.join(directory, "*")
        cases = [os.path.basename(folder) for folder in glob.glob(search_pattern) if os.path.isdir(folder)]
        test_dict = self._load_cases(cases, directory, pattern)
        # Needs to follow as case names will be folders and test_cases the
        # actual cases themselves
        if test_cases:
            test_dict = {k: v for k, v in test_dict.iteritems() if k in test_cases}
        return test_dict

    def _load_cases(self, cases, directory, pattern):
        test_cases = {}
        for example_dir in cases:
            path = os.path.join(directory, example_dir)
            test_module = self.load_module(pattern, [path])
            # Skip anything that's not a valid AMPLE test module
            if not test_module or not hasattr(test_module, 'TEST_DICT'):
                continue
            for k, v in test_module.TEST_DICT.iteritems():
                if k in test_cases:
                    raise RuntimeError("Duplicate key: {0}".format(k))
                test_cases[k] = v
        return test_cases

    def load_module(self, mod_name, paths):
        try:
            mfile, pathname, desc = imp.find_module(mod_name, paths)
        except ImportError:
            logger.critical("Cannot find test module in {0}".format(mod_name))
            return None
        try:
            test_module = imp.load_module(mod_name, mfile, pathname, desc)
        except Exception as e:
            logger.critical("Error loading test case from directory: {0}\n {1}\n".format(paths, e))
            raise Exception(e)
        finally:
            mfile.close()
        return test_module
