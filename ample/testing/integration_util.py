

from unittest import TestCase, TestLoader, TextTestRunner, TestSuite
import glob
import imp
import logging
import os
import shutil
import sys

from ample.constants import SHARE_DIR
from ample.testing.constants import CLUSTER_ARGS, EXTRA_ARGS
from ample.util import ample_util
from ample.util import workers_util

__author__ = "Felix Simkovic and Jens Thomas"
__date__ = "25-Mar-2016"

logger = logging.getLogger(__name__)

# Available packages. Hard-coded for now to show visually what we have in
# argparse module. Not needed otherwise
PACKAGES = ['from_existing_models', 'from_quark_models', 'from_single_model', 
            'homologs', 'ideal_helices',  'import_cluster', 'import_ensembles', 
            'import_models', 'missing_domain', 'nmr_truncate']
if not sys.platform.startswith("win"):
    PACKAGES += ['nmr_remodel', 'rosetta_contacts', 'rosetta_modelling', 
                 'rosetta_restraints']

def add_cmd_options(parser):
    parser.add_argument('-clean', action='store_true', default=False,
                        help="Clean up all test files/directories")
    parser.add_argument('-nproc', type=int, default=1,
                        help="Number of processors to run on (1 per job)")
    parser.add_argument('-dry_run', action='store_true', default=False,
                        help="Don\'t actually run the jobs")
    parser.add_argument('-rosetta_dir',
                        help="Location of rosetta installation directory")
    parser.add_argument('test_cases', nargs='*',
                        help="[ {0} ]".format(" | ".join(PACKAGES)))
    parser.add_argument('-run_dir', type=str, default=None,
                        help="directory to run jobs in")


class AMPLEBaseTest(TestCase):
    RESULTS_PKL = None
    AMPLE_DICT = None
    def setUp(self):
        self.assertTrue(os.path.isfile(self.RESULTS_PKL), "Missing pkl file: {0}".format(self.RESULTS_PKL))
        self.AMPLE_DICT = ample_util.read_amoptd(self.RESULTS_PKL)


class AMPLEIntegrationFramework(object):
    """Framework to run Ample integration testing"""
    
    def __init__(self, test_cases=None, run_dir=None):
        examples_dir = os.path.join(SHARE_DIR, "examples")
        self.test_dict = SuiteLoader().load_cases(examples_dir, 
                                                  test_cases=test_cases)
        if not len(self.test_dict):
            if len(test_cases):
                msg = 'Could not find test cases {0} in directory {1}'.format(test_cases,examples_dir)
            else:
                msg = "Could not find any test cases in directory: {0}".format(examples_dir)
            raise RuntimeError(msg)
        
        # Make a directory to keep all files together
        _root = os.path.abspath(run_dir) if run_dir else self.get_run_dir()
        self.run_dir = os.path.join(_root, "ample_testing")
        if not os.path.isdir(self.run_dir): os.mkdir(self.run_dir)
    
    def get_run_dir(self):
#         # OS X has problems with relative paths - symlink in root messes
#         # up the search in theseus - therefore default set to $HOME directory
#         if "darwin" in sys.platform.lower():
#             return os.environ["HOME"]
#         else:
#             return tempfile.gettempdir()
        return os.getcwd()
    
    def clean(self, clean_all=True, clean_dir=False):
        for name in self.test_dict.keys():
            os.chdir(self.run_dir)
            logger.info("Cleaning {0} in directory {1}".format(name, self.run_dir))
            work_dir = os.path.join(self.run_dir, name)
            if os.path.isdir(work_dir): shutil.rmtree(work_dir)
            logfile = work_dir + '.log'
            if os.path.isfile(logfile): os.unlink(logfile)  
            if clean_all:
                script = work_dir + ample_util.SCRIPT_EXT
                if os.path.isfile(script): os.unlink(script)
        if clean_dir and os.path.isdir(self.run_dir): shutil.rmtree(self.run_dir)
    
    def run(self,
            nproc=1,
            dry_run=False,
            clean_up=True, 
            rosetta_dir=None,
            **kwargs):
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
        
        if dry_run: clean_up = False
        
        if rosetta_dir and not os.path.isdir(rosetta_dir):
            print "Cannot find rosetta_dir: {0}".format(rosetta_dir)
            sys.exit(1)
        
        if clean_up: self.clean()
            
        scripts = self._create_scripts(rosetta_dir, **kwargs)
        if not len(scripts):
            raise RuntimeError("Could not find any test cases to run!")
        
        logger.info("The following test cases will be run:")
        for name in self.test_dict.keys():
            logger.info("{0}: {1}".format(name, self.run_dir))
        
        ## Run all the jobs
        # If we're running on a cluster, we run on as many processors as there are jobs, 
        # as the jobs are just sitting and monitoring the queue
        if kwargs['submit_cluster']:
            logger.info("Jobs will be submitted to a cluster queueing system")
            nproc = len(scripts)
        
        if not dry_run:
            workers_util.run_scripts(job_scripts=scripts,
                                     monitor=None,
                                     chdir=True,
                                     nproc=nproc,
                                     job_name='test')
        
        # Now check the results using the unittesting framework
        self.run_unittest_suite()
        return 
    
    def _create_scripts(self, rosetta_dir, **kwargs):
        """Create scripts and path to resultsd"""
        scripts = []
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
            ensembler = True if name.startswith('ensembler') else False
            if ensembler and sys.platform.startswith('win'):
                logger.critical("Cannot run ensemble module on windows due to multiprocessing bug")
                continue
            
            script = self.write_script(work_dir,  args + [['-work_dir', work_dir]], ensembler=ensembler)
            scripts.append(script)
            # Set path to the results pkl file we will use to run the tests
            self.test_dict[name]['resultsd'] = os.path.join(work_dir,'resultsd.pkl')
            os.chdir(owd)            # Back to where we started
        return scripts
    
    def _is_in_args(self, argt, args):
        if type(argt) is str:
            key = argt
        else:
            key = argt[0]
        return key in [ a[0] for a in args ]

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
        if not kwargs['submit_cluster']: return args
        for k, v in kwargs.iteritems():
            value = None
            if k in CLUSTER_ARGS.keys():
                if v is not None:
                    value = v
                elif CLUSTER_ARGS[k] is not None:
                    value = CLUSTER_ARGS[k]
            if value:
                # Need to add the hypen on to the key so it can be used as a command-line arg
                args.append(["-"+k, value])
        return args
    
    def run_unittest_suite(self):
        suite = TestSuite()
        for name in self.test_dict.keys():
            testClass = self.test_dict[name]['test']
            testClass.RESULTS_PKL = self.test_dict[name]['resultsd']
            _suite = TestLoader().loadTestsFromTestCase(testClass)
            suite.addTests(_suite)  
        TextTestRunner(verbosity=2).run(suite)
    
    def write_script(self, work_dir, args, ensembler):
        """Write script"""
        linechar = "^" if sys.platform.startswith('win') else "\\"
        script = work_dir + ample_util.SCRIPT_EXT

        test_exe = os.path.join(os.environ["CCP4"], "bin", "ample")
        test_exe = test_exe + ample_util.SCRIPT_EXT if sys.platform.startswith("win") else test_exe
        if ensembler:
            if sys.platform.startswith("win"): raise RuntimeError("Cannot run ensemble module on windows due to multiprocessing bug")
            test_exe = '{0} -m ample.ensembler'.format(os.path.join(os.environ["CCP4"], "bin", "ccp4-python"))

        with open(script, 'w') as f:
            f.write(ample_util.SCRIPT_HEADER + os.linesep)
            f.write(os.linesep)
            f.write("{0} {1}".format(test_exe, linechar + os.linesep))
            for argt in args:
                f.write(" ".join(argt) + " " + linechar + os.linesep)
            f.write(os.linesep)
            f.write(os.linesep)
        os.chmod(script, 0o777)
        return os.path.abspath(script)


class SuiteLoader(object):
    """Loader designed to obtain all test cases in a package"""
    
    def load_cases(self, directory, test_cases=None, pattern="test_cases"):
        """function to load a integration test suite"""
        search_pattern = os.path.join(directory, "*")
        cases = [ os.path.basename(folder) for folder in \
                        glob.glob(search_pattern) if os.path.isdir(folder) ]
        test_dict = self._load_cases(cases, directory, pattern)
        # Needs to follow as case names will be folders and test_cases the 
        # actual cases themselves
        if test_cases:
            test_dict = {k:v for k,v in test_dict.iteritems() if k in test_cases}
        return test_dict
        
    def _load_cases(self, cases, directory, pattern):
        test_cases = {}
        for example_dir in cases:
            path = os.path.join(directory, example_dir)
            test_module = self.load_module(pattern, [path])
            # Skip anything that's not a valid AMPLE test module
            if not test_module or not hasattr(test_module, 'TEST_DICT'): continue
            for k, v in test_module.TEST_DICT.iteritems():
                if k in test_cases:
                    raise RuntimeError("Duplicate key: {0}".format(k))
                test_cases[k] = v
        return test_cases
    
    def load_module(self, mod_name, paths):
        try:
            mfile, pathname, desc = imp.find_module(mod_name, paths)
        except ImportError:
            logger.critical("Cannot find test module in {1}".format(mod_name, paths))
            return None
        try:
            test_module = imp.load_module(mod_name, mfile, pathname, desc)
        except Exception as e:
            logger.critical("Error loading test case from directory: {0}\n {1}\n".format(paths, e))
            raise Exception(e)
        finally:
            mfile.close()
        return test_module
