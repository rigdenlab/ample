
import cPickle
import glob
import imp
import os
import shutil
import sys
import tempfile

from ample.constants import SHARE_DIR
from ample.util.ample_util import SCRIPT_EXT, SCRIPT_HEADER
from ample.util import workers_util
from ample.testing.constants import CLUSTER_ARGS, EXTRA_ARGS
from unittest import TestCase, TestLoader, TextTestRunner, TestSuite

__author__ = "Felix Simkovic and Jens Thomas"
__date__ = "25-Mar-2016"

# Available packages. Hard-coded for now to show visually what we have in
# argparse module. Not needed otherwise
PACKAGES = ['from_existing_models', 'from_quark_models', 'from_single_model', 
            'homologs', 'ideal_helices',  'import_cluster', 'import_ensembles', 
            'import_models', 'missing_domain', 'nmr_truncate']
if not sys.platform.startswith("win"):
    PACKAGES += ['nmr_remodel', 'rosetta_contacts', 'rosetta_modelling', 
                 'rosetta_restraints']

def add_cmd_options(parser):
    parser.add_argument('-c', '--clean', action='store_true', default=False,
                        help="Clean up all test files/directories")
    parser.add_argument('-n', '--nproc', type=int, default=1,
                        help="Number of processors to run on (1 per job)")
    parser.add_argument('-d', '--dry_run', action='store_true', default=False,
                        help="Don\'t actually run the jobs")
    parser.add_argument('-r', '--rosetta_dir',
                        help="Location of rosetta installation directory")
    parser.add_argument('-s', '--submit_cluster', action='store_true', default=False,
                        help="Submit to a cluster queueing system")
    parser.add_argument('test_cases', nargs='*',
                        help="[ {0} ]".format(" | ".join(PACKAGES)))
    parser.add_argument('-w', '--run_dir', type=str, default=None,
                        help="directory to run jobs in")

class AMPLEBaseTest(TestCase):
    RESULTS_PKL = None
    AMPLE_DICT = None
    def setUp(self):
        self.assertTrue(os.path.isfile(self.RESULTS_PKL),"Missing pkl file: {0}".format(self.RESULTS_PKL))
        with open(self.RESULTS_PKL) as f: self.AMPLE_DICT = cPickle.load(f)

class AMPLEIntegrationFramework(object):
    """Framework to run Ample integration testing"""
    
    def __init__(self, test_cases=None, run_dir=None):
        self.test_dict = SuiteLoader().load_cases(os.path.join(SHARE_DIR, "examples"), 
                                                  test_cases=test_cases)
        # Make a directory to keep all files together
        _root = os.path.abspath(run_dir) if run_dir else self.get_run_dir()
        self.run_dir = os.path.join(_root, "ample_testing")
        if not os.path.isdir(self.run_dir): os.mkdir(self.run_dir)
    
    def get_run_dir(self):
        # OS X has problems with relative paths - symlink in root messes
        # up the search in theseus - therefore default set to $HOME directory
        if "darwin" in sys.platform.lower():
            return os.environ["HOME"]
        else:
            return tempfile.gettempdir()
    
    def clean(self, clean_all=True, clean_dir=False):
        for name in self.test_dict.keys():
            os.chdir(self.run_dir)
            print "Cleaning {0} in directory {1}".format(name, self.run_dir)
            work_dir = os.path.join(self.run_dir, name)
            if os.path.isdir(work_dir): shutil.rmtree(work_dir)
            logfile = work_dir + '.log'
            if os.path.isfile(logfile): os.unlink(logfile)  
            if clean_all:
                script = work_dir + SCRIPT_EXT
                if os.path.isfile(script): os.unlink(script)
        if clean_dir and os.path.isdir(self.run_dir): shutil.rmtree(self.run_dir)
    
    def run(self, nproc=1, submit_cluster=False, dry_run=False, clean_up=True, 
            rosetta_dir=None, extra_args=None, **kw):
        
        print "Writing files to: {0}".format(self.run_dir)
        
        if dry_run: 
            clean_up = False
        
        if rosetta_dir and not os.path.isdir(rosetta_dir):
            print "Cannot find rosetta_dir: {0}".format(rosetta_dir)
            sys.exit(1)
        
        if clean_up: 
            self.clean()
        
        scripts = self._create_scripts(rosetta_dir, submit_cluster)
        
        print "The following test cases will be run:"
        for name in self.test_dict.keys():
            print "{0}: {1}".format(name, self.run_dir )
        
        ## Run all the jobs
        # If we're running on a cluster, we run on as many processors as there are jobs, 
        # as the jobs are just sitting and monitoring the queue
        if submit_cluster:
            nproc = len(scripts)
        
        if not dry_run:
            workers_util.run_scripts(job_scripts=scripts,
                                     monitor=None,
                                     chdir=True,
                                     nproc=nproc,
                                     job_name='test')

        self.run_unittest_suite()
        return 
    
    def _create_scripts(self, rosetta_dir, submit_cluster):
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
            if EXTRA_ARGS:
                args = self._update_args(args, EXTRA_ARGS)
            if submit_cluster:
                args = self._update_args(args, CLUSTER_ARGS)
            script = self.write_script(work_dir,  args + [['-work_dir', work_dir]])
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
    
    def run_unittest_suite(self):
        suite = TestSuite()
        for name in self.test_dict.keys():
            testClass = self.test_dict[name]['test']
            testClass.RESULTS_PKL = self.test_dict[name]['resultsd']
            _suite = TestLoader().loadTestsFromTestCase(testClass)
            suite.addTests(_suite)  
        TextTestRunner(verbosity=2).run(suite)
    
    def write_script(self, path, args):
        """Write script - ARGS MUST BE IN PAIRS"""
        linechar = "^" if sys.platform.startswith('win') else "\\"
        script = path + SCRIPT_EXT
        with open(script, 'w') as f:
            f.write(SCRIPT_HEADER + os.linesep)
            f.write(os.linesep)
            f.write("ccp4-python -m ample " + linechar + os.linesep)
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
            if not test_module: continue
            for k, v in test_module.TEST_DICT.iteritems():
                if k in test_cases:
                    raise RuntimeError("Duplicate key: {0}".format(k))
                test_cases[k] = v
        return test_cases
    
    def load_module(self, mod_name, paths):
        try:
            mfile, pathname, desc = imp.find_module(mod_name, paths)
        except ImportError:
            print "Cannot find test module in {1}".format(mod_name, paths)
            return None
        try:
            test_module = imp.load_module(mod_name, mfile, pathname, desc)
        finally:
            mfile.close()
        return test_module

        
