'''
Created on 29 Dec 2015

@author: jmht
'''
import argparse
import os
import shutil
import sys

# Our imports
from ample_util import SCRIPT_EXT, SCRIPT_HEADER
import workers

AMPLE_DIR = os.sep.join(os.path.abspath(os.path.dirname(__file__)).split(os.sep)[ :-1 ])

class AmpleException(Exception): pass

def write_script(name, args):
    ample = os.path.join(AMPLE_DIR,'bin', 'ample.py')
    script = name + SCRIPT_EXT
    with open(script, 'w') as f:
        f.write(SCRIPT_HEADER + os.linesep)
        f.write(os.linesep)
        f.write(ample + " \\" + os.linesep)
        # Assumption is all arguments are in pairs
        arg_list = [ " ".join(args[i:i+2]) for i in range(0, len(args), 2) ]
        f.write(" \\\n".join(arg_list))
        f.write(os.linesep)
        f.write(os.linesep)
    
    os.chmod(script, 0o777)
    return os.path.abspath(script)

def clean(test_dict):
    for name in test_dict.keys():
        if os.path.isdir(name): shutil.rmtree(name)
        logfile = name+'.log'
        if os.path.isfile(logfile): os.unlink(logfile)  
        script = name+SCRIPT_EXT
        if os.path.isfile(script): os.unlink(script)  

def run(test_dict,
        nproc=1,
        submit_cluster=False,
        dry_run=False,
        clean_up=True,
        rosetta_dir=None,
        extra_args=None,
        **kw):

    # Create scripts and path to resultsd
    scripts = []
    for name in test_dict.keys():
        args = test_dict[name]['args']
        # Rosetta is the only think likely to change between platforms so we update the entry
        if rosetta_dir:
            args = args['rosetta_dir'] = rosetta_dir
        if extra_args:
            args += extra_args
        script = write_script(name,  args + ['-work_dir', name])
        scripts.append(script)
        test_dict[name]['resultsd'] = os.path.join(name,'resultsd.pkl')
        if clean_up:
            if os.path.isdir(name): shutil.rmtree(name)
            logfile = name+'.log'
            if os.path.isdir(logfile): os.unlink(logfile)
    
    # Run all the jobs
    nproc = nproc
    submit_cluster = submit_cluster
    submit_qtype = 'SGE'
    submit_array = True
    if not dry_run:
        workers.run_scripts(job_scripts=scripts,
                            monitor=None,
                            chdir=True,
                            nproc=nproc,
                            job_time=3600,
                            job_name='test',
                            submit_cluster=submit_cluster,
                            submit_qtype=submit_qtype,
                            submit_queue=None,
                            submit_array=submit_array,
                            submit_max_array=None)
    
    # Now run the tests
    for name in test_dict.keys():
        try:
            test_dict[name]['test'](test_dict[name]['resultsd'])
            print "Job \'{0}\' succeeded".format(name)
        except AmpleException as ae:
            print "* Job \'{0}\' failed a test: {1}".format(name, ae)
        except Exception as e:
            print "*** Job \'{0}\' generated an exeption: {1}".format(name, e)


def parse_args(test_dict=None, extra_args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-clean', action='store_true', default=False,
                        help="Clean up all test files/directories")
    parser.add_argument('-nproc', type=int, default=1,
                        help="Number of processors to run on (1 per job)")
    parser.add_argument('-dry_run', action='store_true', default=False,
                        help="Don\'t actually run the jobs")
    parser.add_argument('-rosetta_dir',
                        help="Location of rosetta installation directory")
    parser.add_argument('-submit_cluster', action='store_true', default=False,
                        help="Submit to a cluster queueing system")
    
    args = parser.parse_args()
    if args.rosetta_dir and not os.path.isdir(args.rosetta_dir):
        print "Cannot find rosetta_dir: {0}".format(args.rosetta_dir)
        sys.exit(1)
    
    argd = vars(args)
    if test_dict:
        if args.clean:
            clean(test_dict)
        else:
            run(test_dict, extra_args=extra_args, **argd)
    else:
        return argd

