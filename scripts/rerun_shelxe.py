#!/usr/bin/env ccp4-python


'''
18.01.2017

@author: jmht
'''

# System
import argparse
import logging
import os
import cPickle
import shutil
import sys

# Horribleness to find MRBUMP
if not "CCP4" in sorted(os.environ.keys()):
    sys.stderr.write('CCP4 not found' + os.linesep)
    sys.exit(1)

# Update the python path with the mrbump folders
if os.path.isdir(os.path.join(os.environ["CCP4"], "share", "mrbump")):
    mrbump = os.path.join(os.environ["CCP4"], "share", "mrbump")
elif os.path.isdir(os.path.join(os.environ["MRBUMP"], "share", "mrbump")):
    mrbump = os.path.join(os.environ["MRBUMP"], "share", "mrbump")
else:
    sys.stderr.write("Error: MrBUMP installation not found" + os.linesep)
    sys.exit(1)

mrbump_incl = os.path.join(mrbump, "include")
sys.path.append(os.path.join(mrbump_incl, 'parsers'))
sys.path.append(os.path.join(mrbump_incl, 'building'))
import MRBUMP_Shelxe

# Custom
from ample.util import workers_util

sys.dont_write_bytecode = True
BK_SUFFIX = '.orig'

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# String for Python script to run shelxe and arpwarp/buccanner depending on SHELXE results
run_shelxe_script_str = """#!/usr/bin/env ccp4-python
import os
import subprocess
import sys

shelxe_script = "{0}"
arp_script = "{1}"
bucc_script = "{2}"

# Horribleness to find MRBUMP
if not "CCP4" in sorted(os.environ.keys()):
    sys.stderr.write('CCP4 not found' + os.linesep)
    sys.exit(1)

# Update the python path with the mrbump folders
if os.path.isdir(os.path.join(os.environ["CCP4"], "share", "mrbump")):
   mrbump = os.path.join(os.environ["CCP4"], "share", "mrbump")
elif os.path.isdir(os.path.join(os.environ["MRBUMP"], "share", "mrbump")):
   mrbump = os.path.join(os.environ["MRBUMP"], "share", "mrbump")
else:
    sys.stderr.write("Error: MrBUMP installation not found" + os.linesep)
    sys.exit(1)

mrbump_incl = os.path.join(mrbump, "include")
sys.path.append(os.path.join(mrbump_incl, 'parsers'))
sys.path.append(os.path.join(mrbump_incl, 'building'))
import MRBUMP_Shelxe

cmd = [ shelxe_script ]
wdir = os.path.dirname(shelxe_script)
shelxe_logfile = os.path.join(wdir, 'shelxe_run.log')
log = open(shelxe_logfile, 'w')
sys.stdout.write("Running SHELXE cmd: {{0}} Logfile is: {{1}}{{2}}".format(" ".join(cmd), shelxe_logfile, os.linesep))
p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=wdir)
rtn = p.wait()
log.close()

if rtn != 0:
    sys.stderr.write("SHELXE returned non-zero return code" + os.linesep)
    sys.exit(1)

if not (arp_script or bucc_script):
    # No rebuiling so we are finished
    sys.exit(0)

# We can do a rebuild - so first check SHELXE scores
shelxe_job = MRBUMP_Shelxe.Shelxe()
shelxe_job.shelxeLogfile = shelxe_logfile
if not shelxe_job.succeeded():
    sys.stdout.write("SHELXE did not succeed" + os.linesep)
    sys.exit(0)
else:
    sys.stdout.write("SHELXE WORKED" + os.linesep)
    

# Shelxe Worked so run arpwarp and then buccaneer
if arp_script:
    sys.stdout.write("Running ARPWARP" + os.linesep)
    cmd = [ arp_script ]
    wdir = os.path.dirname(arp_script)
    arp_logfile = os.path.join(wdir, 'arpwarp.log')
    log = open(arp_logfile, 'w')
    p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=wdir)
    rtn = p.wait()
    log.close()
    
if bucc_script:
    sys.stdout.write("Running BUCCANEER" + os.linesep)
    cmd = [ bucc_script ]
    wdir = os.path.dirname(bucc_script)
    bucc_logfile = os.path.join(wdir, 'buccaneer.log')
    log = open(bucc_logfile, 'w')
    p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=wdir)
    rtn = p.wait()
    log.close() 

sys.exit(0)
    
"""
#run_shelxe_script_str = """#!/bin/bash
#echo hello
#"""

def create_scripts(amoptd):
    # For each mrbump job in amoptd
    job_scripts = []
    for d in amoptd['mrbump_results']:
        #print sorted(d.keys())
        
        # Define the various directories         
        mr_dir = d['Job_directory']
        build_dir = os.path.join(mr_dir, 'build')
        if not os.path.isdir(build_dir): continue # Can't do owt if it doesn't exist
        
        # Backup old build directory and create new ones
        build_dir_bk = build_dir + BK_SUFFIX
        shutil.move(build_dir, build_dir_bk)
        
        # Create new directory tree
        shelxe_dir = os.path.join(build_dir, 'shelxe')
        arp_dir = os.path.join(shelxe_dir, 'rebuild', 'arpwarp')
        bucc_dir = os.path.join(shelxe_dir, 'rebuild', 'buccaneer')
        os.makedirs(arp_dir)
        os.makedirs(bucc_dir)
        
        shelxe_script_old = os.path.join(build_dir_bk, 'shelxe', 'shelx-script.sh')
        shelxe_script_new = os.path.join(build_dir, 'shelxe', 'shelx-script.sh')
        
        arp_script_old = os.path.join(build_dir_bk, 'shelxe', 'rebuild', 'arpwarp', 'arpwarp-script.sh')
        arp_script_new = os.path.join(build_dir, 'shelxe', 'rebuild', 'arpwarp', 'arpwarp-script.sh')
        
        bucc_script_old = os.path.join(build_dir_bk, 'shelxe', 'rebuild', 'buccaneer', 'buccaneer-script.sh')
        bucc_script_new = os.path.join(build_dir, 'shelxe', 'rebuild', 'buccaneer', 'buccaneer-script.sh')
        
        # Copy in the run scripts, amending if necessary - add any keywords required here for time being
        if not os.path.isfile(shelxe_script_old):
            raise RuntimeError("Cannot find shelxe_script: {0}".format(shelxe_script_old))
        kwargs = {
                  'shelxe_exe' : '/home/jmht/bin/shelxe.2014.4.george'
        }
        if 'native_pdb' in amoptd and amoptd['native_pdb'] and os.path.isfile(amoptd['native_pdb']):
            kwargs['native_pdb'] = amoptd['native_pdb']
        manipulate_shelxe_script(shelxe_script_old, shelxe_script_new, **kwargs)
        
        if os.path.isfile(arp_script_old):
            manipulate_arp_script(arp_script_old, arp_script_new)
        else:
            arp_script_new = None

        if os.path.isfile(bucc_script_old):
            manipulate_bucc_script(bucc_script_old, bucc_script_new)
        else:
            bucc_script_new = None

        # Create script to run the steps for this job
        run_script = os.path.join(build_dir, "run_{0}.py".format(d['name']))
        with open(run_script, 'w') as w:
            w.write(run_shelxe_script_str.format(shelxe_script_new, arp_script_new, bucc_script_new))
        os.chmod(run_script, 0o777)
        
        # Add script to run as bash script so that we can submit to the queieing system? - TEST IF NEEDED
        # Add to list of scripts
        job_scripts.append(run_script)
        break
    return job_scripts

def manipulate_shelxe_script(old_script, new_script, **kwargs):
    """Copy old_script to new_script, if necessary updating according to given kwargs"""
    if len(kwargs) == 0:
        shutil.copy2(old_script, new_script)
    else:
        cmd = "#!/bin/bash\n\n"
        with open(old_script, 'r') as f:
            for line in iter(f.readline, ''):
                if line[:2] != 'cp' and 'shelxe-input.pda' in line:
                    l = line.strip().split()
                    if 'shelxe_exe' in kwargs:
                        l[0] = kwargs['shelxe_exe']
                    if 'n_cycles' in kwargs:
                        pass
                    line = " ".join(l + [os.linesep])
                    if 'native_pdb' in kwargs:
                        # Hack in copy of native before command to run shelxe
                        assert os.path.isfile(kwargs['native_pdb'])
                        cline = "cp {0} shelxe-input.ent\n\n".format(kwargs['native_pdb'])
                        line = cline + line
                cmd += line
        with open(new_script, 'w') as o: o.write(cmd)
        os.chmod(new_script, 0o777)
    return new_script

def manipulate_arp_script(old_script, new_script, **kwargs):
    """Copy old_script to new_script, if necessary updating according to given kwargs"""
    if len(kwargs) == 0:
        shutil.copy2(old_script, new_script)
    else:
        raise NotImplementedError()
    
def manipulate_bucc_script(old_script, new_script, **kwargs):
    """Copy old_script to new_script, if necessary updating according to given kwargs"""
    if len(kwargs) == 0:
        shutil.copy2(old_script, new_script)
    else:
        raise NotImplementedError()
    
def shelxe_results(shelxe_logfile, d):
    shelxe_job = MRBUMP_Shelxe.Shelxe()
    shelxe_job.shelxeLogfile = shelxe_logfile
    shelxe_job.parseLog()
    dkeys = [ 'SHELXE_CC', 'SHELXE_ACL', 'SHELXE_MCL', 'SHELXE_NC', 'SHELXE_time', 'SHELXE_version',
             'SHELXE_wMPE', 'SHELXE_os' ]
    for k in dkeys: d[k] = shelxe_job.resultsDict[k]
    return d

def arp_results(bucc_logfile, d):
    ap = parse_arpwarp.ArpwarpLogParser(bucc_logfile)
    d["SXRARP_version"] = ap.version
    d["SXRARP_final_Rfact"] = ap.finalRfact
    d["SXRARP_final_Rfree"] = ap.finalRfree
    return d

def bucc_results(bucc_logfile, d):
    bp = parse_buccaneer.BuccaneerLogParser(bucc_logfile)
    d["SXRBUCC_version"] = bp.version
    d["SXRBUCC_final_Rfact"] = bp.finalRfact
    d["SXRBUCC_final_Rfree"] = bp.finalRfree
    return d

def rerun_shelxe(args):
    logger.info('Preparing scripts')

    # Unpickle dictionary
    amopt_pkl = args.ample_pkl
    with open(args.ample_pkl) as f: amoptd = cPickle.load(f)
    assert 'mrbump_results' in amoptd, "No MRBUMP results in: %s" % amopt_pkl
    
    # Back up old AMPLE pkl file - preserve metadata
    shutil.copy2(amopt_pkl, amopt_pkl + BK_SUFFIX)

    # Get list of jobs to rerun the SHELXE pipeline
    job_scripts = create_scripts(amoptd)
    
    # Run the jobs
    if True:
        logger.info("Running scripts:\n{0}".format(os.linesep.join(job_scripts)))
        ok = workers_util.run_scripts(job_scripts=job_scripts,
                                      nproc=args.nproc,
                                      chdir=True,
                                      submit_cluster=args.submit_cluster,
                                      submit_qtype="SGE",
                                      submit_queue="all.q",
                                      submit_array=True,
                                      submit_max_array=10)
    
    # Collect results from completed jobs
    for oldd in amoptd['mrbump_results']:
        #print sorted(d.keys())
        
        # Add SHELXE, ARPWARP and BUCCANEER results to amoptd
        # We can use the path to the old log files as the paths should be the same
        newd = {}
        if oldd['SHELXE_logfile'] and os.path.isfile(oldd['SHELXE_logfile']):
            newd = shelxe_results(oldd['SHELXE_logfile'], newd)
        if oldd['SXRARP_logfile'] and os.path.isfile(oldd['SXRARP_logfile']):
            newd = arp_results(oldd['SXRARP_logfile'], newd)
        if oldd['SXRARP_logfile'] and os.path.isfile(oldd['SXRBUCC_logfile']):
            newd = arp_results(oldd['SXRBUCC_logfile'], newd)
        
        # Update AMPLE and MRBUMP dictionaries with new values
        mrb_pkl = os.path.join(oldd['Search_directory'], 'results', 'resultsTable.pkl')
        shutil.copy2(mrb_pkl, mrb_pkl + BK_SUFFIX) # Backup old MRBUMP results
        with open(mrb_pkl) as w:  mrb_dict = cPickle.load(w)
       
        # Update values in dictionaries
        for k in newd.keys():
            oldd[k] = newd[k]
            mrb_dict[oldd['name']][k] = newd[k]

        # Writ out updated mrbump dict
        with open(mrb_pkl, 'w') as w: cPickle.dump(mrb_dict,w)
        break
    
    # Write out the updated amoptd
    with open(amopt_pkl, 'w') as w: cPickle.dump(amoptd,w)

    return job_scripts
##End rerun_shelxe()

def main():
    p = argparse.ArgumentParser(description='Rerun SHELXE', prefix_chars="-")
    p.add_argument('-nproc', default=1, type=int,
                    help='Number of processors to run jobs in parallel')
    p.add_argument('-submit_cluster', action="store_true")
    #p.add_argument('-rebuild', action="store_true", help="Rerun rebuild of SHELXE trace with ArpWARP/BUCCANEER")
    p.add_argument('ample_pkl', help='AMPLE pkl file')
    args=p.parse_args()

    args.ample_pkl = os.path.abspath(args.ample_pkl)
    assert os.path.isfile(args.ample_pkl), "Cannot find pkl file: %s" % args.ample_pkl
    logger.info("Re-running SHELXE using:\n%s" % args.ample_pkl)
    
    # Get all search models based on the MrBump run scripts
    #search_models = [re.findall(r'MRBUMP/(.*).sh', x)[0] for x in glob.glob(os.path.join(mrbump_d,'*.sh'))]
    #logger.info("Imported %d search models" % len(search_models))
    rerun_shelxe(args)

if __name__ == "__main__":
    main()
