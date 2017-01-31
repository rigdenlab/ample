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
sys.path.append(os.path.join(mrbump_incl, 'building'))
sys.path.append(os.path.join(mrbump_incl, 'parsers'))
import MRBUMP_ARPwARP
import MRBUMP_Buccaneer
import MRBUMP_Shelxe
import parse_arpwarp
import parse_buccaneer

# Custom
from ample.util import workers_util

sys.dont_write_bytecode = True
BK_SUFFIX = '.orig'

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# String for Python script to run shelxe and arpwarp/buccanner depending on SHELXE results
run_shelxe_script_str = """#!/usr/bin/env ccp4-python
import logging
import os
import subprocess
import sys

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

shelxe_script = "{shelxe_script}"
arp_script = "{arp_script}"
if arp_script == "None": arp_script = None
bucc_script = "{bucc_script}"
if bucc_script == "None": bucc_script = None

# Horribleness to find MRBUMP
if not "CCP4" in sorted(os.environ.keys()):
    logger.critical('CCP4 not found' + os.linesep)
    sys.exit(1)

# Update the python path with the mrbump folders
if os.path.isdir(os.path.join(os.environ["CCP4"], "share", "mrbump")):
   mrbump = os.path.join(os.environ["CCP4"], "share", "mrbump")
elif os.path.isdir(os.path.join(os.environ["MRBUMP"], "share", "mrbump")):
   mrbump = os.path.join(os.environ["MRBUMP"], "share", "mrbump")
else:
    logger.critical("Error: MrBUMP installation not found" + os.linesep)
    sys.exit(1)

mrbump_incl = os.path.join(mrbump, "include")
sys.path.append(os.path.join(mrbump_incl, 'building'))
sys.path.append(os.path.join(mrbump_incl, 'parsers'))
sys.path.append(os.path.join(mrbump_incl, 'tools'))
import MRBUMP_Shelxe
import MRBUMP_phs2mtz
import MRBUMP_ARPwARP
import MRBUMP_Buccaneer

cmd = [ shelxe_script ]
shelxe_dir = os.path.dirname(shelxe_script)
os.chdir(shelxe_dir)
shelxe_logfile = os.path.join(shelxe_dir, 'shelxe_run.log')
log = open(shelxe_logfile, 'w')
logger.info("Running SHELXE cmd: {{0}} Logfile is: {{1}}{{2}}".format(" ".join(cmd), shelxe_logfile, os.linesep))
p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=shelxe_dir)
rtn = p.wait()
log.close()

if rtn != 0:
    logger.critical("SHELXE returned non-zero return code" + os.linesep)
    sys.exit(1)

# Check SHELXE scores
if not os.path.isfile(shelxe_logfile):
    logger.critical("No SHELXE logfile: {{0}}{{1}}".format(shelxe_logfile,os.linesep))
    sys.exit(1)

shelxe_job = MRBUMP_Shelxe.Shelxe()
shelxe_job.shelxeLogfile = shelxe_logfile
if not shelxe_job.succeeded():
    logger.info("SHELXE did not succeed" + os.linesep)
    sys.exit(0)
else:
    logger.info("SHELXE WORKED" + os.linesep)
    
# Convert the phs file to an mtz
phsin = "{shelxe_phs}"
pdbin = "{shelxe_pdb}"
if os.path.isfile(phsin) and os.path.isfile(pdbin):
    sys.stdout.write("Running phs2mtz" + os.linesep)
    p2m = MRBUMP_phs2mtz.PHS2MTZ()
    p2m.phs2mtz(phsin,
                pdbin,
                "{shelxe_mtz}",
                "{shelxe_dir}",
                hklref="{hklref}",
                freeLabel="{FREE}",
                resolution={resolution}
                )

rebuild_dir = os.path.join(shelxe_dir,'rebuild')
if not os.path.isdir(rebuild_dir): os.mkdir(rebuild_dir)

# Shelxe Worked so run arpwarp and then buccaneer
if arp_script:
    logger.info("Running ARPWARP from existing script" + os.linesep)
    cmd = [ arp_script ]
    wdir = os.path.dirname(arp_script)
    os.chdir(wdir)
    arp_logfile = os.path.join(wdir, 'arpwarp.log')
    log = open(arp_logfile, 'w')
    p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=wdir)
    rtn = p.wait()
    log.close()
else:
    logger.info("Running ARPWARP from scratch" + os.linesep)
    arpwarp_dir = os.path.join(rebuild_dir, "arpwarp")
    if not os.path.isdir(arpwarp_dir): os.mkdir(arpwarp_dir)
    arpwarpSX = MRBUMP_ARPwARP.Arpwarp()
    arpwarpSX.runARPwARP("{fasta}",
                         "{shelxe_mtz}",
                         arpwarp_dir,
                         'F',
                         'SIGF',
                         '{FREE}',
                         "PHI_SHELXE",
                         "FOM_SHELXE",
                         '{nmasu}',
                         cycles=5,
                         pdbinFile="{shelxe_pdb}")
    
    
if bucc_script:
    logger.info("Running BUCCANEER from existing script" + os.linesep)
    cmd = [ bucc_script ]
    wdir = os.path.dirname(bucc_script)
    os.chdir(wdir)
    bucc_logfile = os.path.join(wdir, 'buccaneer.log')
    log = open(bucc_logfile, 'w')
    p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=wdir)
    rtn = p.wait()
    log.close()
else:
    logger.info("Running BUCCANEER from scratch" + os.linesep)
    bucc_dir = os.path.join(rebuild_dir, "buccaneer")
    if not os.path.isdir(bucc_dir): os.mkdir(bucc_dir)
    buccaneerSX = MRBUMP_Buccaneer.Buccaneer()
    buccaneerSX.runBuccaneer("{fasta}",
                             "{shelxe_mtz}",
                             os.path.join(bucc_dir, "buccSX_output.pdb"),
                             bucc_dir,
                             'F',
                             'SIGF',
                             '{FREE}',
                             "PHI_SHELXE",
                             "FOM_SHELXE",
                             cycles=5,
                             pdbinFile="{shelxe_pdb}")

sys.exit(0)
    
"""
#run_shelxe_script_str = """#!/bin/bash
#echo hello
#"""

def get_nmasu(amoptd):
    """Parse mrbump logfile to get nmasu"""
    d = amoptd['mrbump_results'][0]
    logfile = os.path.join(amoptd['mrbump_dir'],d['ensemble_name'] + '.log')
    with open(logfile) as f:
        for line in f:
            if line.startswith('Estimated number of molecules to search for in a.s.u.:'):
                nmasu = line.split()[9]
                logger.info('Got nmasu: {0}'.format(nmasu))
                return nmasu
    assert False

def create_scripts(amoptd, args):
    
    nmasu = get_nmasu(amoptd)
    
    # For each mrbump job in amoptd
    job_scripts = []
    for d in amoptd['mrbump_results']:
        #print sorted(d.keys())
        
        # Define the various directories
        if 'MR_directory' in d:  
            mr_dir = d['MR_directory']
        elif 'Job_directory' in d:
            mr_dir = d['Job_directory']
        else: assert "No valid MR_directory key in dict!"
            
        build_dir = os.path.join(mr_dir, 'build')
        if not os.path.isdir(build_dir): continue # Can't do owt if it doesn't exist
        logger.info("Processing directory: {0}".format(mr_dir))
        
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
        
        # Need to set arp and bucc logfile paths if they haven't already been set or we can't parse any newly generated results
        d['SXRARP_logfile'] =  os.path.join(arp_dir,"arpwarp.log")
        d['SXRBUCC_logfile'] =  os.path.join(bucc_dir,"buccaneer.log")
        
        
        # Copy in the run scripts, amending if necessary - add any keywords required here for time being
        if not os.path.isfile(shelxe_script_old):
            #raise RuntimeError("Cannot find shelxe_script: {0}".format(shelxe_script_old))
            continue
        kwargs = {
                  'shelxe_exe' : '/home/jmht/bin/shelxe.2014.4.george'
        }
        if 'native_pdb' in amoptd and amoptd['native_pdb'] and os.path.isfile(amoptd['native_pdb']):
            kwargs['native_pdb'] = amoptd['native_pdb']
        manipulate_shelxe_script(shelxe_script_old, shelxe_script_new, 'shelxe_phaser_{0}'.format(d['name']), **kwargs)
        
        if os.path.isfile(arp_script_old):
            manipulate_arp_script(arp_script_old, arp_script_new)
        else:
            arp_script_new = None

        if os.path.isfile(bucc_script_old):
            manipulate_bucc_script(bucc_script_old, bucc_script_new)
        else:
            bucc_script_new = None
            
        # Get required data for phs2mtz and to run shelxe/arp/bucc
        #base_name = 'shelxe_phaser_loc0_ALL_{0}_UNMOD'.format(d['name'])
        base_name = 'shelxe_phaser_{0}'.format(d['name'])
        shelxd = { 'shelxe_phs' : os.path.join(shelxe_dir,base_name + '.phs'),
                   'shelxe_pdb' : os.path.join(shelxe_dir,base_name + '.pdb'),
                   'shelxe_mtz' : os.path.join(shelxe_dir,base_name + '.mtz'),
                   'hklref' : amoptd['mtz'],
                   #'resolution' : "{0:1.2F}".format(amoptd['mtz_min_resolution']),
                   'resolution' : amoptd['mtz_min_resolution'],
                   'fasta' : amoptd['fasta'],
                   'F' : amoptd['F'],
                   'SIGF' : amoptd['SIGF'],
                   'FREE' : amoptd['FREE'],
                   'shelxe_dir' : shelxe_dir,
                   'shelxe_script' : shelxe_script_new,
                   'arp_script' : arp_script_new,
                   'bucc_script' : bucc_script_new,
                   'nmasu' : nmasu
                  }
        # Create script to run the steps for this job
        run_script = os.path.join(build_dir, "run_{0}.py".format(d['name']))
        with open(run_script, 'w') as w:
            w.write(run_shelxe_script_str.format(**shelxd))
        os.chmod(run_script, 0o777)
        
        # Horrible - the submission script needs to be a shell script, so we create a wrapper bash script    
        if args.submit_cluster:
            run_script_sh = os.path.join(build_dir, "run_{0}.sh".format(d['name']))
            with open(run_script_sh, 'w') as w:
                w.write("#!/bin/bash\n\n{0}\n".format(run_script))
            os.chmod(run_script_sh, 0o777)
            run_script = run_script_sh
        
        # Add script to run as bash script so that we can submit to the queieing system? - TEST IF NEEDED
        # Add to list of scripts
        job_scripts.append(run_script)
        #break
    return job_scripts

def manipulate_shelxe_script(old_script, new_script, sbase, **kwargs):
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
                # Old job failed, this one might not...
                if line.startswith('# No output from SHELXE'):
                    nline = "[[ -f shelxe-input.pdb ]] && mv shelxe-input.pdb {0}.pdb\n".format(sbase)
                    nline += "[[ -f shelxe-input.phs ]] && mv shelxe-input.phs {0}.phs\n\n".format(sbase)
                    line = nline

                cmd += line
        with open(new_script, 'w') as o: o.write(cmd)
        os.chmod(new_script, 0o777)
    return new_script

def manipulate_arp_script(old_script, new_script, **kwargs):
    """Copy old_script to new_script, if necessary updating according to given kwargs"""
    shebang = '#!/bin/bash\n'
    with open(old_script) as f: s = f.read()
    s = shebang + s
    with open(new_script,'w') as w: w.write(s)
    #if len(kwargs) == 0:
    #   shutil.copy2(old_script, new_script)
    #else:
    #    raise NotImplementedError()
    os.chmod(new_script, 0o777)
    
def manipulate_bucc_script(old_script, new_script, **kwargs):
    """Copy old_script to new_script, if necessary updating according to given kwargs"""
    shebang = '#!/bin/bash\n'
    with open(old_script) as f: s = f.read()
    s = shebang + s
    with open(new_script,'w') as w: w.write(s)
    #if len(kwargs) == 0:
    #    shutil.copy2(old_script, new_script)
    #else:
    #    raise NotImplementedError()
    os.chmod(new_script, 0o777)
    
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
    #assert not os.path.isfile(amopt_pkl + BK_SUFFIX)
    shutil.copy2(amopt_pkl, amopt_pkl + BK_SUFFIX)

    if True:
        # Get list of jobs to rerun the SHELXE pipeline
        job_scripts = create_scripts(amoptd, args)

        # Run the jobs
        logger.info("Running scripts:\n{0}".format(os.linesep.join(job_scripts)))
        ok = workers_util.run_scripts(job_scripts=job_scripts,
                                      nproc=args.nproc,
                                      submit_cluster=args.submit_cluster,
                                      submit_qtype="SGE",
                                      submit_queue="all.q",
                                      submit_array=True,
                                      submit_max_array=10)
    
    # Collect results from completed jobs
    for i, oldd in enumerate(amoptd['mrbump_results']):
        #if i == 1: continue
        #print("CHECKING ",oldd['Search_directory'])
        #print sorted(d.keys())
        
        # Add SHELXE, ARPWARP and BUCCANEER results to a dictionary
        # We need to manually set the path to the arp and bucc logfiles as 
        # they may not have been set in the previous run.
        
        newd = {}
        if oldd['SHELXE_logfile'] and os.path.isfile(oldd['SHELXE_logfile']):
            newd = shelxe_results(oldd['SHELXE_logfile'], newd)
        if oldd['SXRARP_logfile'] and os.path.isfile(oldd['SXRARP_logfile']):
            newd = arp_results(oldd['SXRARP_logfile'], newd)
        if oldd['SXRBUCC_logfile'] and os.path.isfile(oldd['SXRBUCC_logfile']):
            newd = bucc_results(oldd['SXRBUCC_logfile'], newd)
        
        # Update AMPLE and MRBUMP dictionaries with new values
        mrb_pkl = os.path.join(oldd['Search_directory'], 'results', 'resultsTable.pkl')
        assert not os.path.isfile(mrb_pkl + BK_SUFFIX)
        #shutil.copy2(mrb_pkl, mrb_pkl + BK_SUFFIX) # Backup old MRBUMP results
        with open(mrb_pkl) as w:  mrb_dict = cPickle.load(w)
       
        # Update values in dictionaries
        for k in newd.keys():
            oldd[k] = newd[k]
            mrb_dict[oldd['name']][oldd['MR_program']][k] = newd[k]

        # Writ out updated mrbump dict
        with open(mrb_pkl, 'w') as w: cPickle.dump(mrb_dict,w)
        #break
    
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
