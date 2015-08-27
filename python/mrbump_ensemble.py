'''
Created on Feb 28, 2013

@author: jmht
'''

# python imports
import logging
import os
import sys
import unittest

# our imports
import mrbump_cmd


"""
for ensemble in ensemble_pdbs:
    create dictionary with all options
    write_keyword_file
    write_script_file
"""

def write_mrbump_files(ensemble_pdbs, amoptd, job_time=86400, ensemble_options=None, directory=None):
    """Write the MRBUMP job files for all the ensembles.

    Args:
    ensemble_pdbs -- list of the ensembles, each a single pdb file
    amoptd -- dictionary with job options
    """
    if not directory: directory = os.getcwd()
    
    job_scripts = []
    keyword_options = {}
    for ensemble_pdb in ensemble_pdbs:
        name = os.path.splitext(os.path.basename(ensemble_pdb))[0] # Get name from pdb path
        
        # Get any options specific to this ensemble
        if ensemble_options and name in ensemble_options: keyword_options = ensemble_options[name]
        
        # Generate dictionary with all the options for this job and write to keyword file
        keyword_dict = mrbump_cmd.keyword_dict(ensemble_pdb, name, amoptd, keyword_options)
        keyword_file = os.path.join(directory,name+'.mrbump')
        keyword_str = mrbump_cmd.mrbump_keyword_file(keyword_dict)
        with open(keyword_file,'w') as f: f.write(keyword_str)
        
        script = write_jobscript(name,
                                 keyword_file,
                                 amoptd,
                                 job_time = job_time
                                 )
        job_scripts.append(script)
            
    if not len(job_scripts):
        msg = "No job scripts created!"
        logging.critical(msg)
        raise RuntimeError, msg
    
    return job_scripts

def write_jobscript(name, keyword_file, amoptd, directory=None, job_time=86400, extra_options={}):
    """
    Create the script to run MrBump for this PDB.
    """
    if not directory: directory = os.getcwd()
        
    # Next the script to run mrbump
    ext='.bat' if sys.platform.startswith("win") else '.sh'
    script_path = os.path.join(directory,name+ext)
    with open(script_path, "w") as job_script:
        # Header
        if not sys.platform.startswith("win"):
            script_header = '#!/bin/sh\n'
            script_header += '[[ ! -d $CCP4_SCR ]] && mkdir $CCP4_SCR\n\n'
            job_script.write(script_header)
        
        # Get the mrbump command-line
        jobcmd = mrbump_cmd.mrbump_cmd(name, amoptd['mtz'], amoptd['mr_sequence'], keyword_file)
        job_script.write(jobcmd)
        
    # Make executable
    os.chmod(script_path, 0o777)
    
    return script_path

class Test(unittest.TestCase):

    def XtestLocal(self):
        
        import glob
        
        d = {
             'mtz' : '/opt/ample-dev1/examples/toxd-example/1dtx.mtz',
             'fasta' : '/opt/ample-dev1/examples/toxd-example/toxd_.fasta',
             'mrbump_programs' : ' molrep ',
             'use_buccaneer' : False ,
             'buccaneer_cycles' : 5,
             'use_arpwarp' : False,
             'arpwarp_cycles': 10,
             'use_shelxe' : True,
             'shelx_cycles' : 10,
             'FREE' : 'FreeR_flag',
             'F' : 'FP',
             'SIGF' : 'SIGFP',
             'nproc' : 3,
             'domain_all_chains_pdb' : None,
             'ASU': None,
             'mr_keys': [],
             'early_terminate' : False
             }
        
        ensemble_dir = "/home/Shared/ample-dev1/examples/toxd-example/ROSETTA_MR_0/ensembles_1"
        ensembles = []
        for infile in glob.glob( os.path.join( ensemble_dir, '*.pdb' ) ):
            ensembles.append(infile)
            
        mrbump_ensemble_local( ensembles, d )
            
    
    def XtestSuccess(self):
        
        pdbid = "/home/jmht/t/ample-dev1/examples/toxd-example/ROSETTA_MR_3/MRBUMP/cluster_2/search_All_atom_trunc_0.551637_rad_3_mrbump"
        
        self.assertTrue( check_success( pdbid ) )
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    #unittest.main()
    # Nothing here - see notes in worker
    pass
    
