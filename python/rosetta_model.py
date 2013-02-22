'''
Created on 21 Feb 2013

@author: jmht
'''

# Python modules
import glob
import logging
import os
import random
import re
import shutil
import subprocess
import time
import unittest

# Our modules
import add_sidechains_SCWRL
import ample_options


class RosettaModel(object):
    """
    Class to run Rosetta modelling
    """

    def __init__(self):
        
        self.nproc = None
        self.nmodels = None
        self.work_dir = None
        self.models_dir = None
        self.rosetta_path = None
        self.rosetta_db = None
        self.frags3mers = None
        self.frags9mers = None
        self.fasta =None
        self.allatom = None
        self.use_scwrl = None
        self.scwrl_exe = None
        
        # Extra options
        self.domain_termini_distance = None
        self.rad_gyr_reweight = None
        self.improve_template = None
        #CCline
        
        self.logger = logging.getLogger()
    
    def generate_seeds(self):
        """
        Generate the list of seeds and return
        """
        
        previous_seeds = [0]  #make random seeds  (1000000, 6000000) ! Must be unique seeds!!!!
        proc = 1
        while proc < self.nproc + 1:
            new_seed = random.randint(1000000, 4000000)
            seed_present = False

            for seed in previous_seeds:
                if new_seed == seed:
                    seed_present = True
                    break

            if seed_present == False:
                previous_seeds.append(new_seed)
                proc += 1
                
        # Keep a log of the seeds
        seedlog = open( self.work_dir + os.sep +'seedlist', "w")
        for seed in  previous_seeds:
            seedlog.write(str(seed) + '\n')
        seedlog.close()
                
        return previous_seeds
    ##End generate_seeds
    
    def split_jobs(self):
        """
        Return a list of number of jobs to run on each processor
        """
        split_jobs = self.nmodels / self.nproc  # ## split jobs between processors
        remainder = self.nmodels % self.nproc
        jobs = [0]
        proc = 1
        while proc < self.nproc + 1:
            jobs.insert(proc, split_jobs)
            proc += 1
        jobs[-1] = jobs[-1] + remainder
        
        return jobs
    ##End split_jobs
    
    def rosetta_cmd(self, wdir, nstruct, seed):
        """
        Return the command to run rosetta as a list suitable for subprocess
        wdir: directory to run in
        nstruct: number of structures to process
        seed: seed for this processor
        """
        
        cmd = [ self.rosetta_path,
               '-database', self.rosetta_db,
               '-in::file::fasta', self.fasta,
               '-in:file:frag3', self.frags3mers,
               '-in:file:frag9', self.frags9mers,
               '-out:path', wdir,
               '-out:pdb',
               '-out:nstruct', nstruct,
               '-out:file:silent', wdir + '/OUT',
               '-run:constant_seed',
               '-run:jran', seed 
                ]
            
        if self.allatom:
            cmd += [ '-return_full_atom true', '-abinitio:relax' ]
        else:
            cmd += [ '-return_full_atom false' ]
            
        # Domain constraints
        if self.domain_termini_distance > 0:
            dcmd = self.setup_domain_constraints()
            cmd += dcmd
            
        # Radius of gyration reweight
        if self.rad_gyr_reweight:
            if "none" in self.rad_gyr_reweight.lower():
                cmd+= ['-rg_reweight', '0']
                
        # Improve Template
        if self.improve_template:
            tcmd = ['-in:file:native',
                    self.improve_template,
                    '-abinitio:steal_3mers',
                    'True',
                    '-abinitio:steal9mers',
                    'True',
                    '-abinitio:start_native',
                    'True',
                    '-templates:force_native_topology',
                    'True' ]
            cmd += tcmd


        return cmd
    ##End make_rosetta_cmd
    
    
    def doModelling(self):
        """
        Run the modelling and return the path to the models directory
        """

       

        # Make modelling directory and get stuff required for run
        os.mkdir(self.models_dir)
        seeds = self.generate_seeds()
        jobs = self.split_jobs()

        # List of processes so we can check when they are done
        processes = []
        # dict mapping process to directories
        directories = {}
        for proc in range(1,self.nproc+1):

            # Get directory to run job in 
            wdir = self.work_dir + os.sep + 'models_' + str(proc)
            directories[wdir] = proc
            os.mkdir(wdir)
            
            # Generate the command for this processor
            seed = str(seeds[proc])
            nstruct = str(jobs[proc])
            cmd = self.rosetta_cmd( wdir, nstruct, seed )
            
            self.logger.debug('Making {} models in directory: {}'.format(nstruct,wdir) )
            self.logger.debug('Executing cmd: {}'.format( " ".join(cmd)) )
            
            print 'Executing cmd: {}'.format( " ".join(cmd))
            
            logf = open(wdir+os.sep+"rosetta_{}.log".format(proc),"w")
            p = subprocess.Popen( cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=wdir )
            processes.append(p)

            #proc += 1
            
        #End spawning loop
        
        # Check to see if they have finished
        done=False
        completed=0
        retcodes = [None]*len(processes) # To hold return codes
        while not done:
            time.sleep(5)
            for i, p in enumerate(processes):
                ret = p.poll()
                retcodes[i] = ret
                if ret != None:
                    completed+=1
            
            if completed == len(processes):
                break
    
        # Check the return codes
        for i, ret in enumerate(retcodes):
            if ret != 0:
                msg = "Error generating models with Rosetta!Got error for processor: {}".format(i)
                logging.critical( msg )  
                raise RuntimeError, msg
        
        if self.use_scwrl:
            # Add sidechains using SCRWL - loop over each directory and output files into the models directory
            for wdir,proc in directories.iteritems():
                add_sidechains_SCWRL.add_sidechains_SCWRL(self.scwrl_exe, wdir, self.models_dir, str(proc), False)
        else:
        # Just copy all modelling files into models directory
            for wd in directories.keys():
                proc = directories[wd]
                for pfile in glob.glob( os.path.join(wd, '*.pdb') ):
                    pdbname = os.path.split(pfile)[1]
                    shutil.copyfile( wd + os.sep + pdbname, self.models_dir + os.sep + str(proc) + '_' + pdbname)

        return self.models_dir
    ##End doModelling
    
    def setup_domain_constraints(self):
        """
        Create the file for restricting the domain termini and return a list suitable
        for adding to the rosetta command
        """
        
        fas = open(self.fasta)
        seq = ''
        for line in fas:
            if not re.search('>', line):
                seq += line.rstrip('\n')
        length = 0
        for x in seq:
            if re.search('\w', x):
                length += 1
    
    
        self.logger.info('restricting termini distance: {}'.format( self.domain_termini_distance ))
        constraints_file = os.path.join(self.work_dir, 'constraints')
        conin = open(constraints_file, "w")
        conin.write('AtomPair CA 1 CA ' + str(length) + ' GAUSSIANFUNC ' + str(self.domain_termini_distance) + ' 5.0 TAG')
        cmd = '-constraints:cst_fa_file', constraints_file, '-constraints:cst_file', constraints_file
        
        return cmd

    def set_from_amopt(self, amopt ):
        """
        Set the values from the amopt object
        """
        
        self.nproc = amopt.d['nproc']
        self.nmodels = amopt.d['nmodels']
        self.work_dir = amopt.d['work_dir']
        
        # Set models directory
        if not amopt.d['models_dir']:
            self.models_dir = amopt.d['work_dir'] + os.sep + "models"
        else:
            self.models_dir = amopt.d['models_dir']
            
        self.rosetta_path = amopt.d['rosetta_path']
        self.rosetta_db = amopt.d['rosetta_db']
        self.frags3mers = amopt.d['frags3mers']
        self.frags9mers = amopt.d['frags9mers']
        self.fasta = amopt.d['fasta']
        self.use_scwrl = amopt.d['use_scwrl']
        self.scwrl_exe = amopt.d['scwrl_exe']
        
        # Extra options
        self.allatom = amopt.d['allatom']
        self.domain_termini_distance = amopt.d['domain_termini_distance']
        self.rad_gyr_reweight = amopt.d['CC']
        self.improve_template = amopt.d['improve_template']


class Test(unittest.TestCase):


    def setUp(self):
        """
        Create a dummy script for creating pdb files
        """
        
        script = "dummy_rosetta.sh"
        f = open(script,"w")
        content = """#!/usr/bin/env python
for i in range(10):
    f = open( "rosy_{}.pdb".format(i), "w")
    f.write( "rosy_{}.pdb".format(i) )
    f.close()"""
        f.write(content)
        os.chmod(script, 0o777)


    def tearDown(self):
        pass


    def testNoRosetta(self):
        """
        Test without Rosetta
        """
        amopt = ample_options.AmpleOptions()
        amopt.d['nproc'] = 3
        amopt.d['nmodels'] = 30
        amopt.d['work_dir'] = os.getcwd()
        amopt.d['models_dir'] = os.getcwd() + os.sep + "models"
        amopt.d['rosetta_path'] = os.getcwd() + os.sep + "dummy_rosetta.sh"
        amopt.d['rosetta_db'] = "ROSETTA_DB"
        amopt.d['frags3mers'] = '3mers'
        amopt.d['frags9mers'] = '9mers'
        amopt.d['fasta'] = "FASTA"
        amopt.d['allatom'] = True
        amopt.d['use_scwrl'] = False
        amopt.d['scwrl_exe'] = ""
        
        amopt.d['domain_termini_distance'] = None
        amopt.d['CC'] = None
        amopt.d['improve_template'] = None


        rm = RosettaModel()
        rm.set_from_amopt( amopt )
        mdir = rm.doModelling()
        print "models in: {}".format(mdir)



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
