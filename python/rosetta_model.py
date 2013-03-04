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
import sys
import time
import unittest

# Our modules
import add_sidechains_SCWRL
import ample_options
import ample_util


class RosettaModel(object):
    """
    Class to run Rosetta modelling
    """

    def __init__(self):
        
        self.nproc = None
        self.nmodels = None
        self.work_dir = None
        self.models_dir = None
        self.rosetta_dir = None
        self.rosetta_path = None
        self.rosetta_cm = None
        self.rosetta_db = None
        self.rosetta_version = None
        
        # Not used yet
        self.make_models = None
        self.make_fragments = None

        self.fasta = None
        self.allatom = None
        self.use_scwrl = None
        self.scwrl_exe = None
        
        # Fragment variables
        self.pdb_code = None
        self.frags3mers = None
        self.frags9mers = None
        self.usehoms = None
        self.fragments_directory = None
        self.fragments_exe = None
        
        # List of seeds
        self.seeds = None
        
        # Extra options
        self.domain_termini_distance = None
        self.rad_gyr_reweight = None
        self.improve_template = None
        
        self.logger = logging.getLogger()
    
    def generate_seeds(self, nseeds):
        """
        Generate a list of nseed seeds
        """
        
        seed_list = []

        # Generate the list of random seeds
        while len(seed_list) < nseeds:
            seed = random.randint(1000000, 4000000)
            if seed not in seed_list:
                seed_list.append(seed)
                
        # Keep a log of the seeds
        seedlog = open( self.work_dir + os.sep +'seedlist', "w")
        for seed in  seed_list:
            seedlog.write(str(seed) + '\n')
        seedlog.close()
        
        self.seeds = seed_list
        return
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
    
    def fragment_cmd(self):
        """
        Return the command to make the fragments as a list
        """
        # Set path to script
        if not self.fragments_exe: 
            if self.rosetta_version == '3.3':
                self.fragments_exe = self.rosetta_dir + '/rosetta_fragments/make_fragments.pl'
            elif self.rosetta_version  == '3.4':
                self.fragments_exe = self.rosetta_dir + '/rosetta_tools/fragment_tools/make_fragments.pl'
        
        # It seems that the script can't tolerate "-" in the directory name leading to the fasta file,
        # so we need to copy the fasta file into the fragments directory and just use the name here
        fasta = os.path.split(  self.fasta )[1]
         
        cmd = [ self.fragments_exe,
               '-rundir', self.fragments_directory,
               '-id', self.pdb_code, fasta ] 
        
        # version dependent flags
        if self.rosetta_version == '3.3':
            # jmht the last 3 don't seem to work with 3.4
            cmd += ['-noporter', '-nojufo', '-nosam','-noprof' ]
        elif self.rosetta_version == '3.4':
            cmd += ['-old_name_format']
            
        # Whether to exclude homologs
        if self.usehoms:
            cmd.append('-nohoms')
        
        return cmd
    ##End fragment_cmd
            

    def generate_fragments(self):
        """
        Run the script to generate the fragments
        
        """
    
        self.logger.info('----- making fragments--------')
        #RUNNING.write('----- making fragments--------\n')
        #RUNNING.flush()
        
        if not os.path.exists( self.fragments_directory ):
            os.mkdir( self.fragments_directory )
            
        # It seems that the script can't tolerate "-" in the directory name leading to the fasta file,
        # so we need to copy the fasta file into the fragments directory
        fasta = os.path.split(  self.fasta )[1]
        shutil.copy2( self.fasta, self.fragments_directory + os.sep + fasta )
            
        cmd = self.fragment_cmd()
        self.logger.info('Executing cmd: {}'.format( " ".join(cmd)) )

        try:
            output = subprocess.check_output( cmd, stderr=subprocess.STDOUT, cwd=self.fragments_directory )
            f = open( self.fragments_directory + os.sep + "make_fragments.log","w")
            f.writelines( output )
            f.close()
        except subprocess.CalledProcessError,e:
            self.logger.critical("Error generating fragments:\n{}".format(e.output))
            raise RuntimeError,e
            
        # old_name_format
        self.frags3mers = self.fragments_directory + os.sep + 'aa' + self.pdb_code + '03_05.200_v1_3'
        self.frags9mers = self.fragments_directory + os.sep + 'aa' + self.pdb_code + '09_05.200_v1_3'
    
        if not os.path.exists( self.frags3mers ) or not os.path.exists( self.frags9mers ):
            raise RuntimeError, "Error making fragments - could not find fragment files:\n{}\n{}\n".format(self.frags3mers,self.frags9mers)
        
    
        #RUNNING.write('Fragments done\n3mers at: ' + frags_3_mers + '\n9mers at: ' + frags_9_mers + '\n\n')
        self.logger.info('Fragments Done\n3mers at: ' + self.frags3mers + '\n9mers at: ' + self.frags9mers + '\n\n')
    
        if os.path.exists( self.fragments_directory + os.sep + self.fragments_directory + '.psipred'):
            ample_util.get_psipred_prediction( self.fragments_directory + os.sep + self.pdb_code + '.psipred')
       
    ##End fragment_cmd
    
    def get_version(self):
        """ Return the Rosetta version as a string"""
        
        if not self.rosetta_version:
            # Get version
            version_file = self.rosetta_dir + '/README.version'
            if not os.path.exists(version_file):
                self.logger.critical('version file for Rosetta not found')
                sys.exit()
                
            version = '3.2'
            try:
                for line in open(version_file,'r'):
                    line.strip()
                    if line.startswith('Rosetta'):
                        version = line.split()[1].strip()
                self.logger.info( 'Your Rosetta version is: {}'.format( version ) )
            except Exception,e:
                print e
                self.logger.critical("Error determining rosetta version")
                sys.exit(1)
            self.rosetta_version = version
            
        return version
    #End get_version 
    
    def modelling_cmd(self, wdir, nstruct, seed):
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
        # Now generate the seeds
        self.generate_seeds( self.nproc )
        jobs = self.split_jobs()

        # List of processes so we can check when they are done
        processes = []
        # dict mapping process to directories
        directories = {}
        for proc in range(1,self.nproc+1):
            
            print "proc ",proc

            # Get directory to run job in 
            wdir = self.work_dir + os.sep + 'models_' + str(proc)
            directories[wdir] = proc
            os.mkdir(wdir)
            
            # Generate the command for this processor
            seed = str(self.seeds[proc-1])
            nstruct = str(jobs[proc])
            cmd = self.modelling_cmd( wdir, nstruct, seed )
            
            self.logger.debug('Making {} models in directory: {}'.format(nstruct,wdir) )
            self.logger.debug('Executing cmd: {}'.format( " ".join(cmd) ) )
            
            logf = open(wdir+os.sep+"rosetta_{}.log".format(proc),"w")
            p = subprocess.Popen( cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=wdir )
            processes.append(p)
            
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


        if not amopt.d['rosetta_dir'] or not os.path.isdir( amopt.d['rosetta_dir'] ):
            msg = "You need to set the rosetta_dir variable to where rosetta is installed"
            self.logger.critical(msg)
            raise RuntimeError,msg
        
        self.rosetta_dir = amopt.d['rosetta_dir']
    
        # Determin version
        amopt.d['rosetta_version'] = self.get_version()

        # Common variables
        self.fasta = amopt.d['fasta']
        self.work_dir = amopt.d['work_dir']
        self.pdb_code = amopt.d['pdb_code']
        
        # Fragment variables
        self.fragments_exe = amopt.d['rosetta_fragments_exe']
        self.frags3mers = amopt.d['frags3mers']
        self.frags9mers = amopt.d['frags9mers']
        self.usehoms = amopt.d['usehoms']
        self.fragments_directory = amopt.d['work_dir'] + os.sep + "rosetta_fragments"


        # Modelling variabls
        if amopt.d['make_models']:
            
            if not amopt.d['make_frags']:
                if not os.path.exists(self.frags3mers) or not os.path.exists(self.frags9mers):
                    msg = "Cannot find both fragment files:\n{}\n{}\n".format(self.frags3mers,self.frags9mers)
                    self.logger.critical(msg)
                    raise RuntimeError,msg
                    
            import platform
            if not amopt.d['rosetta_path']: 
                if platform.mac_ver() == ('', ('', '', ''), ''):
                    amopt.d['rosetta_path'] = self.rosetta_dir + '/rosetta_source/bin/AbinitioRelax.linuxgccrelease'
                else:
                    amopt.d['rosetta_path'] = self.rosetta_dir + '/rosetta_source/bin/AbinitioRelax'
            
            # Always save everything back to the amopt object so we can print it out
            self.rosetta_path = amopt.d['rosetta_path']
            
            if not os.path.exists(self.rosetta_path):
                logger.critical(' cant find Rosetta abinitio: {}'.format(self.rosetta_path) )
                raise RuntimeError,msg
    
            #jmht not used
            # ROSETTA_cluster = rosetta_dir + '/rosetta_source/bin/cluster.linuxgccrelease'
            if not amopt.d['rosetta_db']:
                amopt.d['rosetta_db'] = self.rosetta_dir + '/rosetta_database' 
            self.rosetta_db = amopt.d['rosetta_db']
            
            if not os.path.exists(self.rosetta_db):
                msg = ' cant find Rosetta DB: {}'.format(self.rosetta_db)
                self.logger.critical( msg )
                raise RuntimeError,msg
            
            #if not amopt.d['rosetta_cm']:
            #    amopt.d['rosetta_cm'] = self.rosetta_dir + '/rosetta_source/bin/idealize_jd2.default.linuxgccrelease'
            #self.rosetta_cm = amopt.d['rosetta_cm']
            #if not os.path.exists(ROSETTA_cluster):
            #    logger.critical(' cant find Rosetta cluster, check path names')
            #    sys.exit()
            
            self.nproc = amopt.d['nproc']
            self.nmodels = amopt.d['nmodels']
            # Set models directory
            if not amopt.d['models_dir']:
                self.models_dir = amopt.d['work_dir'] + os.sep + "models"
            else:
                self.models_dir = amopt.d['models_dir']
            
            # Extra modelling options
            self.allatom = amopt.d['allatom']
            self.domain_termini_distance = amopt.d['domain_termini_distance']
            self.rad_gyr_reweight = amopt.d['CC']
            
            if amopt.d['improve_template'] and not os.path.exists( amopt.d['improve_template'] ):
                msg = 'cant find template to improve'
                logger.critical( msg)
                raise RuntimeError(msg)
            self.improve_template = amopt.d['improve_template']
                
            self.use_scwrl = amopt.d['use_scwrl']
            self.scwrl_exe = amopt.d['scwrl_exe']        


class Test(unittest.TestCase):


    def setUp(self):
        """
        Set paths
        """
        logging.basicConfig()
        logging.getLogger().setLevel(logging.DEBUG)
        thisdir = os.getcwd()
        self.ampledir = os.path.abspath( thisdir+os.sep+"..")


    def testMakeFragments(self):
        """See we can create fragments"""
        
        print "testing FragmentGenerator"
        
        amopt = ample_options.AmpleOptions()
        amopt.d['rosetta_dir'] = "/opt/rosetta3.4"
        amopt.d['pdb_code'] = "TOXD_"
        amopt.d['work_dir'] =  os.getcwd()
        amopt.d['usehoms'] =  True
        amopt.d['make_frags'] = True
        amopt.d['rosetta_db'] = None
        amopt.d['rosetta_fragments_exe'] =  "/tmp/make_fragments.pl"
        #amopt.d['rosetta_fragments_exe'] =  None
        amopt.d['fasta'] = self.ampledir + "/examples/toxd-example/toxd_.fasta"
        
        amopt.d['make_models'] = False
        amopt.d['frags3mers'] = None
        amopt.d['frags9mers'] = None
        amopt.d['improve_template'] = None
        
        m = RosettaModel()
        m.set_from_amopt( amopt )
        m.generate_fragments()


    def testNoRosetta(self):
        """
        Test without Rosetta
        """
        
        ## Create a dummy script
        script = "dummy_rosetta.sh"
        f = open(script,"w")
        content = """#!/usr/bin/env python
for i in range(10):
    f = open( "rosy_{}.pdb".format(i), "w")
    f.write( "rosy_{}.pdb".format(i) )
    f.close()"""
        f.write(content)
        f.close()
        os.chmod(script, 0o777)
        
          
        # Set options
        amopt = ample_options.AmpleOptions()
        amopt.d['nproc'] = 3
        amopt.d['nmodels'] = 30
        amopt.d['work_dir'] = os.getcwd()
        amopt.d['models_dir'] = os.getcwd() + os.sep + "models"
        amopt.d['rosetta_dir'] = "/opt/rosetta3.4"
        amopt.d['rosetta_path'] = os.getcwd() + os.sep + "dummy_rosetta.sh"
        amopt.d['rosetta_db'] = None
        amopt.d['frags3mers'] = '3mers'
        amopt.d['frags9mers'] = '9mers'
        amopt.d['rosetta_fragments_exe'] = None
        amopt.d['usehoms'] = None
        amopt.d['make_models'] = True
        amopt.d['make_frags'] =  True
        amopt.d['fasta'] = "FASTA"
        amopt.d['pdb_code'] = "TOXD_"
        amopt.d['improve_template'] = None
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
