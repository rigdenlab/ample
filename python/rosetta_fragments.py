'''
Created on 19 Feb 2013

Class to create the fragments for ROSETTA

@author: jmht
'''

# Python imports
import logging
import os
import shutil
import subprocess
import unittest

# Our imports
import ample_options
import ample_util

class FragmentGenerator(object):
    """
    Generate the 3 & 9 fragments required by Rosetta
    """
    
    def __init__(self):
        
        self.make_fragments_exe = None
        self.fasta = None
        self.pdb_code = None
        self.fragments_directory = None
        self.frags3mers = None
        self.frags9mers = None
        self.extraflags = []
        
    def set_from_amopt(self, amopt ):
        """
        Set the values from the amopt object
        """
        
        if not "rosetta_dir" in amopt.d.keys() and "rosetta_version" in amopt.d.keys():
            raise RuntimeError,"FragmentGenerator: need rosetta_dir and rosetta_version to be set"
        
        # Set path to script
        if not amopt.d['rosetta_fragments_exe']: 
            if amopt.d['rosetta_version'] == '3.3':
                self.make_fragments_exe = amopt.d['rosetta_dir'] + '/rosetta_fragments/make_fragments.pl'
            elif amopt.d['rosetta_version']  == '3.4':
                self.make_fragments_exe = amopt.d['rosetta_dir'] + '/rosetta_tools/fragment_tools/make_fragments.pl'
        else:
            self.make_fragments_exe = amopt.d['rosetta_fragments_exe']
        
        # version dependent flag
        if amopt.d['rosetta_version'] == '3.3':
            # jmht the last 3 don't seem to work with 3.4
            self.extraflags += ['-noporter', '-nojufo', '-nosam','-noprof' ]
        elif amopt.d['rosetta_version'] == '3.4':
            self.extraflags += ['-old_name_format']
            
        # Whether to exclude homologs
        if not amopt.d['usehoms']:
            self.extraflags.append('-nohoms')
            
        self.fasta = amopt.d['fasta']
        self.pdb_code = amopt.d['pdb_code']
        self.fragments_directory = amopt.d['work_dir'] + os.sep + "rosetta_fragments"
        
        # Run our checks
        if not os.path.exists( self.make_fragments_exe ):
            raise RuntimeError, "Can't find make_fragments script: {}\nCheck the path or set the rosetta_fragments_exe flag.".format( self.make_fragments_exe )
    
    
    def generate_fragments(self):
        """
        Run the script to generate the fragments
        
        """
    
        logger = logging.getLogger()
        logger.info('----- making fragments--------')
        #RUNNING.write('----- making fragments--------\n')
        #RUNNING.flush()
        
        if not os.path.exists( self.fragments_directory ):
            os.mkdir( self.fragments_directory )
            
        # It seems that the script can't tolerate "-" in the directory name leading to the fasta file,
        # so we need to copy the fasta file into the fragments directory
        fasta = os.path.split(  self.fasta )[1]
        shutil.copy2( self.fasta, self.fragments_directory + os.sep + fasta )
            
        cmd = [ self.make_fragments_exe, '-rundir', self.fragments_directory, '-id', self.pdb_code ] + self.extraflags + [ fasta ] 
        logger.info('Executing cmd: {}'.format( " ".join(cmd)) )

        os.chdir( self.fragments_directory )
        try:
            #output = subprocess.check_output( cmd, shell=True, stderr=subprocess.STDOUT )
            #cmd = " ".join(cmd)
            output = subprocess.check_output( cmd, stderr=subprocess.STDOUT )
            f = open("make_fragments.log","w")
            f.writelines( output )
            f.close()
        except subprocess.CalledProcessError,e:
            logger.critical("Error generating fragments:\n{}".format(e.output))
            #raise RuntimeError,e
            sys.exit(1)
            
        # Move back out of the directory when we're finished
        os.chdir( ".." )
            
        # old_name_format
        self.frags3mers = self.fragments_directory + os.sep + 'aa' + self.pdb_code + '03_05.200_v1_3'
        self.frags9mers = self.fragments_directory + os.sep + 'aa' + self.pdb_code + '09_05.200_v1_3'
    
        if not os.path.exists( self.frags3mers ) or not os.path.exists( self.frags9mers ):
            raise RuntimeError, "Error in making fragments"
        
    
        #RUNNING.write('Fragments done\n3mers at: ' + frags_3_mers + '\n9mers at: ' + frags_9_mers + '\n\n')
        print' Fragments Done\n3mers at: ' + self.frags3mers + '\n9mers at: ' + self.frags9mers + '\n\n'
    
        if os.path.exists( self.fragments_directory + os.sep + self.fragments_directory + '.psipred'):
            ample_util.get_psipred_prediction( self.fragments_directory + os.sep + self.pdb_code + '.psipred')
        

class TestUtil(unittest.TestCase):
    """
    Unit test
    """
    
    def setUp(self):
        """
        Get paths need to think of a sensible way to do this
        """
        thisdir = os.getcwd()
        self.ampledir = os.path.abspath( thisdir+os.sep+"..")
    
    def testFragmentGenerator(self):
        """See we can create fragments"""
        
        print "testing FragmentGenerator"
        
        amopt = ample_options.AmpleOptions()
        amopt.d['rosetta_dir'] = "/opt/rosetta3.4"
        amopt.d['rosetta_version'] = "3.4"
        amopt.d['pdb_code'] = "TOXD_"
        amopt.d['work_dir'] =  os.getcwd()
        amopt.d['usehoms'] =  True
        amopt.d['rosetta_fragments_exe'] =  "/tmp/make_fragments.pl"
        amopt.d['rosetta_fragments_exe'] =  None
        amopt.d['fasta'] = self.ampledir + "/examples/toxd-example/toxd_.fasta"
        
        fg = FragmentGenerator()
        fg.set_from_amopt( amopt )
        fg.generate_fragments()

# Run unit tests
if __name__ == "__main__":
    unittest.main()
