#!/usr/bin/env ccp4-python

'''
22.03.2015

@author: hlfsimko
'''

import logging
import os
import sys

# Renamed in python 3
try:
    import configparser as ConfigParser
except:
    import ConfigParser

if "CCP4_AMPLE_ROOT" in os.environ and "CCP4" in os.environ:
    root = os.environ["CCP4_AMPLE_ROOT"]
elif "CCP4" in os.environ:
    root = os.path.join(os.environ["CCP4"], "share", "ample")
else:
    raise RuntimeError("Cannot locate CCP4 install")

from ample_ensemble import SIDE_CHAIN_TREATMENTS
import version

##############################################################
# The sections and options within need to be stored
# otherwise we cannot manage interplay between 
# ConfigParser and AMPLE settings dictionary.
# Only default non-dynamic part are the files as they
# should never be stored in the config file to avoid errors

_SECTIONS_REFERENCE = {"Databases" : ['nr',
                                      'rosetta_db'],
                       "Executables" : ['blast_dir',
                                        'fast_protein_cluster_exe',
                                        'gesamt_exe',
                                        'maxcluster_exe',
                                        'mustang_exe',
                                        'rosetta_dir',
                                        'rosetta_fragments_exe',
                                        'rosetta_AbinitioRelax',
                                        'scwrl_exe',
                                        'shelxe_exe',
                                        'spicker_exe',
                                        'theseus_exe'],
                       
                       "Files" : ['alignment_file',
                                  'bbcontacts_file',
                                  'cluster_dir',
                                  'config_file',
                                  'contact_file',
                                  'disulfide_constraints_file',
                                  'domain_all_chains_pdb',
                                  'ensembles',
                                  'fasta',
                                  'frags_3mers',
                                  'frags_9mers',
                                  'import_cluster',
                                  'models',
                                  'mrbump_dir',
                                  'mr_sequence',
                                  'mtz',
                                  'native_pdb',
                                  'nmr_model_in',
                                  'nmr_remodel_fasta',
                                  'psipred_ss2',
                                  'restart_pkl',
                                  'restraints_file',
                                  'rosetta_db',
                                  'score_matrix',
                                  'score_matrix_file_list',
                                  'sf_cif',
                                  'transmembrane_octopusfile',
                                  'transmembrane_lipofile',
                                  'transmembrane_spanfile',
                                  'work_dir']
}


class AMPLEConfigOptions(object):
    
    def __init__(self, config_opts=None):
        # Define some options
        self.d = {}
        self.logger = logging.getLogger()
    
        # Identify which config file to use
        config_file = os.path.abspath(config_opts["config_file"]) \
            if config_opts["config_file"] else \
                os.path.join(root, "include", "ample.config")
        
        self._read_config_file(config_file) # Read the configuration file
        self._read_config_opts(config_opts) # Read the command line arguments
        
        self.d['ample_version'] = version.__version__
        
        if not self.d["rcdir"]:
            self.d["rcdir"] = os.path.join(os.path.expanduser("~"), ".ample")
        
        if not self.d["run_dir"]:
            self.d["run_dir"] = os.getcwd()
            
        if not self.d["side_chain_treatments"]:
            self.d["side_chain_treatments"] = SIDE_CHAIN_TREATMENTS
        
        # Set further options
        self._process_options()
        return
        
    def prettify_parameters(self):
        """Return the parameters nicely formated as a list of strings suitable 
        for writing out to a file"""
        pstr ='Parameters Used in this Run\n\n'
        for k, v in sorted(self.d.items()):
            pstr += "{0} : {1}\n".format(k, v)
        return pstr
    
    def write_config_file(self, filename):
        config = ConfigParser.SafeConfigParser()
    
        # Add all sections to the configparser
        for section in sorted(_SECTIONS_REFERENCE.keys()):
            config.add_section(section)
            
        # Place all entries in our dictionary in the corresponding section in
        # the configparser
        for option, value in self.d.iteritems():
            value = str(value)
            
            # Extract the section in which the entry needs to go
            for k, v in _SECTIONS_REFERENCE.iteritems():
                if option in v:
                    section = k
            
            # We do not want to re-use files or at least not by default.
            # Comment those specifically out to avoid any errors
            if section.lower() == "files":
                config.set(section, "#"+option, value)
            else:
                config.set(section, option, value)
        
        with open(filename, "w") as out: config.write(out)
        return
     
    def _process_options(self):
        self._full_file_paths()
        self._exec_extensions()
        return
        
    def _exec_extensions(self):
        return
        
    def _full_file_paths(self):
        for f in _SECTIONS_REFERENCE["Files"]:
            if f in self.d and self.d[f]:
                file = os.path.abspath(self.d[f])
                if os.path.isfile(file) and os.path.getsize(file) > 0:
                    self.d[f] = file
                elif os.path.isfile(file):
                    raise RuntimeError('%s exists but is empty' % file)
                else:
                    raise RuntimeError('%s does not exist' % file)
            else:
                continue
        return
        
    def _read_config_file(self, config_file):
        config = ConfigParser.SafeConfigParser()
        config.read(config_file)
        
        for section in config.sections():
            
            if not section in _SECTIONS_REFERENCE:
                _SECTIONS_REFERENCE[section] = []
            
            # Basic switch statement to determine the type of the variable
            for k, v in config.items(section):
                if v == "None":
                    self.d[k] = None
                    
                elif v == "True":
                    self.d[k] = True
                    
                elif v == "False":
                    self.d[k] = False
                    
                elif section.lower() == "executables": 
                    self.d[k] = os.path.abspath(v)
                    
                elif section.lower() == "files":
                    self.d[k] = os.path.abspath(v)
                    
                elif v.isdigit():
                    self.d[k] = int(v)
                    
                elif self._isfloat(v):
                    self.d[k] = float(v)
                    
                else: 
                    self.d[k] = v
                    
                _SECTIONS_REFERENCE[section].append(k)
                
        return
    
    def _read_config_opts(self, config_opts):
        for k, v in config_opts.items():            
            self.d[k]=v if config_opts[k] else None
        return
    
    def _isfloat(self, value):
        try: 
            float(value)
            return True
        except:
            return False


if __name__ == "__main__":
    import argparse
    options = argparse.ArgumentParser()
    options.add_argument("config_file")
    optd = vars(options.parse_args())
    
    cp = AMPLEConfigOptions(optd)
    cp.write_config_file("check.config")
