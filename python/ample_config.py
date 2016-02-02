#!/usr/bin/env ccp4-python

'''
30.01.2016

@author: hlfsimko
'''

import logging
import os
import sys

try:
    import configparser as ConfigParser
except:
    import ConfigParser

# Need the root directory to find config file if no other is provided
if "CCP4_AMPLE_ROOT" in os.environ and "CCP4" in os.environ:
    root = os.environ["CCP4_AMPLE_ROOT"]
elif "CCP4" in os.environ:
    root = os.path.join(os.environ["CCP4"], "share", "ample")
else:
    raise RuntimeError("Cannot locate CCP4 install")

from ample_ensemble import SIDE_CHAIN_TREATMENTS
import version

_logger=logging.getLogger(__name__)

##############################################################
# The sections and options within need to be stored
# otherwise we cannot manage interplay between 
# ConfigParser and AMPLE settings dictionary.
# Only default non-dynamic part are the files as they
# should never be stored in the config file to avoid errors

_SECTIONS_REFERENCE = {"AMPLE_info" : ["ample_version"],
                       "Databases" : ['nr',
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
                                  'models',
                                  'mrbump_dir',
                                  'mr_sequence',
                                  'mtz',
                                  'native_pdb',
                                  'nmr_model_in',
                                  'nmr_remodel_fasta',
                                  'out_config_file',
                                  'psipred_ss2',
                                  'restart_pkl',
                                  'restraints_file',
                                  'score_matrix',
                                  'score_matrix_file_list',
                                  'sf_cif',
                                  'transmembrane_octopusfile',
                                  'transmembrane_lipofile',
                                  'transmembrane_spanfile',
                                  'work_dir'],
                        "Unspecified" : []
}


class AMPLEConfigOptions(object):
    
    def __init__(self):
        
        self.d = {} # store all options here
        self.debug = False
        
        self.quick_mode = {
                           'max_ensemble_models' : 10,
                           'nmodels' : 200,
                           'percent' : 20,
                           'shelx_cycles' : 5,
                           'use_arpwarp' : False,
                           'use_buccaneer' : False,
                           'phaser_kill' : 15
        }

        # Test use scrwl
        self.devel_mode = {
                           'early_terminate': False,
                           'benchmark_mode': True,
                           'shelxe_rebuild' : True,
                           'shelxe_rebuild_arpwarp' : True,
                           'shelxe_rebuild_buccaneer' : True,
                           'use_arpwarp' : False,
                           'use_buccaneer' : False,
                           #'mr_keys' : [ [ 'PKEY', 'KILL','TIME','360'  ] ],
        }
        
        self.webserver_uri = {
                               'purge': True,
                               'shelxe_rebuild_buccaneer': True,
                               'submit_cluster' : True,
                               'submit_max_array' : 10,
                               'submit_qtype' : "SGE",
                               'submit_queue' : "all.q",
        }
        
    def populate(self, config_opts):
        
        # Convert Namespace to Dictionary
        config_opts = vars(config_opts)

        # Identify which config file to use
        config_file = os.path.abspath(config_opts["config_file"]) \
            if config_opts["config_file"] else \
                os.path.join(root, "include", "ample.ini")
        _logger.debug("Using configuration file: {0}".format(config_file))

         # Read the configuration file
        self._read_config_file(config_file)
        # Read the command line arguments
        self._read_config_opts(config_opts) 

        # Set further options
        self._process_options()
        return
     
    def _process_options(self):
        
        self.d['ample_version'] = version.__version__
        
        if "rcdir" in self.d and not self.d["rcdir"]:
            self.d["rcdir"] = os.path.join(os.path.expanduser("~"), ".ample")
        
        if "run_dir" in self.d and not self.d["run_dir"]:
            self.d["run_dir"] = os.getcwd()
            
        if "side_chain_treatments" in self.d and not self.d["side_chain_treatments"]:
            self.d["side_chain_treatments"] = SIDE_CHAIN_TREATMENTS
        
        # Set full file paths
        for k, v in self.d.iteritems():
            if k in _SECTIONS_REFERENCE["Files"] and v:
                self.d[k] = os.path.abspath(v)
        
        # Any changes here
        if self.d['submit_qtype']:
            self.d['submit_qtype'] = self.d['submit_qtype'].upper()
        
        if self.d['shelxe_rebuild']:
            self.d['shelxe_rebuild_arpwap']=True
            self.d['shelxe_rebuild_buccaneer']=True
        
        # Check if using any preset options
        if self.d['devel_mode']: self._preset_options('devel_mode')
        if self.d['quick_mode']: self._preset_options('quick_mode')
        if self.d['webserver_uri']: self._preset_options('webserver_uri')
        
        # Write config to job specific directory
        self.d['out_config_file'] = os.path.join(self.d['work_dir'], self.d['name']+".ini")
        
        return
    
    def _preset_options(self ,mode):
        assert hasattr(self, mode),"Unknown mode: {0}".format(mode)
        for k, v in getattr(self, mode).iteritems():
            # Set any that haven't been set
            if self.d[k] == None:
                self.d[k] = v
            else:
                # Already set - only overwrite if it's set to a default value, otherwise we
                # let the user go with what they've chosen but warn
                if self.d[k] != self.defaults[k] and self.d[k] != v  :
                    print "WARNING! Overriding {0} setting: {1} : {2} with user setting {3}".format(mode, k, v, self.d[k])
                else:
                    # We overwrite the default with our value
                    if self.debug:
                        print "Overriding default setting: {0} : {1} with {2} setting {3}".format(k, mode, self.defaults[k], v)
                    self.d[k] = v
        return
        
    def _read_config_file(self, config_file):
        config = ConfigParser.SafeConfigParser()
        config.read(config_file)
        
        for section in config.sections():
            
            if not section in _SECTIONS_REFERENCE:
                _SECTIONS_REFERENCE[section] = []
            
            # Basic switch statement to determine the type of the variable
            for k, v in config.items(section):
                if v.lower() == "none":
                    self.d[k] = None
                    
                elif v.lower() == "true":
                    self.d[k] = True
                    
                elif v.lower() == "false":
                    self.d[k] = False
                
                elif section.lower() == "databases":
                    self.d[k] = os.path.abspath(v)    
                
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
        tmpv = None
        cmdline_flags = []
        
        for k, v in config_opts.iteritems():
            if v is not None: cmdline_flags.append(k)
            
            tmpv = v[0] if isinstance(v, list) else v
                        
            if isinstance(tmpv, str):
                if tmpv.lower() == "true":
                    tmpv = True
                elif tmpv.lower() == "false":
                    tmpv = False
                elif tmpv.lower() == "none":
                    tmpv = None

            if k not in self.d:
                self.d[k] = tmpv
                _SECTIONS_REFERENCE["Unspecified"].append(k)
            elif tmpv != None: 
                _logger.debug("Changing {0}: {1} => {2}".format(k, 
                                                                self.d[k], 
                                                                tmpv))
                self.d[k] = tmpv
            
            
        self.d['cmdline_flags'] = cmdline_flags
        return
    
    def _isfloat(self, value):
        try: 
            float(value)
            return True
        except:
            return False
            
    def prettify_parameters(self):
        """Return the parameters nicely formated as a list of strings suitable 
        for writing out to a file"""
        pstr ='Parameters Used in this Run\n\n'
        for k, v in sorted(self.d.items()):
            pstr += "{0} : {1}\n".format(k, v)
        return pstr
    
    def write_config_file(self):
        config = ConfigParser.SafeConfigParser()
        self._write_config_file(config)
        _logger.info("AMPLE configuration written to: {0}".format(self.d['out_config_file']))
        with open(self.d['out_config_file'], "w") as out: 
            config.write(out)
        return
    
    def _write_config_file(self, config_parser):
        # Add all sections to the configparser
        for section in sorted(_SECTIONS_REFERENCE.keys()):
            config_parser.add_section(section)
            
        # Place all entries in our dictionary in the corresponding section in
        # the configparser
        for option, value in self.d.iteritems():
            # Extract the section in which the entry needs to go
            sections = [k for (k, v) in _SECTIONS_REFERENCE.items() if option in v]
            
            # Make sure we only have each option assigned to a single section
            section = None if len(sections) > 1 else sections[0]    
            assert section, "Uncertainty about option: {0}".format(option)
            
            # We do not want to re-use files or at least not by default.
            # Comment those specifically out to avoid any errors
            if section.lower() == "ample_info":
                config_parser.set(section, "#"+option, str(value))
            elif section.lower() == "files":
                config_parser.set(section, "#"+option, str(value))
            else:
                config_parser.set(section, option, str(value))

        return

if __name__ == "__main__":
    import argparse
    options = argparse.ArgumentParser()
    options.add_argument("config_file")
    optd = vars(options.parse_args())
    
    cp = AMPLEConfigOptions(optd)
    cp.write_config_file("check.config")
