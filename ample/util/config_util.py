#!/usr/bin/env ccp4-python

'''
30.01.2016

@author: hlfsimko
'''

import logging
import os

from ample.constants import SHARE_DIR
from ample.ensembler.constants import SIDE_CHAIN_TREATMENTS
from ample.util import version

# Python 3.x --> ConfigParser renamed to configparser
try:
    import configparser as ConfigParser
except ImportError:
    import ConfigParser 

LOGGER = logging.getLogger(__name__)

##############################################################
# The sections and options within need to be stored
# otherwise we cannot manage interplay between 
# ConfigParser and AMPLE settings dictionary.
# Some default non-dynamic parts are stored below to avoid errors

_SECTIONS_REFERENCE = {"AMPLE_info" : ["ample_version",
                                       "ccp4_version",
                                       "cmdline_flags"],
                       
                       "Databases" : ['nr',
                                      'rosetta_db'],
                       
                       "Executables" : ['blast_dir',
                                        'cluster_exe',
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
                                  'ample_log',
                                  'bbcontacts_file',
                                  'cluster_dir',
                                  'config_file',
                                  'contact_file',
                                  'disulfide_constraints_file',
                                  'domain_all_chains_pdb',
                                  'ensembles',
                                  'ensembles_directory',
                                  'ensemble_ok',
                                  'fasta',
                                  'frags_3mers',
                                  'frags_9mers',
                                  'models',
                                  'models_dir',
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
                                  'results_path',
                                  'score_matrix',
                                  'score_matrix_file_list',
                                  'sf_cif',
                                  'single_model',
                                  'transmembrane_octopusfile',
                                  'transmembrane_lipofile',
                                  'transmembrane_spanfile',
                                  'truncation_scorefile',
                                  'work_dir'],
                        # Data stored in amopt.d but not really part of AMPLE's configuration
                        "No_config" : ["benchmark_results",
                                       "ensembles_data",
                                       "fasta_length",
                                       "mrbump_results",
                                       "sequence",
                                       "truncation_variances",
                                       "truncation_levels",
                                       "truncation_nresidues"],
                       # In case we haven't specified anything or it is new
                        "Unspecified" : [],
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
                           'benchmark_mode': True,
                           'early_terminate': False,
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
        self.config_opts = config_opts = vars(config_opts)

        # Identify which config file to use
        config_file = self._get_config_file(config_opts['config_file'])

        # Read the configuration file
        self._read_config_file(config_file)
        # Read the command line arguments
        self._read_config_opts(config_opts) 

        # Set further options
        self._process_options()
        return
     
    def _get_config_file(self, cmd_file=None):
        config_file = os.path.abspath(cmd_file) if cmd_file else \
                            os.path.join(SHARE_DIR, "include", "ample.ini")
        if not os.path.isfile(config_file):
            msg = "Cannot find configuration file - terminating..."
            LOGGER.critical(msg)
            raise RuntimeError(msg)
        LOGGER.debug("Using configuration file: {0}".format(config_file))
        return config_file
     
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
        
        return
    
    def _preset_options(self, mode):
        assert hasattr(self, mode),"Unknown mode: {0}".format(mode)
        for k, v in getattr(self, mode).iteritems():
            # Set any that haven't been set
            if not k in self.d or self.d[k] == None:
                self.d[k] = v
            else:
                # Already set - only overwrite if it's set to a default value, otherwise we
                # let the user go with what they've chosen but warn
                if k in self.config_opts.keys() and self.d[k] != self.config_opts[k] and self.d[k] != v:
                    print "WARNING! Overriding {0} setting: {1} : {2} with {3}".format(mode, k, self.d[k], v)
                    self.d[k] = v
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
            elif tmpv != None: 
                print("Changing {0}: {1} => {2}".format(k, self.d[k], tmpv))
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
        # Write config to job specific directory
        self.d["out_config_file"] = f = os.path.join(self.d['work_dir'], 
                                                     self.d['name']+".ini")
        LOGGER.info("AMPLE configuration written to: {0}".format(f))
        with open(f, "w") as out: config.write(out)
        return
    
    def _write_config_file(self, config_parser):
        # Add all sections to the configparser
        for section in sorted(_SECTIONS_REFERENCE.keys()):
            if section.lower() == "no_config": continue
            config_parser.add_section(section)
        
        # Place all entries in our dictionary in the corresponding section in
        # the configparser
        for option in sorted(self.d.keys()):
            # Extract the section in which the entry needs to go
            sections = [k for (k, v) in _SECTIONS_REFERENCE.items() \
                            if any(entry.lower() == option.lower() for entry in v)]

            # Make sure we only have each option assigned to a single section
            section = "Unspecified" if len(sections) != 1 else sections[0]

            # We do not want to re-use files or at least not by default.
            # Comment those specifically out to avoid any errors
            if section.lower() == "no_config":
                continue
            elif section.lower() == "ample_info" or \
                 section.lower() == "files" or \
                 section.lower() == "unspecified":
                config_parser.set(section, "#" + option, str(self.d[option]))
            else:
                config_parser.set(section, option, str(self.d[option]))

        return

