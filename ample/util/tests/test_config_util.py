"""Tests for util.config_util"""

import os
import tempfile
import unittest

from ample import constants
from ample.util import config_util
from ample.util import version

__author__ = "Adam Simpkin & Felix Simkovic"

class Test(unittest.TestCase):
    MAX_DIFF = None
    
    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share,'testfiles')
    
    def test_process_options(self):
        #Test the process_options
        options = config_util.AMPLEConfigOptions()
        options.d = {'fasta' : 'foo',
                     'side_chain_treatments' : [],
                     'rcdir' : os.path.join(os.sep, 'foo', 'bar'),
                     'nproc' : 1,
                     'submit_qtype' : False,
                     'shelxe_rebuild' : True,
                     'shelxe_rebuild_arpwarp' : False,
                     'shelxe_rebuild_buccaneer' : False,
                     'classic_mode' : False,
                     'devel_mode' : False,
                     'quick_mode' : False,
                     'thin_clusters' : False,
                     'webserver_uri' : False,
        }
        expected = {'ample_version' : version.__version__,
                    'fasta' : os.path.join(os.getcwd(), 'foo'),
                    'rcdir' : os.path.join(os.sep, 'foo', 'bar'),
                    'side_chain_treatments' : ['polyAla', 'reliable', 'allatom'],
                     'nproc' : 1,
                    'submit_qtype' : False,
                    'shelxe_rebuild' : True,
                    'shelxe_rebuild_arpwarp' : True,
                    'shelxe_rebuild_buccaneer' : True,
                    'classic_mode' : False,
                    'devel_mode' : False,
                    'quick_mode' : False,
                    'thin_clusters' : False,
                    'webserver_uri' : False,
        }
        options._process_options()
        self.assertItemsEqual(options.d, expected)
        
    def test_preset_options_quick_mode(self):
        options = config_util.AMPLEConfigOptions()
        options.d = {'ensemble_max_models' : 1000,
                     'nmodels' : 1000,
                     'percent' : 5,
                     'shelx_cycles' : 15,
                     'use_arpwarp' : True,
                     'use_buccaneer' : True,
                     'phaser_kill' : 360,
        }
        options.cmdline_opts = {}
        expected = {'ensemble_max_models' : 10,
                    'nmodels' : 200,
                    'percent' : 20,
                    'shelx_cycles' : 5,
                    'use_arpwarp' : False,
                    'use_buccaneer' : False,
                    'phaser_kill' : 15,
        }
        options._preset_options("quick_mode")
        self.assertEqual(options.d, expected)
    
    def test_preset_options_devel_mode(self):
        options = config_util.AMPLEConfigOptions()
        options.d = {'just_to_check' : True,
                     'early_terminate': True,
                     'benchmark_mode': False,
                     'shelxe_rebuild' : None,
                     'shelxe_rebuild_arpwarp' : False,
                     'shelxe_rebuild_buccaneer' : False,
                     'use_arpwarp' : None,
                     'use_buccaneer' : False,
        }
        options.cmdline_opts = {}
        expected = {'just_to_check' : True,
                    'benchmark_mode': True,
                    'early_terminate': False,
                    'shelxe_rebuild' : True,
                    'shelxe_rebuild_arpwarp' : True,
                    'shelxe_rebuild_buccaneer' : True,
                    'use_arpwarp' : False,
                    'use_buccaneer' : False,
        }
        options._preset_options("devel_mode")
        self.assertEqual(options.d, expected)
        
    def test_preset_options_webserver_mode(self):
        options = config_util.AMPLEConfigOptions()
        options.d = {'purge': False,
                     'cluster_method' : 'spicker',
                     'nproc': None,
                     'shelxe_rebuild_buccaneer': None,
                     'shelxe_rebuild_arpwarp': None,
                     'submit_cluster' : False,
                     'submit_max_array' : 20,
                     'submit_qtype' : None,
        }
        options.cmdline_opts = {}
        expected = {'purge': True,
                    'cluster_method' : 'spicker_tm',
                    'nproc': 1,
                    'shelxe_rebuild_buccaneer': True,
                    'shelxe_rebuild_arpwarp': False,
                    'submit_cluster' : True,
                    'submit_max_array' : 10,
                    'submit_qtype' : "SGE",
                    'submit_queue' : "all.q",
        }
        options._preset_options("webserver_uri")
        self.assertEqual(options.d, expected)
    
    def test_read_config_file(self):
        #Test reading the config file
        
        options = config_util.AMPLEConfigOptions()
         
        f = tempfile.NamedTemporaryFile("w", delete=False)
        f.write("[Databases]" + os.linesep)
        f.write("nr : nr_database" + os.linesep)
        f.write("rosetta_db : rosetta_database" + os.linesep + os.linesep)
        f.write("[Executables]" + os.linesep)
        f.write("blast_dir : blast_dir" + os.linesep)
        f.write("fast_protein_cluster_exe : foo.exe" + os.linesep)
        f.write("gesamt_exe : gesamt.exe" + os.linesep)
        f.write("maxcluster_exe : maxcluster.exe" + os.linesep)
        f.write("mustang_exe : mustang.exe" + os.linesep)
        f.write("rosetta_dir : rosetta_directory" + os.linesep)
        f.write("rosetta_fragments_exe : bar.exe" + os.linesep)
        f.write("rosetta_AbinitioRelax :rosetta_ab_relax" + os.linesep)
        f.write("scwrl_exe : scwrl.exe" + os.linesep)
        f.write("shelxe_exe : shelxe.exe" + os.linesep)
        f.write("spicker_exe : spicker.exe" + os.linesep)
        f.write("theseus_exe : theseus.exe" + os.linesep + os.linesep)
        f.write("[Files]" + os.linesep)
        f.write("alignment_file : foo.ali" + os.linesep)
        f.write("ample_log : foo.log" + os.linesep)
        f.write("bbcontacts_file : foo.bar" + os.linesep)
        f.write("cluster_dir : foo/bar" + os.linesep)
        f.write("config_file : foo.bar" + os.linesep)
        f.write("fasta : foo.fasta" + os.linesep)
        f.write("frags_3mers : foo.3mers" + os.linesep)
        f.write("frags_9mers : bar.9mers" + os.linesep)
        f.write("models : models" + os.linesep)
        f.write("mtz : foo.mtz" + os.linesep)
        f.write("native_pdb : foo.pdb" + os.linesep)
        f.write("psipred_ss2 : None" + os.linesep)
        f.write("restart_pkl : False" + os.linesep)
        f.write("score_matrix : True")
        f.close()
        
        options._read_config_file(f.name)
        
        expected = {
                    'nr' : os.path.abspath('nr_database'),
                    'rosetta_db' : os.path.abspath('rosetta_database'),
                    'blast_dir' : os.path.abspath('blast_dir'),
                    'fast_protein_cluster_exe' : os.path.abspath('foo.exe'),
                    'gesamt_exe' : os.path.abspath('gesamt.exe'),
                    'maxcluster_exe' : os.path.abspath('maxcluster.exe'),
                    'mustang_exe' : os.path.abspath('mustang.exe'),
                    'rosetta_dir' : os.path.abspath('rosetta_directory'),
                    'rosetta_fragments_exe' : os.path.abspath('bar.exe'),
                    'rosetta_AbinitioRelax' :os.path.abspath('rosetta_ab_relax'),
                    'scwrl_exe' : os.path.abspath('scwrl.exe'),
                    'shelxe_exe' : os.path.abspath('shelxe.exe'),
                    'spicker_exe' : os.path.abspath('spicker.exe'),
                    'theseus_exe' : os.path.abspath('theseus.exe'),
                    'alignment_file' : os.path.abspath('foo.ali'),
                    'ample_log' : os.path.abspath('foo.log'),
                    'bbcontacts_file' : os.path.abspath('foo.bar'),
                    'cluster_dir' : os.path.abspath('foo/bar'),
                    'config_file' : os.path.abspath('foo.bar'),
                    'fasta' : os.path.abspath('foo.fasta'),
                    'frags_3mers' : os.path.abspath('foo.3mers'),
                    'frags_9mers' : os.path.abspath('bar.9mers'),
                    'models' : os.path.abspath('models'),
                    'mtz' : os.path.abspath('foo.mtz'),
                    'native_pdb' : os.path.abspath('foo.pdb'),
                    'psipred_ss2' : None,
                    'restart_pkl' : False,
                    'score_matrix' : True
                    }
        
        self.assertItemsEqual(options.d, expected)
  
    def test_roundtrip_config_file(self):
        """Test we can read in and write out a config file with no changes"""
        options = config_util.AMPLEConfigOptions()
         
        
        # All files are run through abspath on reading by the config parser so we need to convert the names
        # We can't test [Files] as these are hashed out on writing and ignored on reading
        inputd = { 'nr' : os.path.abspath('nr'),
                   'rosetta_db' : os.path.abspath('rosetta_db'),
                   'blast_dir' : os.path.abspath('blast_dir'),
                   'fast_protein_cluster_exe' : os.path.abspath('fpc.exe'),
                  }
        input_str = """[AMPLE_info]

[Databases]
nr = {nr}
rosetta_db = {rosetta_db}

[Executables]
blast_dir = {blast_dir}
fast_protein_cluster_exe = {fast_protein_cluster_exe}

[Files]

[Molecular_Replacement]
F = None
FREE = None
arpwarp_cycles = 10
ensemble_options = None

[Unspecified]

""".format(**inputd)
        f1 = tempfile.NamedTemporaryFile("w", delete=False)
        f1.write(input_str)
        f1.close()
        
        # Read in the file from disk
        options._read_config_file(f1.name)
        
        f2 = tempfile.NamedTemporaryFile("w", delete=False)
        f2.close()
        
        # Write back out to disk
        options.write_config_file(config_file=f2.name)
        
        # Read in and make sure the strings are the same
        with open(f2.name) as f2: output_str = f2.read()
        
        self.assertItemsEqual(output_str, input_str)
        os.unlink(f1.name)
        os.unlink(f2.name)
        
    def test_read_config_opts(self):
        #Test read config options
        
        options = config_util.AMPLEConfigOptions()
        config_opts = {
                       'early_terminate': False,
                       'benchmark_mode': True,
                       'shelxe_rebuild' : True,
                       'shelxe_rebuild_arpwarp' : None,
                       'shelxe_rebuild_buccaneer' : True,
                       'use_arpwarp' : False,
                       'use_buccaneer' : False,
        }
        
        options.d = {
                     'ensemble_max_models' : 10,
                     'nmodels' : 400,
                     'percent' : 35,
                     'shelx_cycles' : 5,
                     'use_arpwarp' : True,
                     'use_buccaneer' : None,
                     'phaser_kill' : 15
        }
        
        expected = {
                     'nmodels' : 400,
                     'shelxe_rebuild_arpwarp' : None,
                     'cmdline_flags' : sorted(['early_terminate',
                                               'benchmark_mode',
                                               'shelxe_rebuild',
                                               'shelxe_rebuild_buccaneer',
                                               'use_arpwarp',
                                               'use_buccaneer']),
                     'percent' : 35,
                     'shelx_cycles' : 5,    
                     'phaser_kill' : 15,
                     'early_terminate': False,
                     'benchmark_mode': True,
                     'shelxe_rebuild' : True,
                     'ensemble_max_models' : 10,
                     'shelxe_rebuild_arpwarp' : None,
                     'shelxe_rebuild_buccaneer' : True,
                     'use_arpwarp' : False,
                     'use_buccaneer' : False,
        }
        
        options._read_cmdline_opts(config_opts)
        self.assertItemsEqual(options.d, expected)
    
    def test_isfloat(self):
        #Test the _isfloat function
        
        options = config_util.AMPLEConfigOptions()
        value = 75
        expected = True
        self.assertEqual(expected, options._isfloat(value))
        
        value = 120840287
        expected = True
        self.assertEqual(expected, options._isfloat(value))
        
        value = 10120.0
        expected = True
        self.assertEqual(expected, options._isfloat(value))
        
        value = 0.0
        expected = True
        self.assertEqual(expected, options._isfloat(value))
        
        value = 190820198409384039285
        expected = True
        self.assertEqual(expected, options._isfloat(value))
        
        value = "not a number"
        expected = False
        self.assertEqual(expected, options._isfloat(value))
        
if __name__ == "__main__":
    unittest.main()
