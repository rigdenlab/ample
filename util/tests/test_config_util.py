"""Tests for util.config_util"""

import os
import tempfile
import unittest
from ample.util import config_util
from ample.util import version

class TestCases(unittest.TestCase):
    MAX_DIFF = None
    
    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        paths = cls.thisd.split(os.sep)
        cls.ample_dir = os.sep.join(paths[ :-2 ])
        cls.tests_dir = os.path.join(cls.ample_dir, "tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir, 'testfiles')
    
    def test_process_options(self):
        #Test the process_options
        options = config_util.AMPLEConfigOptions()
        options.d = {
                     'fasta' : os.path.join(self.testfiles_dir, '2uui.fasta'),
                     'native_pdb' : os.path.join(self.testfiles_dir, '2UUI.pdb'),
                     'rcdir' : 'foo/bar',
                     'side_chain_treatments' : False,
                     'submit_qtype' : "sge",
                     'shelxe_rebuild' : True,
                     'devel_mode' : False,
                     'quick_mode' : False,
                     'webserver_uri' : False
                     }
        
    
        expected = {
                    'ample_version' : version.__version__,
                    'fasta' : os.path.join(self.testfiles_dir, '2uui.fasta'),
                    'native_pdb' : os.path.join(self.testfiles_dir, '2UUI.pdb'),
                    'rcdir' : 'foo/bar',
                    'side_chain_treatments' : ['polyAla', 'reliable', 'allatom'],
                    'submit_qtype' : "SGE",
                    'shelxe_rebuild' : True,
                    'shelxe_rebuild_arpwap' : True,
                    'shelxe_rebuild_buccaneer' : True,
                    'devel_mode' : False,
                    'quick_mode' : False,
                    'webserver_uri' : False
                    }
        
        options._process_options()
        self.assertItemsEqual(options.d, expected)
        
    def test_preset_options(self):
        #Test the preset options
        
        options = config_util.AMPLEConfigOptions()
        #quick_mode test
        options.d = {
                     'max_ensemble_models' : 10,
                     'nmodels' : 400,
                     'percent' : 35,
                     'shelx_cycles' : 5,
                     'use_arpwarp' : True,
                     'use_buccaneer' : None,
                     'phaser_kill' : 15
        }
        
        options.config_opts = {
                               'max_ensemble_models' : 10,
                               'nmodels' : 200,
                               'percent' : 20,
                               'shelx_cycles' : 5,
                               'use_arpwarp' : False,
                               'use_buccaneer' : False,
                               'phaser_kill' : 15
        }
        
        expected = {
                    'max_ensemble_models' : 10,
                    'nmodels' : 400,
                    'percent' : 35,
                    'shelx_cycles' : 5,
                    'use_arpwarp' : True,
                    'use_buccaneer' : False,
                    'phaser_kill' : 15
        }
        
        options._preset_options("quick_mode")
        self.assertEqual(options.d, expected)
        
        #devel_mode test
        options.d = {
                     'early_terminate': True,
                     'benchmark_mode': True,
                     'shelxe_rebuild' : None,
                     'shelxe_rebuild_arpwarp' : True,
                     'shelxe_rebuild_buccaneer' : True,
                     'use_arpwarp' : None,
                     'use_buccaneer' : False,
        }
        
        options.config_opts = {
                              'early_terminate': False,
                              'benchmark_mode': True,
                              'shelxe_rebuild' : True,
                              'shelxe_rebuild_arpwarp' : True,
                              'shelxe_rebuild_buccaneer' : True,
                              'use_arpwarp' : False,
                              'use_buccaneer' : False,
        }
        
        expected = {
                    'early_terminate': True,
                    'benchmark_mode': True,
                    'shelxe_rebuild' : True,
                    'shelxe_rebuild_arpwarp' : True,
                    'shelxe_rebuild_buccaneer' : True,
                    'use_arpwarp' : False,
                    'use_buccaneer' : False,
        }
        
        options._preset_options("devel_mode")
        self.assertEqual(options.d, expected)
        
        #webserver_uri test
        options.d = {
                     'purge': True,
                     'shelxe_rebuild_buccaneer': None,
                     'submit_cluster' : False,
                     'submit_max_array' : 20,
                     'submit_qtype' : None,
                     'submit_queue' : "all.q",
        }
        
        options.config_opts = {
                               'purge': True,
                               'shelxe_rebuild_buccaneer': True,
                               'submit_cluster' : True,
                               'submit_max_array' : 10,
                               'submit_qtype' : "SGE",
                               'submit_queue' : "all.q",
        }
        
        expected = {
                    'purge': True,
                    'shelxe_rebuild_buccaneer': True,
                    'submit_cluster' : False,
                    'submit_max_array' : 20,
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
                    'rosetta_abinitiorelax' :os.path.abspath('rosetta_ab_relax'),
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
                     'max_ensemble_models' : 10,
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
                     'max_ensemble_models' : 10,
                     'shelxe_rebuild_arpwarp' : None,
                     'shelxe_rebuild_buccaneer' : True,
                     'use_arpwarp' : False,
                     'use_buccaneer' : False,
        }
        
        options._read_config_opts(config_opts)
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
