"""Test functions for ensembler.ensembler_util"""

import unittest

from ample.ensembler import ensembler_util

class Test(unittest.TestCase):
        
    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        paths = cls.thisd.split(os.sep)
        cls.ample_dir = os.sep.join(paths[ :-1 ])
        cls.tests_dir = os.path.join(cls.ample_dir, "tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir, 'testfiles')
        cls.theseus_exe = ample_util.find_exe('theseus')
        cls.spicker_exe = ample_util.find_exe('spicker')
        cls.maxcluster_exe = ample_util.find_exe('maxcluster')

        return
    
    def testSummary(self):
        """Not really a unittest but prints a table for 
           visual analysis of correctness
        """

        ensembles_data = []
        for i in xrange(1, 101, 5):
            d = {'truncation_variance' : round(i*3.54*1.25, 1),
                 'truncation_pruning': None,
                 'truncation_method': 'percent',
                 'subcluster_centroid_model': '/foo/bar/TEST_{0}.pdb'.format(i),
                 'ensemble_num_atoms': 100,
                 'truncation_level': i,
                 'truncation_dir': '/foo/bar/TEST',
                 'cluster_method': 'spicker',
                 'num_clusters': 1,
                 'num_residues': 20,
                 'percent_truncation': '10',
                 'name': 'TEST_{0}'.format(i),
                 'side_chain_treatment': 'TEST',
                 'cluster_num_models': 31,
                 'subcluster_num_models': 3,
                 'ensemble_pdb': '/foo/bar/TEST_ENS_{0}.pdb'.format(i),
                 'truncation_residues': [i for i in xrange(1, 21)],
                 'subcluster_radius_threshold': 1,
                 'cluster_num': 1,
                 'pruned_residues': None,
                 'cluster_centroid': '/foo/bar/centroid.pdb'}
            ensembles_data.append(d)

        print ensemble_summary(ensembles_data)

        return

    def testImportModule(self):
        """Test that we import the correct module for ensembling"""

        ########################################################################
        # Test Case 1
        module_handler = find_ensembler_module({'homologs': False,
                                                'single_model_mode': False})
        print module_handler
        self.assertEqual("ample.ensembler.abinitio", module_handler.__name__)

        ########################################################################
        # Test Case 2
        module_handler = find_ensembler_module({'homologs': True,
                                                'single_model_mode': False})
        self.assertEqual("ample.ensembler.homologs", module_handler.__name__)

        ########################################################################
        # Test Case 3
        module_handler = find_ensembler_module({'homologs': False,
                                                'single_model_mode': True})
        self.assertEqual("ample.ensembler.single_model", module_handler.__name__)

        return

    def testGetEnsemblingKwargs(self):
        """Test that the keyword arguments for ensembling are set correctly"""

        # All keys found in return kwargs written here in corresponding sections
        common_keys = ['ensembles_directory', 'nproc', 'percent_truncation',
                       'side_chain_treatments', 'truncation_method']

        abinitio_keys = ['cluster_dir', 'cluster_exe', 'cluster_method',
                         'num_clusters', 'truncation_pruning', 'use_scwrl']

        homolog_keys = ['alignment_file', 'gesamt_exe', 'homolog_aligner',
                        'mustang_exe']

        single_model_keys = ['truncation_pruning', 'truncation_scorefile',
                             'truncation_scorefile_header']

        # Generate a pseudo-AMPLE dictionary - values can be None as not required
        amoptd_fake = { key: None for key in common_keys + abinitio_keys + \
                                             homolog_keys + single_model_keys }
        amoptd_fake.update({'homologs': None,
                            'percent': None,
                            'single_model_mode': None})

        ########################################################################
        # Test Case 1 - ab initio
        ref_keys = common_keys + abinitio_keys
        kwargs = _get_ensembling_kwargs(amoptd_fake)
        self.assertEqual(sorted(ref_keys), sorted(kwargs.keys()))

        ########################################################################
        # Test Case 2 - homologs and no single model
        ref_keys = common_keys + homolog_keys
        amoptd_fake.update({'homologs': True, 'single_model_mode': False})
        kwargs = _get_ensembling_kwargs(amoptd_fake)
        self.assertEqual(sorted(ref_keys), sorted(kwargs.keys()))

        ########################################################################
        # Test Case 3 - no homologs but single model
        ref_keys = common_keys + single_model_keys
        amoptd_fake.update({'homologs': False, 'single_model_mode': True})
        kwargs = _get_ensembling_kwargs(amoptd_fake)
        self.assertEqual(sorted(ref_keys), sorted(kwargs.keys()))

        return

    def test_sortEnsembles1(self):
        """Test sorting for ab initio data"""

        clusters = [1, 2, 3]
        tlevels = [19, 20, 50, 80, 100]
        radii = [1, 2, 3]
        sidechains = ["allatom", "polyAla", "reliable"]

        # Shuffle the data so the order of the combinations is random
        # Although undesirable in a test, exactly what we need as sorted outcome
        # should ALWAYS be the same
        random.shuffle(clusters)
        random.shuffle(tlevels)
        random.shuffle(radii)
        random.shuffle(sidechains)

        ensemble_data = []
        ensemble_pdbs = []

        # Generate ensembles to sort 
        for comb in list(itertools.product(clusters, tlevels, radii, sidechains)):

            name = "c{cluster}_tl{truncation}_r{radius}_{schain}".format(cluster=comb[0],
                                                                         truncation=comb[1],
                                                                         radius=comb[2],
                                                                         schain=comb[3])
            pdb = os.path.join("/foo", "bar", name+".pdb")
            ensemble = {'name' : name,
                        'ensemble_pdb': pdb,
                        'cluster_num': comb[0],
                        'truncation_level': comb[1],
                        'subcluster_radius_threshold': comb[2],
                        'side_chain_treatment': comb[3],
            }

            ensemble_data.append(ensemble)
            ensemble_pdbs.append(pdb)

        ########################################################################       
        # Test Case 1    
        ensemble_pdb_sorted = sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=True)
        self.assertEqual("/foo/bar/c1_tl20_r1_allatom.pdb", ensemble_pdb_sorted[0])
        self.assertEqual("/foo/bar/c1_tl19_r1_allatom.pdb", ensemble_pdb_sorted[18])
        self.assertEqual("/foo/bar/c1_tl100_r3_reliable.pdb", ensemble_pdb_sorted[44])
        self.assertEqual("/foo/bar/c2_tl20_r1_allatom.pdb", ensemble_pdb_sorted[45])
        self.assertEqual("/foo/bar/c2_tl19_r1_allatom.pdb", ensemble_pdb_sorted[63])
        self.assertEqual("/foo/bar/c2_tl100_r3_reliable.pdb", ensemble_pdb_sorted[89])
        self.assertEqual("/foo/bar/c3_tl20_r1_allatom.pdb", ensemble_pdb_sorted[90])
        self.assertEqual("/foo/bar/c3_tl19_r1_allatom.pdb", ensemble_pdb_sorted[108])
        self.assertEqual("/foo/bar/c3_tl100_r3_reliable.pdb", ensemble_pdb_sorted[-1])
        ########################################################################       
        # Test Case 2
        ensemble_pdb_sorted = sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=False)
        self.assertEqual("/foo/bar/c1_tl19_r1_allatom.pdb", ensemble_pdb_sorted[0])
        self.assertEqual("/foo/bar/c3_tl100_r3_reliable.pdb", ensemble_pdb_sorted[-1])
        self.assertEqual("/foo/bar/c2_tl50_r2_polyAla.pdb", ensemble_pdb_sorted[67])
        ########################################################################       
        # Test Case 3
        ensemble_pdb_sorted = sort_ensembles(ensemble_pdbs)
        self.assertEqual("/foo/bar/c1_tl100_r1_allatom.pdb", ensemble_pdb_sorted[0])
        self.assertEqual("/foo/bar/c3_tl80_r3_reliable.pdb", ensemble_pdb_sorted[-1])

        return
    
    def test_sortEnsembles2(self):
        """Test sorting for homolog data"""

        tlevels = [19, 20, 50, 80, 100]
        sidechains = ["allatom", "polyAla", "reliable"]

        # Shuffle the data so the order of the combinations is random
        # Although undesirable in a test, exactly what we need as sorted outcome
        # should ALWAYS be the same
        random.shuffle(tlevels)
        random.shuffle(sidechains)

        ensemble_data = []
        ensemble_pdbs = []

        # Generate ensembles to sort 
        for comb in list(itertools.product(tlevels, sidechains)):

            name = "e{truncation}_{schain}".format(truncation=comb[0],
                                                   schain=comb[1])

            pdb = os.path.join("/foo", "bar", name+".pdb")
            ensemble = {'name' : name,
                        'ensemble_pdb': pdb,
                        'truncation_level': comb[0],
                        'side_chain_treatment': comb[1],
            }

            ensemble_data.append(ensemble)
            ensemble_pdbs.append(pdb)

        ########################################################################       
        # Test Case 1
        ensemble_pdb_sorted = sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=True)
        self.assertEqual("/foo/bar/e20_allatom.pdb", ensemble_pdb_sorted[0])
        self.assertEqual("/foo/bar/e19_polyAla.pdb", ensemble_pdb_sorted[len(ensemble_pdb_sorted)/2])
        self.assertEqual("/foo/bar/e100_reliable.pdb", ensemble_pdb_sorted[-1])
        ########################################################################       
        # Test Case 2
        ensemble_pdb_sorted = sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=False)
        self.assertEqual("/foo/bar/e19_allatom.pdb", ensemble_pdb_sorted[0])
        self.assertEqual("/foo/bar/e100_reliable.pdb", ensemble_pdb_sorted[-1])
        self.assertEqual("/foo/bar/e50_polyAla.pdb", ensemble_pdb_sorted[len(ensemble_pdb_sorted)/2])
        ########################################################################       
        # Test Case 3
        ensemble_pdb_sorted = sort_ensembles(ensemble_pdbs)
        self.assertEqual("/foo/bar/e100_allatom.pdb", ensemble_pdb_sorted[0])
        self.assertEqual("/foo/bar/e80_reliable.pdb", ensemble_pdb_sorted[-1])

        return
    
    def test_sortEnsembles3(self):
        '''Test sorting for single model data'''

        score_keys = ["ntv", "bbc", "kicker"]
        tlevels = [19, 20, 50, 80, 100]
        sidechains = ["allatom", "polyAla", "reliable"]

        # Shuffle the data so the order of the combinations is random
        # Although undesirable in a test, exactly what we need as sorted outcome
        # should ALWAYS be the same
        random.shuffle(score_keys)
        random.shuffle(tlevels)
        random.shuffle(sidechains)

        ensemble_data = []
        ensemble_pdbs = []

        # Generate ensembles to sort 
        for comb in list(itertools.product(score_keys, tlevels, sidechains)):

            name = "{score_key}_tl{truncation}_{schain}".format(score_key=comb[0],
                                                                truncation=comb[1],
                                                                schain=comb[2])
            pdb = os.path.join("/foo", "bar", name+".pdb")
            ensemble = {'name' : name,
                        'ensemble_pdb': pdb,
                        'truncation_score_key': comb[0],
                        'truncation_level': comb[1],
                        'side_chain_treatment': comb[2],
            }

            ensemble_data.append(ensemble)
            ensemble_pdbs.append(pdb)

        ########################################################################       
        # Test Case 1
        ensemble_pdb_sorted = sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=True)
        self.assertEqual("/foo/bar/bbc_tl20_allatom.pdb", ensemble_pdb_sorted[0])
        self.assertEqual("/foo/bar/bbc_tl19_allatom.pdb", ensemble_pdb_sorted[6])
        self.assertEqual("/foo/bar/bbc_tl100_reliable.pdb", ensemble_pdb_sorted[14])
        self.assertEqual("/foo/bar/kicker_tl20_allatom.pdb", ensemble_pdb_sorted[15])
        self.assertEqual("/foo/bar/kicker_tl19_allatom.pdb", ensemble_pdb_sorted[21])
        self.assertEqual("/foo/bar/kicker_tl100_reliable.pdb", ensemble_pdb_sorted[29])
        self.assertEqual("/foo/bar/ntv_tl20_allatom.pdb", ensemble_pdb_sorted[30])
        self.assertEqual("/foo/bar/ntv_tl19_allatom.pdb", ensemble_pdb_sorted[36])
        self.assertEqual("/foo/bar/ntv_tl100_reliable.pdb", ensemble_pdb_sorted[-1])
        ########################################################################       
        # Test Case 2
        ensemble_pdb_sorted = sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=False)
        self.assertEqual("/foo/bar/bbc_tl19_allatom.pdb", ensemble_pdb_sorted[0])
        self.assertEqual("/foo/bar/ntv_tl100_reliable.pdb", ensemble_pdb_sorted[-1])
        self.assertEqual("/foo/bar/kicker_tl50_polyAla.pdb", ensemble_pdb_sorted[len(ensemble_pdb_sorted)/2])
        ########################################################################       
        # Test Case 3
        ensemble_pdb_sorted = sort_ensembles(ensemble_pdbs)
        self.assertEqual("/foo/bar/bbc_tl100_allatom.pdb", ensemble_pdb_sorted[0])
        self.assertEqual("/foo/bar/ntv_tl80_reliable.pdb", ensemble_pdb_sorted[-1])

        return
    