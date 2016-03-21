"""Test functions for util.ensembler_util"""

import itertools
import random
import os
import unittest
from ample.util import ensembler_util

class Test(unittest.TestCase):
        
    def test_summary(self):
        ensembles_data = []
        for i in xrange(1, 101, 5):
            d = {'truncation_variance' : round(i*3.54*1.25, 1),
                 'truncation_pruning': None,
                 'truncation_method': 'percent',
                 'subcluster_centroid_model': os.path.join(os.sep, 'foo', 'bar', 
                                                           'TEST_{0}.pdb'.format(i)),
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
                 'ensemble_pdb': os.path.join(os.sep, 'foo', 'bar', 
                                              'TEST_ENS_{0}.pdb'.format(i)),
                 'truncation_residues': [i for i in xrange(1, 21)],
                 'subcluster_radius_threshold': 1,
                 'cluster_num': 1,
                 'pruned_residues': None,
                 'cluster_centroid': os.path.join(os.sep, 'foo', 
                                                  'bar', 'centroid.pdb')}
            ensembles_data.append(d)

        print ensembler_util.ensemble_summary(ensembles_data)

        return

    def test_importModule(self):
        ########################################################################
        # Test Case 1
        module_handler = ensembler_util.find_ensembler_module({'homologs': False,
                                                               'single_model_mode': False})
        self.assertEqual("ample.ensembler.abinitio", module_handler.__name__)

        ########################################################################
        # Test Case 2
        module_handler = ensembler_util.find_ensembler_module({'homologs': True,
                                                               'single_model_mode': False})
        self.assertEqual("ample.ensembler.homologs", module_handler.__name__)

        ########################################################################
        # Test Case 3
        module_handler = ensembler_util.find_ensembler_module({'homologs': False,
                                                               'single_model_mode': True})
        self.assertEqual("ample.ensembler.single_model", module_handler.__name__)

        return

    def test_getEnsemblingKwargs(self):
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
        kwargs = ensembler_util._get_ensembling_kwargs(amoptd_fake)
        self.assertEqual(sorted(ref_keys), sorted(kwargs.keys()))

        ########################################################################
        # Test Case 2 - homologs and no single model
        ref_keys = common_keys + homolog_keys
        amoptd_fake.update({'homologs': True, 'single_model_mode': False})
        kwargs = ensembler_util._get_ensembling_kwargs(amoptd_fake)
        self.assertEqual(sorted(ref_keys), sorted(kwargs.keys()))

        ########################################################################
        # Test Case 3 - no homologs but single model
        ref_keys = common_keys + single_model_keys
        amoptd_fake.update({'homologs': False, 'single_model_mode': True})
        kwargs = ensembler_util._get_ensembling_kwargs(amoptd_fake)
        self.assertEqual(sorted(ref_keys), sorted(kwargs.keys()))

        return

    def test_sortEnsembles1(self):
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
            pdb = os.path.join(os.sep, "foo", "bar", name+".pdb")
            ensemble = {'name' : name,
                        'ensemble_pdb': pdb,
                        'cluster_num': comb[0],
                        'truncation_level': comb[1],
                        'subcluster_radius_threshold': comb[2],
                        'side_chain_treatment': comb[3],
            }

            ensemble_data.append(ensemble)
            ensemble_pdbs.append(pdb)

        _root = os.path.join(os.sep, "foo", "bar")  
        ########################################################################       
        # Test Case 1    
        ensemble_pdb_sorted = ensembler_util.sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=True)
        self.assertEqual(os.path.join(_root, "c1_tl20_r1_allatom.pdb"), ensemble_pdb_sorted[0])
        self.assertEqual(os.path.join(_root, "c1_tl19_r1_allatom.pdb"), ensemble_pdb_sorted[18])
        self.assertEqual(os.path.join(_root, "c1_tl100_r3_reliable.pdb"), ensemble_pdb_sorted[44])
        self.assertEqual(os.path.join(_root, "c2_tl20_r1_allatom.pdb"), ensemble_pdb_sorted[45])
        self.assertEqual(os.path.join(_root, "c2_tl19_r1_allatom.pdb"), ensemble_pdb_sorted[63])
        self.assertEqual(os.path.join(_root, "c2_tl100_r3_reliable.pdb"), ensemble_pdb_sorted[89])
        self.assertEqual(os.path.join(_root, "c3_tl20_r1_allatom.pdb"), ensemble_pdb_sorted[90])
        self.assertEqual(os.path.join(_root, "c3_tl19_r1_allatom.pdb"), ensemble_pdb_sorted[108])
        self.assertEqual(os.path.join(_root, "c3_tl100_r3_reliable.pdb"), ensemble_pdb_sorted[-1])
        ########################################################################       
        # Test Case 2
        ensemble_pdb_sorted = ensembler_util.sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=False)
        self.assertEqual(os.path.join(_root, "c1_tl19_r1_allatom.pdb"), ensemble_pdb_sorted[0])
        self.assertEqual(os.path.join(_root, "c3_tl100_r3_reliable.pdb"), ensemble_pdb_sorted[-1])
        self.assertEqual(os.path.join(_root, "c2_tl50_r2_polyAla.pdb"), ensemble_pdb_sorted[67])
        ########################################################################       
        # Test Case 3
        ensemble_pdb_sorted = ensembler_util.sort_ensembles(ensemble_pdbs)
        self.assertEqual(os.path.join(_root, "c1_tl100_r1_allatom.pdb"), ensemble_pdb_sorted[0])
        self.assertEqual(os.path.join(_root, "c3_tl80_r3_reliable.pdb"), ensemble_pdb_sorted[-1])

        return
    
    def test_sortEnsembles2(self):
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

            pdb = os.path.join(os.sep, "foo", "bar", name+".pdb")
            ensemble = {'name' : name,
                        'ensemble_pdb': pdb,
                        'truncation_level': comb[0],
                        'side_chain_treatment': comb[1],
            }

            ensemble_data.append(ensemble)
            ensemble_pdbs.append(pdb)

        _root = os.path.join(os.sep, "foo", "bar")
        ########################################################################       
        # Test Case 1
        ensemble_pdb_sorted = ensembler_util.sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=True)
        self.assertEqual(os.path.join(_root, "e20_allatom.pdb"), ensemble_pdb_sorted[0])
        self.assertEqual(os.path.join(_root, "e19_polyAla.pdb"), ensemble_pdb_sorted[len(ensemble_pdb_sorted)/2])
        self.assertEqual(os.path.join(_root, "e100_reliable.pdb"), ensemble_pdb_sorted[-1])
        ########################################################################       
        # Test Case 2
        ensemble_pdb_sorted = ensembler_util.sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=False)
        self.assertEqual(os.path.join(_root, "e19_allatom.pdb"), ensemble_pdb_sorted[0])
        self.assertEqual(os.path.join(_root, "e100_reliable.pdb"), ensemble_pdb_sorted[-1])
        self.assertEqual(os.path.join(_root, "e50_polyAla.pdb"), ensemble_pdb_sorted[len(ensemble_pdb_sorted)/2])
        ########################################################################       
        # Test Case 3
        ensemble_pdb_sorted = ensembler_util.sort_ensembles(ensemble_pdbs)
        self.assertEqual(os.path.join(_root, "e100_allatom.pdb"), ensemble_pdb_sorted[0])
        self.assertEqual(os.path.join(_root, "e80_reliable.pdb"), ensemble_pdb_sorted[-1])

        return
    
    def test_sortEnsembles3(self):
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
            pdb = os.path.join(os.sep, "foo", "bar", name+".pdb")
            ensemble = {'name' : name,
                        'ensemble_pdb': pdb,
                        'truncation_score_key': comb[0],
                        'truncation_level': comb[1],
                        'side_chain_treatment': comb[2],
            }

            ensemble_data.append(ensemble)
            ensemble_pdbs.append(pdb)

        _root = os.path.join(os.sep, "foo", "bar")
        ########################################################################       
        # Test Case 1
        ensemble_pdb_sorted = ensembler_util.sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=True)
        self.assertEqual(os.path.join(_root, "bbc_tl20_allatom.pdb"), ensemble_pdb_sorted[0])
        self.assertEqual(os.path.join(_root, "bbc_tl19_allatom.pdb"), ensemble_pdb_sorted[6])
        self.assertEqual(os.path.join(_root, "bbc_tl100_reliable.pdb"), ensemble_pdb_sorted[14])
        self.assertEqual(os.path.join(_root, "kicker_tl20_allatom.pdb"), ensemble_pdb_sorted[15])
        self.assertEqual(os.path.join(_root, "kicker_tl19_allatom.pdb"), ensemble_pdb_sorted[21])
        self.assertEqual(os.path.join(_root, "kicker_tl100_reliable.pdb"), ensemble_pdb_sorted[29])
        self.assertEqual(os.path.join(_root, "ntv_tl20_allatom.pdb"), ensemble_pdb_sorted[30])
        self.assertEqual(os.path.join(_root, "ntv_tl19_allatom.pdb"), ensemble_pdb_sorted[36])
        self.assertEqual(os.path.join(_root, "ntv_tl100_reliable.pdb"), ensemble_pdb_sorted[-1])
        ########################################################################       
        # Test Case 2
        ensemble_pdb_sorted = ensembler_util.sort_ensembles(ensemble_pdbs, ensemble_data, prioritise=False)
        self.assertEqual(os.path.join(_root, "bbc_tl19_allatom.pdb"), ensemble_pdb_sorted[0])
        self.assertEqual(os.path.join(_root, "ntv_tl100_reliable.pdb"), ensemble_pdb_sorted[-1])
        self.assertEqual(os.path.join(_root, "kicker_tl50_polyAla.pdb"), ensemble_pdb_sorted[len(ensemble_pdb_sorted)/2])
        ########################################################################       
        # Test Case 3
        ensemble_pdb_sorted = ensembler_util.sort_ensembles(ensemble_pdbs)
        self.assertEqual(os.path.join(_root, "bbc_tl100_allatom.pdb"), ensemble_pdb_sorted[0])
        self.assertEqual(os.path.join(_root, "ntv_tl80_reliable.pdb"), ensemble_pdb_sorted[-1])

        return

if __name__ == "__main__":
    unittest.main()
