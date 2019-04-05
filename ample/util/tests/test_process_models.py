"""Test facility for ample/util/check_models.py"""

__author__ = "Jens Thomas"

import glob
import iotbx.pdb
import os
import shutil
import tempfile
import unittest

from ample import constants
from ample.util import process_models

import logging
logging.basicConfig(level=logging.DEBUG)


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath(os.path.dirname(__file__)) 
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share,'testfiles')

    def test_check_models_nmr(self):
        pdbin = os.path.join(self.ample_share, 'examples', 'nmr-truncate', 'input', '2LC9.pdb')
        results = process_models.CheckModelsResult()
        results.models_dir = os.curdir
        process_models.check_models([pdbin], results)
        
        self.assertIsNone(results.error)
        self.assertTrue(results.ensemble)
        self.assertEquals(results.sequence,
                          'MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAAAINMVFQMGETGVAAFTNSLPMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL')
        self.assertEqual(results.num_models, 10)
        self.assertEqual(results.num_structures, 1)
        self.assertEqual(results.single_ensemble, True)
        
        self.assertFalse(results.created_updated_models)
        self.assertFalse(results.homologs)
        self.assertFalse(results.merged_chains)

    def test_check_models_single_structure(self):
        pdbin = os.path.join(self.testfiles_dir, '1K33_S_00000001.pdb')
        results = process_models.CheckModelsResult()
        results.models_dir = os.curdir
        process_models.check_models([pdbin], results)
        
        self.assertIsNone(results.error)
        self.assertEqual(results.num_models, 1)
        self.assertEqual(results.num_structures, 1)
        
        self.assertFalse(results.created_updated_models)
        self.assertFalse(results.homologs)
        self.assertFalse(results.merged_chains)
        self.assertFalse(results.ensemble)

    def test_check_models_single_structure_fail(self):
        """Should fail as has a non-protein chain present"""
        pdbin = os.path.join(self.testfiles_dir, '1K33.pdb')
        results = process_models.CheckModelsResult()
        results.models_dir = os.curdir
        process_models.check_models([pdbin], results)
        self.assertIsNotNone(results.error)
        

    def test_check_models_single_structure_merge_chains(self):
        pdbin = os.path.join(self.testfiles_dir, '1K33_S_00000001.pdb')
        outdir = tempfile.mkdtemp(dir=os.path.abspath(os.curdir))
        
        # write out pdb with no chain
        h = iotbx.pdb.pdb_input(pdbin).construct_hierarchy()
        for m in h.models():
            for c in m.chains():
                c.id = ''
        pdbout = os.path.join(outdir, 'nochain_id.pdb')
        with open(pdbout, 'w') as f:
            f.write("REMARK Original file:{}\n".format(pdbin))
            f.write(h.as_pdb_string(anisou=False))
        pdbin = pdbout
        
        results = process_models.CheckModelsResult()
        results.models_dir = outdir
        process_models.check_models([pdbin], results)
        
        self.assertIsNone(results.error)
        self.assertEqual(results.num_models, 1)
        self.assertEqual(results.num_structures, 1)
        self.assertTrue(results.created_updated_models)
        
        self.assertFalse(results.homologs)
        self.assertFalse(results.merged_chains)
        self.assertFalse(results.ensemble)
        shutil.rmtree(outdir)


    def test_check_models_multi_abinitio(self):
        pdbdir = os.path.join(self.testfiles_dir, 'models')
        outdir = tempfile.mkdtemp(dir=os.path.abspath(os.curdir))
        
        results = process_models.check_models_dir(pdbdir, outdir)
        self.assertEqual(results.num_models, 30)
        self.assertEqual(results.num_structures, 30)
        self.assertEquals(results.sequence, 'QPRRKLCILHRNPGRCYDKIPAFYYNQKKKQCERFDWSGCGGNSNRFKTIEECRRTCIG')

        self.assertIsNone(results.error)
        self.assertFalse(results.created_updated_models)
        self.assertFalse(results.homologs)
        self.assertFalse(results.merged_chains)
        self.assertFalse(results.ensemble)
        os.rmdir(outdir)
        

    def test_check_models_multi_abinitio_missing_chainid(self):
        pdbdir = os.path.join(self.testfiles_dir, 'models')
        indir = tempfile.mkdtemp(dir=os.path.abspath(os.curdir))
        outdir = tempfile.mkdtemp(dir=os.path.abspath(os.curdir))

        # write out pdbs with no chain
        for pdbin in glob.glob(os.path.join(pdbdir, '*.pdb')):
            h = iotbx.pdb.pdb_input(pdbin).construct_hierarchy()
            for m in h.models():
                for c in m.chains():
                    c.id = ''
            name = os.path.basename(pdbin)
            pdbout = os.path.join(indir, name)
            with open(pdbout, 'w') as f:
                f.write("REMARK Original file:{}\n".format(pdbin))
                f.write(h.as_pdb_string(anisou=False))

        results = process_models.check_models_dir(indir, outdir)
 
        self.assertIsNone(results.error)
        self.assertEqual(results.num_models, 30)
        self.assertEqual(results.num_structures, 30)
        self.assertTrue(results.created_updated_models)

        self.assertFalse(results.homologs)
        self.assertFalse(results.merged_chains)
        self.assertFalse(results.ensemble)
        shutil.rmtree(indir)        
        shutil.rmtree(outdir)        


    def test_check_models_multi_merge_chains(self):
        pdbname = '6gvl'
        pdbroot = os.path.join(self.testfiles_dir, pdbname + '.pdb')
        indir = tempfile.mkdtemp(dir=os.path.abspath(os.curdir))
        outdir = tempfile.mkdtemp(dir=os.path.abspath(os.curdir))

        for i in range(3):
            name = "{}_{}.pdb".format(pdbname, i)
            pdb = os.path.join(indir, name)
            shutil.copy(pdbroot, pdb)

        results = process_models.check_models_dir(indir, outdir)
   
        self.assertIsNone(results.error)
        self.assertTrue(results.created_updated_models)
        self.assertTrue(results.merged_chains)
        self.assertEqual(results.num_models, 3)
        self.assertEqual(results.num_structures, 3)

        self.assertFalse(results.homologs)
        self.assertFalse(results.ensemble)
        shutil.rmtree(indir)        
        shutil.rmtree(outdir) 


    def test_check_models_homologs(self):
        pdbdir = os.path.join(self.ample_share, 'examples', 'homologs', 'input')
        outdir = os.curdir
        
        results = process_models.check_models_dir(pdbdir, outdir)
        
        self.assertIsNone(results.error)
        self.assertTrue(results.homologs)
        self.assertEqual(results.num_models, 3)
        self.assertEqual(results.num_structures, 3)

        self.assertFalse(results.created_updated_models)
        self.assertFalse(results.merged_chains)
        self.assertFalse(results.ensemble)


    def test_extract_and_validate_models(self):
        mdir = os.path.join(self.testfiles_dir, 'models')
        amoptd = {'models': mdir,
                  'models_dir': mdir}
#         sequence = 'QPRRKLCILHRNPGRCYDKIPAFYYNQKKKQCERFDWSGCGGNSNRFKTIEECRRTCIG'
        process_models.extract_and_validate_models(amoptd)
        self.assertEqual(len(amoptd['processed_models']), 30)
        self.assertEqual(amoptd['models_dir'], mdir)


    def test_extract_and_validate_models_quark(self):
        temp_dir = tempfile.mkdtemp(dir=os.getcwd())
        amoptd = {'models' : os.path.join(self.testfiles_dir, 'result.tar.bz2'),
                  'models_dir' : temp_dir,
                  'work_dir' : '.'}
#         sequence = 'QPRRKLCILHRNPGRCYDKIPAFYYNQKKKQCERFDWSGCGGNSNRFKTIEECRRTCIG'
        process_models.extract_and_validate_models(amoptd)
        self.assertTrue(amoptd['quark_models'])
        self.assertEqual(len(amoptd['processed_models']), 200)
        shutil.rmtree(temp_dir)

        
    def test_extract_and_validate_models_quark_zip(self):
        temp_dir = tempfile.mkdtemp(dir=os.getcwd())
        amoptd = {'models' : os.path.join(self.testfiles_dir, 'result.zip'),
                  'models_dir' : temp_dir,
                  'work_dir' : '.'}
#         sequence = 'QPRRKLCILHRNPGRCYDKIPAFYYNQKKKQCERFDWSGCGGNSNRFKTIEECRRTCIG'
        process_models.extract_and_validate_models(amoptd)
        self.assertTrue(amoptd['quark_models'])
        self.assertEqual(len(amoptd['processed_models']), 200)
        shutil.rmtree(temp_dir)

        
    def test_extract_and_validate_models_quark_file(self):
        temp_dir = tempfile.mkdtemp(dir=os.getcwd())
        amoptd = {'models' : os.path.join(self.testfiles_dir, 'alldecoy.pdb'),
                  'models_dir' : temp_dir,
                  'work_dir' : '.'}
#         sequence = 'QPRRKLCILHRNPGRCYDKIPAFYYNQKKKQCERFDWSGCGGNSNRFKTIEECRRTCIG'
        process_models.extract_and_validate_models(amoptd)
        self.assertTrue(amoptd['quark_models'])
        self.assertEqual(len(amoptd['processed_models']), 10)
        shutil.rmtree(temp_dir)


    def test_extract_and_validate_models_quark_old(self):
        temp_dir = tempfile.mkdtemp(dir=os.getcwd())
        amoptd = {'models' : os.path.join(self.testfiles_dir, 'decoys.tar.gz'),
                  'models_dir' : temp_dir,
                  'work_dir' : '.'}
#         sequence = 'QPRRKLCILHRNPGRCYDKIPAFYYNQKKKQCERFDWSGCGGNSNRFKTIEECRRTCIG'
        process_models.extract_and_validate_models(amoptd)
        self.assertTrue(amoptd['quark_models'])
        self.assertEqual(len(amoptd['processed_models']), 200)
        shutil.rmtree(temp_dir)


if __name__ == "__main__":
    unittest.main(verbosity=2)
