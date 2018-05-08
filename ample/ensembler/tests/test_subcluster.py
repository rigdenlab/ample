
import glob
import os
import unittest
from ample import constants
from ample.ensembler import subcluster
from ample.util import ample_util
from ample.testing import test_funcs

class Test_1(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.thisd = os.path.abspath(os.path.dirname( __file__ ))
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share, 'testfiles')

    def test_radius_cctbx(self):
        # Test we can reproduce the original thresholds
        radius = 8
        clusterer = subcluster.CctbxClusterer()
        # Only select a few as is very slow
        pdb_list = [os.path.join(self.testfiles_dir,"models",pdb) for pdb in ['1_S_00000001.pdb',
                                                                              '1_S_00000002.pdb',
                                                                              '1_S_00000003.pdb',
                                                                              '1_S_00000004.pdb'] ]
        clusterer.generate_distance_matrix(pdb_list)
        cluster_files1 = [os.path.basename(x) for x in clusterer.cluster_by_radius(radius)]
        ref = ['1_S_00000002.pdb', '1_S_00000004.pdb']
        self.assertItemsEqual(ref, cluster_files1)

    @unittest.skipUnless(test_funcs.found_exe("gesamt" + ample_util.EXE_EXT), "gesamt exec missing")
    def test_gesamt_matrix_generic(self):
        # Test we can reproduce the original thresholds
        gesamt_exe = ample_util.find_exe("gesamt" + ample_util.EXE_EXT)
        clusterer = subcluster.GesamtClusterer(executable=gesamt_exe)
        pdb_list = glob.glob(os.path.join(self.testfiles_dir, "models",'*.pdb'))
        clusterer._generate_distance_matrix_generic(pdb_list, purge_all=True)
        # Test two files manually
        index1 = 2
        index2 = 25
        f1 = pdb_list[index1]
        f2 = pdb_list[index2]
        # Run gesamt to get the score between the two
        logfile = 'gesamt.log'
        ample_util.run_command([gesamt_exe, f1, f2], logfile=logfile)
        qscore = None
        with open(logfile) as f:
            for l in f.readlines():
                if l.startswith(' Q-score'):
                    qscore = float(l.split()[2])

        self.assertIsNotNone(qscore, "No q-score found")
        # read score matrix
        matrix = []
        with open(subcluster.SCORE_MATRIX_NAME) as f:
            for l in f.readlines():
                if not l.strip(): continue
                fields = l.split()
                matrix.append((int(fields[0]),int(fields[1]), float(fields[2])))
        # Make sure the score matches
        for l in matrix:
            if l[0] == index1 and l[1] == index2:
                # Gesamt log and out file formats have different precisions
                self.assertAlmostEqual(l[2], qscore, 3, "Q-scores differ: {0} - {1}".format(l[2], qscore))
        os.unlink(logfile)
        os.unlink(subcluster.SCORE_MATRIX_NAME)
        os.unlink(subcluster.FILE_LIST_NAME)
        return

    @unittest.skipUnless(test_funcs.found_exe("gesamt" + ample_util.EXE_EXT), "gesamt exec missing")
    def test_gesamt_radius(self):
        # Test we can reproduce the original thresholds
        gesamt_exe = ample_util.find_exe("gesamt" + ample_util.EXE_EXT)
        clusterer = subcluster.GesamtClusterer(executable=gesamt_exe)
        pdb_list = glob.glob(os.path.join(self.testfiles_dir, "models",'*.pdb'))

        radius = 4
        clusterer.generate_distance_matrix(pdb_list, purge=True)
        cluster_files1 = [os.path.basename(x) for x in clusterer.cluster_by_radius(radius)]
        ref = ['1_S_00000002.pdb', '1_S_00000004.pdb', '1_S_00000005.pdb', '2_S_00000001.pdb',
              '2_S_00000005.pdb', '3_S_00000003.pdb', '3_S_00000004.pdb', '3_S_00000006.pdb',
              '4_S_00000002.pdb', '4_S_00000005.pdb', '5_S_00000004.pdb', '5_S_00000005.pdb']
        self.assertItemsEqual(ref,cluster_files1)
        #clusterer.dump_pdb_matrix('foo')
        return

    def test_radius_lsqkab(self):
        # Test we can reproduce the original thresholds
        clusterer = subcluster.LsqkabClusterer()
        pdb_list = glob.glob(os.path.join(self.testfiles_dir, "models",'*.pdb'))
        clusterer.generate_distance_matrix(pdb_list)
        clusterer.dump_pdb_matrix('lsqkab.matrix')
        os.unlink('lsqkab.matrix')
        return

    @unittest.skipUnless(test_funcs.found_exe("maxcluster" + ample_util.EXE_EXT), "maxcluster exec missing")
    def test_radius_maxcluster(self):
        # Test we can reproduce the original thresholds
        maxcluster_exe = ample_util.find_exe('maxcluster' + ample_util.EXE_EXT)
        radius = 4
        clusterer = subcluster.MaxClusterer(maxcluster_exe)
        pdb_list = glob.glob(os.path.join(self.testfiles_dir, "models",'*.pdb'))
        clusterer.generate_distance_matrix(pdb_list)
        cluster_files1 = [os.path.basename(x) for x in clusterer.cluster_by_radius(radius)]
        ref=['4_S_00000003.pdb', '2_S_00000005.pdb', '2_S_00000001.pdb', '3_S_00000006.pdb',
             '5_S_00000005.pdb', '3_S_00000003.pdb', '1_S_00000004.pdb', '4_S_00000005.pdb',
             '3_S_00000004.pdb', '1_S_00000002.pdb', '5_S_00000004.pdb', '4_S_00000002.pdb', '1_S_00000005.pdb']
        self.assertItemsEqual(ref,cluster_files1)
        os.unlink('files.list')
        os.unlink('maxcluster.log')

    @unittest.skipUnless(test_funcs.found_exe("maxcluster" + ample_util.EXE_EXT), "maxcluster exec missing")
    def test_cluster_score(self):
        maxcluster_exe = ample_util.find_exe('maxcluster' + ample_util.EXE_EXT)
        radius = 4
        clusterer = subcluster.MaxClusterer(maxcluster_exe)
        pdb_list = glob.glob(os.path.join(self.testfiles_dir, "models",'*.pdb'))

        clusterer.generate_distance_matrix(pdb_list)
        clusterer.cluster_by_radius(radius)
        variance = clusterer.cluster_score
        ref = 4.748
        self.assertLessEqual(abs(ref-variance), 0.001, "Incorrect variance: {0} -> {1}".format(variance, ref))


@unittest.skipUnless(test_funcs.found_exe("fast_protein_cluster" + ample_util.EXE_EXT), "fast_protein_cluster exec missing")
class Test_2(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ))
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share, 'testfiles')
        cls.fpc_exe = ample_util.find_exe("fast_protein_cluster" + ample_util.EXE_EXT)

    def test_indices_fpc(self):
        # Test we can reproduce the original thresholds
        radius = 4
        clusterer = subcluster.FpcClusterer(self.fpc_exe)
        pdb_list = glob.glob(os.path.join(self.testfiles_dir, "models",'*.pdb'))
        clusterer.generate_distance_matrix(pdb_list)
        indices = clusterer._cluster_indices(radius)
        ref=[2, 4, 9, 10, 11, 14, 15, 18, 19, 21, 23, 25, 28]
        self.assertEqual(ref, indices[0])
        os.unlink('files.list')
        os.unlink('cluster_output.names')
        os.unlink('cluster_output.cluster.stats')
        os.unlink('cluster_output.clusters')
        os.unlink('fpc.matrix')
        os.unlink('fast_protein_cluster.log')

    def test_radius_fpc(self):
        # Test we can reproduce the original thresholds
        radius = 4
        clusterer = subcluster.FpcClusterer(self.fpc_exe)
        pdb_list = glob.glob(os.path.join(self.testfiles_dir, "models",'*.pdb'))
        clusterer.generate_distance_matrix(pdb_list)
        cluster_files1 = [os.path.basename(x) for x in clusterer.cluster_by_radius(radius)]
        ref=['4_S_00000003.pdb', '2_S_00000005.pdb', '2_S_00000001.pdb', '3_S_00000006.pdb',
             '5_S_00000005.pdb', '3_S_00000003.pdb', '1_S_00000004.pdb', '4_S_00000005.pdb',
             '3_S_00000004.pdb', '1_S_00000002.pdb', '5_S_00000004.pdb', '4_S_00000002.pdb', '1_S_00000005.pdb']
        self.assertItemsEqual(ref,cluster_files1)
        os.unlink('files.list')
        os.unlink('cluster_output.names')
        os.unlink('cluster_output.cluster.stats')
        os.unlink('cluster_output.clusters')
        os.unlink('fpc.matrix')
        os.unlink('fast_protein_cluster.log')

if __name__ == "__main__":
    unittest.main()
