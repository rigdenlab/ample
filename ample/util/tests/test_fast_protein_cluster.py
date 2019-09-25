import glob
import os
import shutil
import unittest
from ample import constants
from ample.testing import test_funcs
from ample.util import ample_util
from ample.util import fast_protein_cluster


@unittest.skip("fast_protein_cluster deprecated")
@unittest.skipUnless(
    test_funcs.found_exe("fast_protein_cluster" + ample_util.EXE_EXT), "fast_protein_cluster exec missing"
)
class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        cls.ample_share = constants.SHARE_DIR
        cls.testfiles_dir = os.path.join(cls.ample_share, 'testfiles')
        cls.fpc_exe = ample_util.find_exe("fast_protein_cluster" + ample_util.EXE_EXT)

    def test_fpc_kmeans_rmsd(self):
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
        mdir = os.path.join(self.testfiles_dir, "models")
        models = glob.glob(mdir + os.sep + "*.pdb")

        wdir = 'fpc_test'
        if not os.path.isdir(wdir):
            os.mkdir(wdir)
        fpc = fast_protein_cluster.FPC()
        num_clusters = 3
        score_type = 'rmsd'
        cluster_method = 'kmeans'
        clusters, cluster_data = fpc.cluster(
            models=models,
            num_clusters=num_clusters,
            nproc=4,
            score_type=score_type,
            cluster_method=cluster_method,
            work_dir=wdir,
            fpc_exe=self.fpc_exe,
            benchmark=True,
        )

        self.assertEqual(len(clusters), num_clusters)
        d = cluster_data[0]
        self.assertEqual(d['cluster_num_models'], 16)
        self.assertEqual(d['cluster_method'], 'kmeans_rmsd')
        self.assertEqual(os.path.basename(d['cluster_centroid']), '4_S_00000005.pdb')
        shutil.rmtree(wdir)

    def test_fpc_hierarch_tm(self):
        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()

        mdir = os.path.join(self.testfiles_dir, "models")
        models = glob.glob(mdir + os.sep + "*.pdb")

        wdir = 'fpc_test'
        if not os.path.isdir(wdir):
            os.mkdir(wdir)
        fpc = fast_protein_cluster.FPC()
        num_clusters = 1
        score_type = 'tm'
        cluster_method = 'hcomplete'
        clusters, cluster_data = fpc.cluster(
            models=models,
            num_clusters=num_clusters,
            score_type=score_type,
            cluster_method=cluster_method,
            work_dir=wdir,
            fpc_exe=self.fpc_exe,
            nproc=4,
            benchmark=True,
        )

        self.assertEqual(len(clusters), num_clusters)
        d = cluster_data[0]
        self.assertEqual(d['cluster_num_models'], 16)
        self.assertEqual(d['cluster_method'], 'hcomplete_tm')
        self.assertEqual(os.path.basename(d['cluster_centroid']), '5_S_00000005.pdb')
        shutil.rmtree(wdir)


if __name__ == "__main__":
    unittest.main()
