
import glob
import os
import shutil
import unittest
from ample.python import fast_protein_cluster
from ample.util import ample_util

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = cls.thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -2 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        return

    def _find_exe(self):
        try: 
            fpc_exe = ample_util.find_exe("fast_protein_cluster")
        except:
            self.assertTrue(False, "Cannot find fast_protein_cluster executable in environment")
        return fpc_exe

    def test_fpc_kmeans_rmsd(self):
        
        fpc_exe = self._find_exe()        
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        mdir=os.path.join(self.testfiles_dir,"models")
        models=glob.glob(mdir+os.sep+"*.pdb")
        
        wdir='fpc_test'
        if not os.path.isdir(wdir): os.mkdir(wdir)
        fpc=fast_protein_cluster.FPC()
        num_clusters=3
        score_type='rmsd'
        cluster_method='kmeans'
        clusters,cluster_data=fpc.cluster(models=models,
                                          num_clusters=num_clusters,
                                          nproc=4,
                                          score_type=score_type,
                                          cluster_method=cluster_method,
                                          work_dir=wdir,
                                          fpc_exe=fpc_exe,
                                          benchmark=True
                                          )
        
        self.assertEqual(len(clusters),num_clusters)
        d=cluster_data[0]
        self.assertEqual(d['cluster_num_models'],17)
        self.assertEqual(d['cluster_method'],'kmeans_rmsd')
        self.assertEqual(os.path.basename(d['cluster_centroid']),'4_S_00000005.pdb')
        shutil.rmtree(wdir)

    def test_fpc_hierarch_tm(self):
        
        fpc_exe = self._find_exe()        
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        
        mdir=os.path.join(self.testfiles_dir,"models")
        models=glob.glob(mdir+os.sep+"*.pdb")
        
        wdir='fpc_test'
        if not os.path.isdir(wdir): os.mkdir(wdir)
        fpc=fast_protein_cluster.FPC()
        num_clusters=1
        score_type='tm'
        cluster_method='hcomplete'
        clusters,cluster_data=fpc.cluster(models=models,
                                          num_clusters=num_clusters,
                                          score_type=score_type,
                                          cluster_method=cluster_method,
                                          work_dir=wdir,
                                          fpc_exe=fpc_exe,
                                          nproc=4,
                                          benchmark=True
                                          )
        
        self.assertEqual(len(clusters),num_clusters)
        d=cluster_data[0]
        self.assertEqual(d['cluster_num_models'],16)
        self.assertEqual(d['cluster_method'],'hcomplete_tm')
        self.assertEqual(os.path.basename(d['cluster_centroid']),'5_S_00000005.pdb')
        shutil.rmtree(wdir)

if __name__ == "__main__":
    unittest.main()