
import glob
import os
import unittest
from ample.python import subcluster
from ample.util import ample_util

from ample.python.subcluster import FILE_LIST_NAME, SCORE_MATRIX_NAME

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -2 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"testing")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        cls.maxcluster_exe=ample_util.find_exe('maxcluster')

    def test_radius_cctbx(self):
        # Test we can reproduce the original thresholds
        radius = 4
        clusterer = subcluster.CctbxClusterer()
        pdb_list = glob.glob(os.path.join(self.testfiles_dir,"models",'*.pdb'))
        clusterer.generate_distance_matrix(pdb_list)
        cluster_files1 = [os.path.basename(x) for x in clusterer.cluster_by_radius(radius)]
        ref=['4_S_00000003.pdb', '2_S_00000005.pdb', '2_S_00000001.pdb', 
             '3_S_00000006.pdb', '5_S_00000005.pdb', '3_S_00000003.pdb', 
             '1_S_00000004.pdb', '4_S_00000005.pdb', '3_S_00000004.pdb', 
             '1_S_00000002.pdb', '5_S_00000004.pdb', '4_S_00000002.pdb', 
             '1_S_00000005.pdb']
        self.assertItemsEqual(ref,cluster_files1)
    
    def test_gesamt_matrix(self):
        # Test we can reproduce the original thresholds
        gesamt_exe = os.path.join(os.environ['CCP4'], "bin", "gesamt")
        clusterer = subcluster.GesamtClusterer(executable = gesamt_exe)
        pdb_list = glob.glob(os.path.join(self.testfiles_dir,"models",'*.pdb'))
        clusterer.generate_distance_matrix(pdb_list, purge_all=True)
        # Test two files manually
        index1=2
        index2=25
        f1 = pdb_list[index1]
        f2 = pdb_list[index2]
        # Run gesamt to get the score between the two
        logfile = 'gesamt.log' 
        ample_util.run_command([gesamt_exe, f1, f2], logfile=logfile)
        qscore=None
        with open(logfile) as f:
            for l in f.readlines():
                if l.startswith(' Q-score'):
                    qscore = float(l.split()[2])
        
        self.assertIsNotNone(qscore, "No q-score found")
        # read score matrix
        matrix = []
        with open(SCORE_MATRIX_NAME) as f:
            for l in f.readlines():
                if not l.strip(): continue
                fields = l.split()
                matrix.append((int(fields[0]),int(fields[1]),float(fields[2])))
        # Make sure the score matches
        for l in matrix:
            if l[0] == index1 and l[1] == index2:
                # Gesamt log and out file formats have different precisions
                self.assertAlmostEqual(l[2], qscore, 3, "Q-scores differ: {0} - {1}".format(l[2], qscore))
        os.unlink(logfile)
        os.unlink(SCORE_MATRIX_NAME)
        os.unlink(FILE_LIST_NAME)
    
    def test_radius_lsqkab(self):
        # Test we can reproduce the original thresholds
        radius = 4
        clusterer = subcluster.LsqkabClusterer()
        pdb_list = glob.glob(os.path.join(self.testfiles_dir,"models",'*.pdb'))
        clusterer.generate_distance_matrix(pdb_list)
        clusterer.dump_pdb_matrix('lsqkab.matrix')
        os.unlink('lsqkab.matrix')

    def test_radius_maxcluster(self):
        # Test we can reproduce the original thresholds
        radius = 4
        clusterer = subcluster.MaxClusterer( self.maxcluster_exe )
        pdb_list = glob.glob(os.path.join(self.testfiles_dir,"models",'*.pdb'))
        clusterer.generate_distance_matrix( pdb_list )
        cluster_files1 = [os.path.basename(x) for x in clusterer.cluster_by_radius( radius )]
        ref=['4_S_00000003.pdb', '2_S_00000005.pdb', '2_S_00000001.pdb', '3_S_00000006.pdb',
             '5_S_00000005.pdb', '3_S_00000003.pdb', '1_S_00000004.pdb', '4_S_00000005.pdb',
             '3_S_00000004.pdb', '1_S_00000002.pdb', '5_S_00000004.pdb', '4_S_00000002.pdb', '1_S_00000005.pdb']
        self.assertItemsEqual(ref,cluster_files1)
        os.unlink('maxcluster.log')
    
    def test_indices_fpc(self):
        # Test we can reproduce the original thresholds
        try: 
            fpc_exe = ample_util.find_exe("fast_protein_cluster")
        except:
            self.assertTrue(False, "Cannot find fast_protein_cluster executable in environment")
        
        radius = 4
        clusterer = subcluster.FpcClusterer( fpc_exe )
        pdb_list = glob.glob(os.path.join(self.testfiles_dir,"models",'*.pdb'))
        clusterer.generate_distance_matrix( pdb_list )
        indices=clusterer._cluster_indices(radius) 
        ref=[2, 3, 5, 6, 13, 14, 15, 16, 22, 23, 25, 26, 27]
        self.assertEqual(ref,indices)
        os.unlink('files.list')
        os.unlink('cluster_output.names')
        os.unlink('cluster_output.cluster.stats')
        os.unlink('cluster_output.clusters')
        os.unlink('fpc.matrix')
        os.unlink('fast_protein_cluster.log')
    
    def test_radius_fpc(self):
        # Test we can reproduce the original thresholds
        try: 
            fpc_exe = ample_util.find_exe("fast_protein_cluster")
        except:
            self.assertTrue(False, "Cannot find fast_protein_cluster executable in environment")
        radius = 4
        clusterer = subcluster.FpcClusterer( fpc_exe )
        pdb_list = glob.glob(os.path.join(self.testfiles_dir,"models",'*.pdb'))
        clusterer.generate_distance_matrix( pdb_list )
        cluster_files1 = [os.path.basename(x) for x in clusterer.cluster_by_radius( radius )]
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
