#!/usr/bin/env ccp4-python

#edit the sidechains to make polyala, all and reliable

import ample_util
import copy
import glob
import logging
import re
import os
import unittest


class SubClusterer(object):
    """Base class for clustering pdbs by distance
    Sub-classes just need to provide a generate_distance_matrix class
    """
    
    def __init__(self,executable):
        if not os.path.exists(executable) and os.access(executable, os.X_OK):
            raise RuntimeError,"Cannot find subclusterer executable: {0}".format(executable) 
        self.executable = executable
        self.distance_matrix = None
        self.index2pdb = []
        return
    
    def generate_distance_matrix(self,pdb_list):
        assert False

    def cluster_by_radius(self, radius):
        """Return a list of pdbs clustered by the given radius"""
        if self.distance_matrix is None:
            raise RuntimeError,"Need to call generate_distance_matrix before cluster_by_radius!"
        return [ self.index2pdb[i] for i in self._cluster_indices(radius)]

    def _cluster_indices(self,thresh):
        """Return the indices of the largest cluster that have distances < thresh.
        We loop through each row of the distance matrix and for each row (pdb) see
        how many pdbs are < thresh to this pdb. We return the largest cluster.
        """
        #self.dump_matrix("maxcluster.csv")
        thresh=float(thresh)
        max_cluster=[]
        m=self.distance_matrix
        for i in range(len(m)):
            cluster=[i]
            for j in range(len(m)):
                if m[i][j] is None or j==i: continue
                if float(m[i][j]) < thresh:
                    cluster.append(j)
            if len(cluster) > len(max_cluster):
                max_cluster=copy.copy(cluster)
        return sorted(max_cluster)
    
    def dump_matrix(self,file_name):
        with open(file_name,'w') as f:
            for row in self.distance_matrix:
                f.write(",".join(map(str,row))+"\n")
            f.write("\n")
        
class MaxClusterer(SubClusterer):
    """Class to cluster files with maxcluster"""
    
    def generate_distance_matrix(self, pdb_list):
        """Run maxcluster to generate the distance distance_matrix"""
        
        no_models = len( pdb_list )
        if not no_models:
            msg = "generate_distance_matrix got empty pdb_list!"
            logging.critical( msg )
            raise RuntimeError, msg
        
        self.index2pdb=[0]*no_models
    
        # Maxcluster arguments
        # -l [file]   File containing a list of PDB model fragments
        # -L [n]      Log level (default is 4 for single MaxSub, 1 for lists)
        # -d [f]      The distance cut-off for search (default auto-calibrate)
        # -bb         Perform RMSD fit using backbone atoms
        #     -C [n]      Cluster method: 0 - No clustering
        # -rmsd ???
        #os.system(MAX + ' -l list  -L 4 -rmsd -d 1000 -bb -C0 >MAX_LOG ')
        #print 'MAX Done'
        
        # Create the list of files for maxcluster
        fname = os.path.join(os.getcwd(), "files.list" )
        with open( fname, 'w' ) as f:
            f.write( "\n".join( pdb_list )+"\n" )
            
        #log_name = "maxcluster_radius_{0}.log".format(radius)
        log_name = "maxcluster.log"
        cmd = [ self.executable, "-l", fname, "-L", "4", "-rmsd", "-d", "1000", "-bb", "-C0" ]
        retcode = ample_util.run_command( cmd, logfile=log_name )
        
        if retcode != 0:
            msg = "non-zero return code for maxcluster in generate_distance_matrix!"
            logging.critical( msg )
            raise RuntimeError, msg
        
        # Create a square distance_matrix no_models in size filled with None
        self.distance_matrix = [[None for col in range(no_models)] for row in range(no_models)]
    
        #jmht Save output for parsing - might make more sense to use one of the dedicated maxcluster output formats
        #max_log = open(cur_dir+'/MAX_LOG')
        max_log = open( log_name, 'r')
        pattern = re.compile('INFO  \: Model')
        for line in max_log:
            if re.match(pattern, line):
    
                # Split so that we get a list with
                # 0: model 1 index
                # 1: path to model 1 without .pdb suffix
                # 2: model 2 index
                # 3: path to model 2 without .pdb suffix
                # 4: distance metric
                split = re.split('INFO  \: Model\s*(\d*)\s*(.*)\.pdb\s*vs\. Model\s*(\d*)\s*(.*)\.pdb\s*=\s*(\d*\.\d*)', line)
                self.distance_matrix[  int(split[1]) -1 ][  int(split[3]) -1]  = split[5]
    
                if split[2]+'.pdb' not  in self.index2pdb:
                    self.index2pdb[int(split[1]) -1]  =  split[2]+'.pdb'
    
                if split[4]+'.pdb' not  in self.index2pdb:
                    self.index2pdb[int(split[3]) -1]  =  split[4]+'.pdb'
    
        # Copy in other half of matrix - we use a full matrix as it's easier to scan for clusters
        for x in range(len(self.distance_matrix)):
            for y in range(len(self.distance_matrix)):
                self.distance_matrix[y][x] = self.distance_matrix[x][y]
        return
    
class FpcClusterer(SubClusterer):
    """Class to cluster files with fast_protein_clusterer"""
    
    def generate_distance_matrix(self,pdb_list):
        
        # Create list of pdb files
        fname = os.path.join(os.getcwd(), "files.list" )
        with open( fname, 'w' ) as f:
            f.write( "\n".join( pdb_list )+"\n" )
            
        # Index is just the order of the pdb in the file
        self.index2pdb=pdb_list
        
        # Run fast_protein_cluster - this is just to generate the distance matrix, but there
        # doesn't seem to be a way to stop it clustering as well - not a problem as it just
        # generates more files
        log_name = "fast_protein_cluster.log"
        matrix_file = "fpc.matrix"
        cmd = [self.executable,
               "--cluster_write_text_matrix",
               matrix_file,
               "-i",
               fname]
               
        retcode = ample_util.run_command( cmd, logfile=log_name )
        if retcode != 0:
            msg = "non-zero return code for fast_protein_cluster in generate_distance_matrix!"
            logging.critical( msg )
            raise RuntimeError, msg

        mlen=0
        data=[]
        with open(matrix_file) as f:
            for l in f:
                l=l.strip().split()
                x=int(l[0])
                y=int(l[1])
                d=float(l[2])
                mlen=max(mlen,x+1) # +1 as we want the length
                data.append((x,y,d))
         
        # create empty matrix - we use None's but this means we need to check for then when
        # looking through the matrix
        # use square matrix to make indexing easier as we're unlikely to be very big
        m=[[None for i in range(mlen)] for j in range(mlen)]
         
        # Fill in all values (upper triangle)
        for i,j,d in data:
            if i > j:
                m[j][i]=d
            else:
                m[i][j]=d
                 
        # Copy to lower
        for x in range(mlen):
            for y in range(mlen):
                if x==y: continue
                m[y][x] = m[x][y]
                
        self.distance_matrix=m
        
        return

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
        cls.ample_dir = os.sep.join( paths[ : -1 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        cls.maxcluster_exe=ample_util.find_exe('maxcluster')
        cls.fpc_exe="/opt/fast_protein_cluster.1.1.2/fast_protein_cluster"
        return

    def testIndicesMaxcluster(self):
        """Test we can reproduce the original thresholds"""

        radius = 4
        clusterer = MaxClusterer( self.maxcluster_exe )
        pdb_list = glob.glob(os.path.join(self.ample_dir,'examples','toxd-example','models','*.pdb'))
        clusterer.generate_distance_matrix( pdb_list )
        indices=clusterer._cluster_indices(radius) 

        ref=[2, 3, 5, 6, 13, 14, 15, 16, 22, 23, 25, 26, 27]
        self.assertEqual(ref,indices)
        os.unlink('files.list')
        os.unlink('maxcluster.log')
        return
    
    def testRadiusMaxcluster(self):
        """Test we can reproduce the original thresholds"""

        radius = 4
        clusterer = MaxClusterer( self.maxcluster_exe )
        pdb_list = glob.glob(os.path.join(self.ample_dir,'examples','toxd-example','models','*.pdb'))
        clusterer.generate_distance_matrix( pdb_list )
        cluster_files1 = [os.path.basename(x) for x in clusterer.cluster_by_radius( radius )]
        
        ref=['4_S_00000003.pdb', '2_S_00000005.pdb', '2_S_00000001.pdb', '3_S_00000006.pdb',
             '5_S_00000005.pdb', '3_S_00000003.pdb', '1_S_00000004.pdb', '4_S_00000005.pdb',
             '3_S_00000004.pdb', '1_S_00000002.pdb', '5_S_00000004.pdb', '4_S_00000002.pdb', '1_S_00000005.pdb']
        
        self.assertEqual(ref,cluster_files1)
        
        os.unlink('maxcluster.log')

        return
    
    def testIndicesFpc(self):
        """Test we can reproduce the original thresholds"""
    
        radius = 4
        clusterer = FpcClusterer( self.fpc_exe )
        pdb_list = glob.glob(os.path.join(self.ample_dir,'examples','toxd-example','models','*.pdb'))
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
        return
    
    def testRadiusFpc(self):
        """Test we can reproduce the original thresholds"""
    
        radius = 4
        clusterer = FpcClusterer( self.fpc_exe )
        pdb_list = glob.glob(os.path.join(self.ample_dir,'examples','toxd-example','models','*.pdb'))
        clusterer.generate_distance_matrix( pdb_list )

        cluster_files1 = [os.path.basename(x) for x in clusterer.cluster_by_radius( radius )]
        
        ref=['4_S_00000003.pdb', '2_S_00000005.pdb', '2_S_00000001.pdb', '3_S_00000006.pdb',
             '5_S_00000005.pdb', '3_S_00000003.pdb', '1_S_00000004.pdb', '4_S_00000005.pdb',
             '3_S_00000004.pdb', '1_S_00000002.pdb', '5_S_00000004.pdb', '4_S_00000002.pdb', '1_S_00000005.pdb']
        
        self.assertEqual(ref,cluster_files1)

        os.unlink('files.list')
        os.unlink('cluster_output.names')
        os.unlink('cluster_output.cluster.stats')
        os.unlink('cluster_output.clusters')
        os.unlink('fpc.matrix')
        os.unlink('fast_protein_cluster.log')
        return

def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testRadiusMaxcluster'))
    suite.addTest(Test('testIndicesMaxcluster'))
    suite.addTest(Test('testIndicesFpc'))
    suite.addTest(Test('testRadiusFpc'))
    return suite
    
#
# Run unit tests
if __name__ == "__main__":
 
    unittest.TextTestRunner(verbosity=2).run(testSuite())



