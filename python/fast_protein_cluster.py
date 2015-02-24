import os
import unittest

class FPC(object):
    """
    Class 
    """

    def __init__(self,pfile):
        pass

    def cluster(self,
                models=None,
                num_clusters=None,
                nproc=1,
                score_type="rmsd",
                cluster_method="kmeans",
                work_dir=None,
                fpc_exe=None):
  
        
        if not os.path.isdir(work_dir): os.mkdir(work_dir)
        os.chdir(work_dir)
        
        if not len(models) or not all([os.path.isfile(m) for m in models]):
            raise RuntimeError,"Missing models: {0}".format(models)
        
        # Create list of files
        flist='files.list'
        with open(flist,'w') as f:
            for m in models:
                f.write("{0}\n".format(os.path.abspath(m)))
        
        if not os.path.isfile(fpc_exe):
            raise RuntimeError,"Cannot find fast_protein_cluster executable: {0}".format(fpc_exe)
        
        # Build up the command-line
        cmd=[fpc_exe]
        if score_type=="rmsd":
            cmd += ['--rmsd']
        elif score_type=="tm":
            cmd += ['--tmscore']
        else:
            raise RuntimeError,"Unrecognised score_type: {0}".format(score_type)
        
        if cluster_method=="kmeans":
            cmd += ['--cluster_kmeans']
        elif cluster_method=="hcomplete":
            cmd += ['--cluster_hcomplete']
        else:
            raise RuntimeError,"Unrecognised cluster_method: {0}".format(cluster_method)
        
        if nproc > 1:
            cmd += ['--nthreads',str(nproc)]
        
        # Always save the distance matrix
        cmd += ['--write_text_matrix','matrix.txt']
        
        # Finally the list of files
        cmd += ['-i',flist]
        
        logfile=os.path.abspath("fast_protein_cluster.log")
        retcode = ample_util.run_command(cmd, logfile=logfile)
        if retcode != 0:
            msg = "non-zero return code for fast_protein_cluster in cluster!\nCheck logfile:{0}".format(logfile)
            #logging.critical(msg)
            raise RuntimeError, msg
    
        cluster_list='cluster_output.clusters'
        cluster_stats='cluster_output.cluster.stats'
        if not os.path.isfile(cluster_list) or os.path.isfile(cluster_stats):
            raise RuntimeError,"Cannot find files: {0} and {1}".format(cluster_list,cluster_stats)
        
        # Check stats and get centroids
        csizes=[]
        ccentroids=[]
        with open(cluster_stats) as f:
            for line in f:
                if line.startswith("Cluster:"):
                    fields=line.split()
                    csizes.append(int(fields[4]))
                    centroids.append(fields[7])
                    assert int(fields[1]) == len(csizes)+1,"Error parsing {0}".format(cluster_stats)
        
        if len(csizes) != num_clusters:
            raise RuntimeError,"Found {0} clusters in {1} but was expecting {2}".format(len(csizes),cluster_stats,num_clusters)
        
        clusters=[[] for i in range(num_clusters)]
        # Read in the clusters
        with open(cluster_list) as f:
            for line in f:
                fields=line.split()
                model=fields[0]
                idxCluster=int(fields[1])
                clusters[idxCluster].append(model)
        
        # Check
        for i,cs in enumerate(csizes):
            if not cs == len(clusters[i]):
                raise RuntimeError,"Cluster {0} size {1} does not match stats size {2}".format(i,len(clusters[i]),cs)
        
        # Create the data
        clusters_data=[]
        for i,csize in enumerate(csizes):
            cluster_data={}
            cluster_data['cluster_num']=i+1
            cluster_data['cluster_centroid']=centroids[i]
            cluster_data['cluster_num_models']=csize
            cluster_data['cluster_method']="{0}_{1}".format(cluster_method,score_type)
            cluster_data['num_clusters']=num_clusters
            clusters_data.append(cluster_data)
      
        return clusters, clusters_data


 
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
        cls.ample_dir = os.sep.join( paths[ : -1 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        return

    def testParse1(self):
        """parse 2bhw"""
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        
        
        return

def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testParse1'))
    return suite
    
#
# Run unit tests
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(testSuite())



