#!/usr/bin/env ccp4-python

#edit the sidechains to make polyala, all and reliable

from collections import namedtuple
import itertools
import logging
import re
import os
import shutil

import mmtbx.superpose
import numpy

from ample.util import ample_util
from ample.util import pdb_edit

_logger = logging.getLogger()
logging.basicConfig(level=logging.DEBUG)

SCORE_MATRIX_NAME = 'score.matrix'
FILE_LIST_NAME = 'files.list'
RMSD_MAX = 50
QSCORE_MIN = 0.01

class SubClusterer(object):
    """Base class for clustering pdbs by distance
    Sub-classes just need to provide a generate_distance_matrix class
    """
    
    def __init__(self,executable=None):
        if executable and not os.path.exists(executable) and os.access(executable, os.X_OK):
            raise RuntimeError,"Cannot find subclusterer executable: {0}".format(executable) 
        self.executable = executable
        self.distance_matrix = None
        self.index2pdb = []
        self.cluster_score = None
        return
    
    def generate_distance_matrix(self,pdb_list):
        assert False

    def cluster_by_radius(self, radius):
        """Return a list of pdbs clustered by the given radius"""
        if self.distance_matrix is None:
            raise RuntimeError("Need to call generate_distance_matrix before cluster_by_radius!")
        cluster_indices, cluster_score = self._cluster_indices(radius)
        self.cluster_score = cluster_score
        if cluster_indices:
            return [ self.index2pdb[i] for i in cluster_indices ]
        else:
            return None

    def _cluster_indices(self,thresh):
        """Return the indices of the largest cluster that have distances < thresh.
        We loop through each row of the distance matrix and for each row (pdb) see
        how many pdbs are < thresh to this pdb. We return the largest cluster.
        """
        #self.dump_matrix("maxcluster.csv")
        thresh = float(thresh)
        
        # get mask of all elements where condition is true. We exclude 0.0 to ensure we don't get the
        # index of the model that the row is compared with, as this needs to be the first model in the 
        # ensemble. This means we would also exclude models that had an rmsd of zero to the centroid, but
        # as these are likely to be identical (and this occurrence rare), this should be ok
        condition = numpy.logical_and(self.distance_matrix <= thresh, self.distance_matrix != 0.0)
        # Array of sums of each row - largest number is a row where most items satisfy condition
        condition_sum =  sum(condition)
        # Find all rows that have the maximum of the condition true and then select the first one
        row_index = numpy.where(condition_sum == numpy.max(condition_sum))[0][0]
        # Select all values from that row where the condition is true and insert the first index so that
        # it becomes the centroid of that cluster 
        max_cluster = numpy.insert(numpy.where(condition[row_index])[0], 0, row_index)
        
#             max_cluster = []
#             len_matrix = len(self.distance_matrix)
#             for i in range(len_matrix):
#                 cluster = [i]
#                 for j in range(len_matrix):
#                     if self.distance_matrix[i][j] is None or j==i: continue
#                     if float(self.distance_matrix[i][j]) < thresh:
#                         cluster.append(j)
#                 if len(cluster) > len(max_cluster):
#                     max_cluster = copy.copy(cluster)

        if len(max_cluster) == 1:
            return None, None
        else:
            cluster_score = self.calculate_score(max_cluster)
            return sorted(max_cluster), cluster_score
    
    def calculate_score(self, cluster):
        """Given a list of indices of a cluster, calculate the rmsd we want to give to phaser 
        """
        ALL_BY_ALL = True
        if ALL_BY_ALL:
            # Calculate all the rmsds of all decoys in the cluster with each other
#             rmsds = []
#             lenc = len(cluster)
#             for i in range(lenc):
#                 for j in range(i+1,lenc):
#                     i1 = cluster[i]
#                     i2 = cluster[j]
#                     rmsds.append(self.distance_matrix[i1][i2])
            rmsds = [ self.distance_matrix[i] for i in itertools.combinations(cluster,2) ]
        else:
            # Just use the rmsds of the decoys to the the cluster centroid - assumes
            # the centroid approximates the native
            row = cluster[0]
            rmsds = [self.distance_matrix[row][j] for j in cluster[1:]]

        #print "rmsds: ",rmsds
        #print "mean: ",statistics_util.mean(rmsds)
        #print "max: ",max(rmsds)
        return max(rmsds)
        #return statistics_util.mean(rmsds)
    
    def dump_raw_matrix(self,file_name):
        with open(file_name,'w') as f:
            for row in self.distance_matrix:
                f.write(",".join(map(str,row))+"\n")
            f.write("\n")
        return
            
    def dump_pdb_matrix(self, file_name=SCORE_MATRIX_NAME, offset=0):
        with open(file_name,'w') as f:
            l = len(self.distance_matrix) + offset
            for i in range(offset, l):
                for j in range(i, l):
                    f.write("{0: > 4d} {1: > 4d} {2: > 8.3F}\n".format(i,j, self.distance_matrix[i-offset][j-offset]))
            f.write("\n")
        return os.path.abspath(file_name)

class CctbxClusterer(SubClusterer):
    """Class to cluster files with maxcluster"""

    def generate_distance_matrix(self, pdb_list):
        """Run cctbx to generate the distance distance_matrix"""
         
        num_models = len(pdb_list)
        if not num_models:
            msg = "generate_distance_matrix got empty pdb_list!"
            logging.critical(msg)
            raise RuntimeError, msg
 
        # Index is just the order of the pdb in the file
        self.index2pdb=pdb_list
     
        # Create a square distance_matrix num_models in size filled with None
        self.distance_matrix = numpy.zeros([num_models, num_models])
        # Set zeros diagonal
         
        for i, m1 in enumerate(pdb_list):
            fixed = mmtbx.superpose.SuperposePDB(m1, preset='ca', log=None, quiet=True)
            for j, m2 in enumerate(pdb_list):
                if j <= i: continue
                moving = mmtbx.superpose.SuperposePDB(m2, preset='ca', log=None, quiet=True)
                rmsd, lsq = moving.superpose(fixed)
                self.distance_matrix[i][j]=float(rmsd)
         
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
        with open( fname, 'w' ) as f: f.write( "\n".join( pdb_list )+"\n" )
            
        # Index is just the order of the pdb in the file
        self.index2pdb = pdb_list
        
        # Run fast_protein_cluster - this is just to generate the distance matrix, but there
        # doesn't seem to be a way to stop it clustering as well - not a problem as it just
        # generates more files
        log_name = os.path.abspath("fast_protein_cluster.log")
        matrix_file = "fpc.matrix"
        cmd = [self.executable,
               "--cluster_write_text_matrix",
               matrix_file,
               "-i",
               fname]
               
        retcode = ample_util.run_command( cmd, logfile=log_name )
        if retcode != 0:
            msg = "non-zero return code for fast_protein_cluster in generate_distance_matrix!\nCheck logfile:{0}".format(log_name)
            logging.critical( msg )
            raise RuntimeError, msg

        mlen=0
        data=[]
        with open(matrix_file) as f:
            for l in f:
                l = l.strip().split()
                x = int(l[0])
                y = int(l[1])
                d = float(l[2])
                mlen = max(mlen,x+1) # +1 as we want the length
                data.append((x,y,d))
         
        # create empty matrix - we use None's but this means we need to check for then when
        # looking through the matrix
        # use square matrix to make indexing easier as we're unlikely to be very big
        m = numpy.zeros([mlen, mlen])
         
        # Fill in all values (upper triangle)
        for i,j,d in data:
            if i > j:
                m[j][i] = d
            else:
                m[i][j] = d
                 
        # Copy to lower
        for x in range(mlen):
            for y in range(mlen):
                if x==y: continue
                m[y][x] = m[x][y]
                
        self.distance_matrix = m
        return
    
class GesamtClusterer(SubClusterer):
    """Class to cluster files with Gesamt"""
    
    def generate_distance_matrix(self, pdb_list, nproc=1, purge=True):
        if True:
            self._generate_pairwise_rmsd_matrix(pdb_list, nproc=nproc, purge=purge)
        else:
            self._generate_distance_matrix_generic(self,
                                                   pdb_list,
                                                   purge=purge,
                                                   purge_all=False,
                                                   metric='qscore',
                                                   nproc=nproc)
        return
    
    def _generate_pairwise_rmsd_matrix(self, models, nproc=1, purge=False):
        """
        Use gesamt to generate an all-by-all pairwise rmsd matrix of a list of pdb models
        
        Notes: 
        gesamt -input-list inp_list.dat -sheaf-x

where inp_list.dat  contains:

1ADZ.pdb -s /1/A
1ADZ.pdb -s /2/A
1ADZ.pdb -s /3/A

        """

        # Index is just the order of the pdb in the file
        self.index2pdb = models

        # Create file with list of pdbs and model/chain
        glist = 'gesamt_models.dat'
        with open(glist, 'w') as w:
            for m in models:
                w.write("{0} -s /1/A \n".format(m))
            w.write('\n')
        
        cmd = [ self.executable, '-input-list', glist, '-sheaf-x', '-nthreads={0}'.format(nproc)]

        logfile = os.path.abspath('gesamt_archive.log')
        rtn = ample_util.run_command(cmd, logfile)
        if rtn != 0: raise RuntimeError("Error running gesamt - check logfile: {0}".format(logfile))
        
        # Create a square distance_matrix no_models in size filled with None
        num_models = len(models)
        self.distance_matrix = numpy.zeros([num_models, num_models])
        
        # Read in the rmsds calculated
        self._parse_gesamt_rmsd_log(logfile, num_models)
        
        if purge:
            os.unlink(glist)
            os.unlink(logfile)
        return
    
    def _parse_gesamt_rmsd_log(self, logfile, num_models):
        found = False
        with open(logfile) as f:
            while True:
                line = f.readline()
                if line.startswith(' ===== CROSS-RMSDs'):
                    line = f.readline() # skip blank line
                    for i in range(num_models):
                        line = f.readline().strip()
                        fields = line.split('|')
                        
                        # paranoid check
                        nmodel = int(fields[0])
                        assert nmodel ==  i+1,"Error parsing gesamt logfile: {0}".format(logfile)
                        
                        rmsd_txt = fields[2].strip()
                        # poke into distance matrix
                        rmsds = [ float(r) for r in rmsd_txt.split() ]
                        for j in range(len(rmsds)):
                            if j == i: continue
                            self.distance_matrix[i][j] = rmsds[j]
                    found = True
                    break
        
        if not found: raise RuntimeError("Could not generate distance matrix with gesamt")
        return
        
    def _generate_distance_matrix_generic(self, models, purge=True, purge_all=False, metric='qscore', nproc=1):
        # Make sure all the files are in the same directory otherwise we wont' work
        mdir = os.path.dirname(models[0])
        if not all([ os.path.dirname(p) == mdir for p in models ]):
            raise RuntimeError("All pdb files are not in the same directory!")
        
        # Create list of pdb files
        fname = os.path.join(os.getcwd(), FILE_LIST_NAME)
        with open(fname, 'w') as f: f.write("\n".join(models)+"\n")
            
        # Index is just the order of the pdb in the file
        self.index2pdb = models
        nmodels = len(models)
        
        # Make the archive
        _logger.debug("Generating gesamt archive from models in directory {0}".format(mdir))
        garchive = 'gesamt.archive'
        if not os.path.isdir(garchive): os.mkdir(garchive)
        logfile = os.path.abspath('gesamt_archive.log')
        cmd = [ self.executable, '--make-archive', garchive, '-pdb', mdir ]
        #cmd += [ '-nthreads=auto' ]
        cmd += [ '-nthreads={0}'.format(nproc) ]
        # HACK FOR DYLD!!!!
        env = None
        #env = {'DYLD_LIBRARY_PATH' : '/opt/ccp4-devtools/install/lib'}
        rtn = ample_util.run_command(cmd, logfile,env = env )
        if rtn != 0: raise RuntimeError("Error running gesamt - check logfile: {0}".format(logfile))
        
        if purge_all: os.unlink(logfile)
        
        # Now loop through each file creating the matrix
        if metric == 'rmsd':
            parity = 0.0
        elif metric == 'qscore':
            parity = 1
        else: raise RuntimeError("Unrecognised metric: {0}".format(metric))
        
        #m = [[parity for _ in range(nmodels)] for _ in range(nmodels)]
        m = numpy.full([nmodels, nmodels], parity)
        
        for i, model in enumerate(models):
            mname = os.path.basename(model)
            gesamt_out = '{0}_gesamt.out'.format(mname)
            logfile = '{0}_gesamt.log'.format(mname)
            cmd = [ self.executable, model, '-archive', garchive, '-o', gesamt_out ]
            cmd += [ '-nthreads={0}'.format(nproc) ]
            rtn = ample_util.run_command(cmd, logfile)
            if rtn != 0: raise RuntimeError("Error running gesamt!")
            else:
                if purge: os.unlink(logfile)
                
            gdata = self._parse_gesamt_out(gesamt_out)
            assert gdata[0].file_name == mname, gdata[0].file_name + " " + mname
            score_dict = { g.file_name: (g.rmsd, g.q_score) for g in gdata  }
            
            for j in range(i + 1, nmodels):
                # Try and get the rmsd and qscore for this model. If it's missing we assume the model was 
                # too divergent for gesamt to find it and we set the rmsd and qscore to fixed values
                model2 = os.path.basename(models[j])
                try:
                    rmsd, qscore = score_dict[model2]
                except KeyError:
                    rmsd = RMSD_MAX
                    qscore = QSCORE_MIN
                
                if metric == 'rmsd':
                    score = rmsd
                elif metric == 'qscore':
                    score = qscore
                else: raise RuntimeError("Unrecognised metric: {0}".format(metric))
                
                m[i][j] = score
                
            if purge_all: os.unlink(gesamt_out)
                    
        # Copy upper half of matrix to lower
        for x in range(nmodels):
            for y in range(nmodels):
                if x == y: continue
                m[y][x] = m[x][y]
                
        self.distance_matrix = m

        # Remove the gesamt archive
        if purge: shutil.rmtree(garchive)
        
        # Write out the matrix in a form spicker can use
        self.dump_pdb_matrix(SCORE_MATRIX_NAME)
        return

    def _parse_gesamt_out(self, out_file):
        # Assumption is there are no pdb_codes
        GesamtData = namedtuple('GesamtData', ['count', 'chain_id', 'q_score', 'rmsd', 'seq_id', 'nalign', 'nres', 'file_name'])
        data = []
        with open(out_file) as f:
            for i, line in enumerate(f):
                if i < 2: continue # First 2 lines are headers
                if not line.strip(): continue # ignore blanks
                try:
                    tmp = GesamtData(*line.split())
                    # Convert from strings to correct types
                    data.append(GesamtData(int(tmp.count),
                                           tmp.chain_id,
                                           float(tmp.q_score),
                                           float(tmp.rmsd),
                                           tmp.seq_id,
                                           int(tmp.nalign),
                                           int(tmp.nres),
                                           os.path.basename(tmp.file_name)))
                except Exception as e:
                    msg = 'Error parsing line {0}: {1}\n{2}'.format(i, line, e.message)
                    logging.critical(msg)
                    raise e  
                
        assert len(data),"Failed to read any data!"
        return data
    
class LsqkabClusterer(SubClusterer):
    """Class to cluster files with Lsqkab"""
    
    def calc_rmsd(self, model1, model2, nresidues=None, logfile='lsqkab.out', purge=False):
        
        if not nresidues:  _, nresidues = pdb_edit.num_atoms_and_residues(model1, first=True)
        
        stdin = """FIT RESIDUE CA 1 TO {0} CHAIN {1}
MATCH 1 to  {0} CHAIN {1}
output  RMS    
end""".format(nresidues, 'A')

        cmd = [ 'lsqkab', 'XYZINM', model1, 'XYZINF', model2 ]
        
        ample_util.run_command(cmd, logfile=logfile, stdin=stdin)
        rmsd =  self.parse_lsqkab_output(logfile)
        
        # cleanup 
        if purge:
            os.unlink(logfile)
            os.unlink('RMSTAB')
        
        return rmsd
                
    def generate_distance_matrix(self, models, purge=True, metric='qscore'):

        # Index is just the order of the pdb in the file
        self.index2pdb = models
        num_models = len(models)
        
        # Assume all models are the same size and only have a single chain
        # We also assume that the chain is called 'A' (not relevant here)
        _, nresidues = pdb_edit.num_atoms_and_residues(models[0], first=True)
        
        # Create a square distance_matrix no_models in size filled with None
        self.distance_matrix = numpy.zeros([num_models, num_models])
        
        logfile='lsqkab.out'
        parity = 0.0
        for i, fixed in enumerate(models):
            for j, model2 in enumerate(models):
                if j < i: continue
                if j == i:
                    rmsd = parity
                elif j > i:
                    rmsd = self.calc_rmsd(fixed, model2, nresidues=nresidues, logfile=logfile)
                self.distance_matrix[i][j] = rmsd
        
        # Clean up output files from lsqkab
        os.unlink(logfile)
        os.unlink('RMSTAB')
         
        # Copy in other half of matrix - we use a full matrix as it's easier to scan for clusters
        for x in range(len(self.distance_matrix)):
            for y in range(len(self.distance_matrix)):
                self.distance_matrix[y][x] = self.distance_matrix[x][y]
        return
    
    def parse_lsqkab_output(self, output_file):
        with open(output_file) as f:
            for  l in f.readlines():
                if l.startswith("          RMS     XYZ DISPLACEMENT ="):
                    return float(l.split()[4])
        assert False

class MaxClusterer(SubClusterer):
    """Class to cluster files with maxcluster"""
    
    def generate_distance_matrix(self, pdb_list):
        """Run maxcluster to generate the distance distance_matrix"""
        
        num_models = len( pdb_list )
        if not num_models:
            msg = "generate_distance_matrix got empty pdb_list!"
            logging.critical( msg )
            raise RuntimeError, msg
        
        self.index2pdb=[0]*num_models
    
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
        fname = os.path.join(os.getcwd(), FILE_LIST_NAME )
        with open( fname, 'w' ) as f:
            f.write( "\n".join( pdb_list )+"\n" )
            
        #log_name = "maxcluster_radius_{0}.log".format(radius)
        log_name = os.path.abspath("maxcluster.log")
        cmd = [ self.executable, "-l", fname, "-L", "4", "-rmsd", "-d", "1000", "-bb", "-C0" ]
        retcode = ample_util.run_command( cmd, logfile=log_name )
        
        if retcode != 0:
            msg = "non-zero return code for maxcluster in generate_distance_matrix!\nSee logfile: {0}".format(log_name)
            logging.critical( msg )
            raise RuntimeError, msg
        
        # Create a square distance_matrix no_models in size filled with None
        parity = 0.0
        self.distance_matrix = numpy.full([num_models, num_models], parity)

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
                self.distance_matrix[  int(split[1]) -1 ][  int(split[3]) -1]  = float(split[5])
    
                if split[2]+'.pdb' not  in self.index2pdb:
                    self.index2pdb[int(split[1]) -1]  =  split[2]+'.pdb'
    
                if split[4]+'.pdb' not  in self.index2pdb:
                    self.index2pdb[int(split[3]) -1]  =  split[4]+'.pdb'
    
        # Copy in other half of matrix - we use a full matrix as it's easier to scan for clusters
        for x in range(len(self.distance_matrix)):
            for y in range(len(self.distance_matrix)):
                self.distance_matrix[y][x] = self.distance_matrix[x][y]
        return
