#!/usr/bin/env python

# python imports
import glob
import logging
import os
import re
import shutil
import sys

# our imports
from ample.util import ample_util
from ample.ensembler._ensembler import Cluster

logger = logging.getLogger(__name__)


class Spickerer(object):
    def __init__(self, spicker_exe=None, run_dir=None):
        """Initialise from a dictionary of options"""

        if not spicker_exe:
            if 'CCP4' in os.environ: spicker_exe = os.path.join(os.environ['CCP4'], 'bin', 'spicker')
            else: raise RuntimeError, "Cannot find a spicker executable!"
        if not (os.path.exists(spicker_exe) and os.access(spicker_exe, os.X_OK)):
            raise RuntimeError, "Cannot find a valid spicker executable: {0}".format(spicker_exe)
        self.spicker_exe = spicker_exe
        self.run_dir = run_dir
        self.results = None
        self.cluster_method = 'spicker'
        self.score_type = 'rmsd'
        return

    def get_length(self, pdb):
        pdb = open(pdb)
        counter = 0
        for line in pdb:
            pdb_pattern = re.compile('^ATOM')
            pdb_result = pdb_pattern.match(line)
            if pdb_result:
                atom = line[13:16]
                if re.search('CA', atom):
                    counter += 1
        return str(counter)

    def create_input_files(self, models, score_type='rmsd', score_matrix=None):
        """
        jmht
        Create the input files required to run spicker
        (See notes in spicker.f FORTRAN file for a description of the required files)
        """
        if not len(models):
            msg = "no models provided!"
            logger.critical(msg)
            raise RuntimeError, msg

        if score_type: score_type = score_type.lower()
        logger.debug("Using score_type: {0}".format(score_type))

        if score_type == 'read_matrix':
            if not (score_matrix and os.path.isfile(score_matrix)):
                msg = 'Cannot find score_matrix: {0}'.format(score_matrix)
                logger.critical(msg)
                raise RuntimeError, msg
            logger.debug("Using score_matrix: {0}".format(score_matrix))
            shutil.copy(score_matrix, os.path.join(self.run_dir, 'score.matrix'))


#         elif score_type == 'tm':
#             # Create file so spicker knows to calculate TM scores
#             with open('TM.score','w') as f: f.write('\n')
#         elif score_type == 'rmsd':
#             pass
#         else:
#             raise RuntimeError,"Unknown score_type: {0}".format(score_type)

# read_out - Input file for spicker with coordinates of the CA atoms for each of the PDB structures
#
# file_list - a list of the full path of all PDBs - used so we can loop through it and copy the selected
# ones to the relevant directory after we have run spicker - the order of these must match the order
# of the structures in the rep1.tra1 file
        list_string = ''
        counter = 0
        with open('rep1.tra1', "w") as read_out, open('file_list', "w") as file_list:
            for infile in models:
                pdbname = os.path.basename(infile)
                file_list.write(infile + '\n')
                list_string = list_string + pdbname + '\n'
                counter += 1
                length = self.get_length(infile)
                # 1st field is length, 2nd energy, 3rd & 4th don't seem to be used for anything
                read_out.write('\t' + length + '\t926.917       ' + str(counter) + '       ' + str(counter) + '\n')
                with open(infile) as read:
                    # Write out the coordinates of the CA atoms
                    for line in read:
                        pattern = re.compile(
                            '^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)'
                        )
                        result = re.match(pattern, line)
                        if result:
                            split = re.split(pattern, line)
                            if split[2] == 'CA':
                                read_out.write('     ' + split[6] + '     ' + split[7] + '     ' + split[8] + '\n')

        # from spicker.f
        # *       'rmsinp'---Mandatory, length of protein & piece for RMSD calculation;
        with open('rmsinp', "w") as rmsinp:
            rmsinp.write('1  ' + length + '\n\n')
            rmsinp.write(length + '\n')

        # make tra.in
        # from spicker.f
        # *       'tra.in'---Mandatory, list of trajectory names used for clustering.
        # *                  In the first line of 'tra.in', there are 3 parameters:
        # *                  par1: number of decoy files
        # *                  par2: 1, default cutoff, best for decoys from template-based
        # *                           modeling;
        # *                       -1, cutoff based on variation, best for decoys from
        # *                           ab initio modeling.
        # *                       -2, use TM scores for clustering
        # *                  par3: 1, closc from all decoys; -1, closc clustered decoys
        # *                  From second lines are the file names which contain coordinates
        # *                  of 3D structure decoys. All these files are mandatory
        par2 = '-2' if score_type == 'tm' else '-1'
        with open('tra.in', "w") as tra:
            tra.write('1 {0} 1 \nrep1.tra1\n'.format(par2))

        # Create the file with the sequence of the PDB structures
        # from spicker.f
        # *       'seq.dat'--Mandatory, sequence file, for output of PDB models.
        with open('seq.dat', "w") as seq, open(models[0], 'r') as a_pdb:
            for line in a_pdb:
                pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s*(\d*)\s')
                result = re.match(pattern, line)
                if result:
                    split = re.split(pattern, line)
                    if split[2] == 'CA':
                        seq.write('\t' + split[5] + '\t' + split[3] + '\n')
        return

    def cluster(self,
                models,
                num_clusters=10,
                max_cluster_size=200,
                run_dir=None,
                score_type='rmsd',
                score_matrix=None,
                nproc=1):
        """Cluster decoys using spicker
    
        Parameters
        ----------
        models : list
           A list containing structure decoys
        cluster_dir : str
           The directory to store the cluster data in
        cluster_method_type : str
           The method to be used to cluster the decoys
        score_type : str
           The scoring metric for clustering
        num_clusters : int
           The number of clusters to produce
        max_cluster_size : int
           The maximum number of decoys per cluster
        cluster_exe : str
           The path to the spicker executable
        nproc : int
           The number of processors to use
        score_matrix : str, optional
           The path to the score matrix to be used
    
        Returns
        -------
        list
           A list containing the clusters
    
        Raises
        ------
        RuntimeError
           No clusters returned by SPICKER
        """
        self._cluster(models, run_dir=run_dir, score_type=score_type, score_matrix=score_matrix, nproc=nproc)

        ns_clusters = len(self.results)
        if ns_clusters == 0: raise RuntimeError('No clusters returned by SPICKER')
        if ns_clusters < int(num_clusters):
            logger.critical('Requested {0} clusters but SPICKER only found {1} so using {1} clusters'.format(
                num_clusters, ns_clusters))
            num_clusters = ns_clusters

        clusters = []
        for i in range(num_clusters):
            cluster = self.results[i]
            # We truncate the list of models to max_cluster_size. This probably
            # needs to be redone, because as the models are ordered by their similarity
            # to the cluster centroid, we automatically select the 200 most similar
            # to the centroid. However if the cluster is large and the models similar
            # then when theseus calculates the variances, the variances will be
            # representative of the 200, but might not show how the models vary throughout
            # the whole cluster, which could provide better information for truncating the models.
            #
            # Note - hlfsimko - 24.02.16: Maybe we can keep the full length clusters
            #                             and slice the list after truncation?
            cluster.num_clusters = ns_clusters
            cluster.models = cluster.models[0:max_cluster_size]
            clusters.append(cluster)

        return clusters

    def _cluster(self, models, run_dir=None, score_type='rmsd', score_matrix=None, nproc=1):
        """
        Run spicker to cluster the models
        """
        owd = os.getcwd()
        if run_dir: self.run_dir = os.path.abspath(run_dir)
        if not self.run_dir: self.run_dir = os.path.join(owd, 'spicker')
        if not os.path.isdir(self.run_dir): os.mkdir(self.run_dir)
        os.chdir(self.run_dir)

        logger.debug("Running spicker with score_type {0} in directory: {1}".format(score_type, self.run_dir))
        logger.debug("Using executable: {0} on {1} processors".format(self.spicker_exe, nproc))

        self.score_type = score_type
        self.create_input_files(models, score_type=score_type, score_matrix=score_matrix)

        # We need special care if we are running with tm scores as we will be using the OPENMP
        # version of spicker which requires increasing the stack size on linux and setting the
        # OMP_NUM_THREADS environment variable on all platforms
        # The stack size on 64-bit linux seems to be 15Mb, so I guess asking for 50 seems reasonable
        # I'm assuming that the limit is in bytes and specified by an integer so 50Mb -> 50000000
        preexec_fn = None
        env = {'OMP_NUM_THREADS': str(nproc)}
        if sys.platform.lower().startswith('linux'):

            def set_stack():
                import resource
                stack_bytes = 50000000  # 50Mb
                resource.setrlimit(resource.RLIMIT_STACK, (stack_bytes, stack_bytes))

            preexec_fn = set_stack

        logfile = os.path.abspath("spicker.log")
        rtn = ample_util.run_command([self.spicker_exe], logfile=logfile, env=env, preexec_fn=preexec_fn)
        if not rtn == 0:
            raise RuntimeError, "Error running spicker, check logfile: {0}".format(logfile)

        # Read the log and generate the results
        self.results = self.process_log()

        # Always go back to where we started
        os.chdir(owd)
        return

    def process_log(self, logfile=None):
        """Read the spicker str.txt file and return a list of SpickerResults for each cluster.
        
        We use the R_nat value to order the files in the cluster
        """
        if not logfile: logfile = os.path.join(self.run_dir, 'str.txt')
        clusterCounts = []
        index2rcens = []

        # File with the spicker results for each cluster
        logger.debug("Processing spicker output file: {0}".format(logfile))
        f = open(logfile, 'r')
        line = f.readline()
        while line:
            line = line.strip()
            if line.startswith("#Cluster"):
                # skip 2 lines to Nstr
                f.readline()
                f.readline()
                line = f.readline().strip()
                if not line.startswith("Nstr="):
                    raise RuntimeError, "Problem reading file: {0}".format(logfile)

                ccount = int(line.split()[1])
                clusterCounts.append(ccount)

                # Loop through this cluster
                i2rcen = []
                line = f.readline().strip()
                while not line.startswith("------"):
                    fields = line.split()
                    #  i_cl   i_str  R_nat   R_cen  E    #str     traj
                    # tuple of: ( index in file , distance from centroid )
                    i2rcen.append((int(fields[5]), float(fields[3])))
                    line = f.readline().strip()
                index2rcens.append(i2rcen)
            line = f.readline()

        # Sort clusters by the R_cen - distance from cluster centroid
        for i, l in enumerate(index2rcens):
            # Sort by the distance form the centroid, so first becomes centroid
            sorted_by_rcen = sorted(l, key=lambda tup: tup[1])
            index2rcens[i] = sorted_by_rcen

        # Now map the indices to their files
        # Get ordered list of the pdb files
        flist = os.path.join(self.run_dir, 'file_list')
        pdb_list = [line.strip() for line in open(flist, 'r')]

        results = []
        # create results
        num_clusters = len(clusterCounts)
        for cluster in range(num_clusters):
            result = Cluster()
            result.cluster_method = self.cluster_method
            result.cluster_score_type = self.score_type
            result.index = cluster + 1
            result.num_clusters = num_clusters
            pdb_file = os.path.join(self.run_dir, "spicker_cluster_{0}.list".format(cluster + 1))
            with open(pdb_file, "w") as f:
                for i, (idx, rcen) in enumerate(index2rcens[cluster]):
                    pdb = pdb_list[idx - 1]
                    #if i == 0: result.cluster_centroid = pdb
                    result.models.append(pdb)
                    result.r_cen.append(rcen)
                    f.write(pdb + "\n")
            results.append(result)
        return results

    def results_summary(self):
        """Summarise the spicker results"""

        if not self.results: raise RuntimeError, "Could not find any results!"
        rstr = "---- Spicker Results ----\n\n"

        for i, r in enumerate(self.results):
            rstr += "Cluster: {0}\n".format(i + 1)
            rstr += "* number of models: {0}\n".format(r.size)
            #rstr += "* files are listed in file: {0}\n".format(r.pdb_file)
            rstr += "* centroid model is: {0}\n".format(r.centroid)
            rstr += "\n"

        return rstr

if __name__ == "__main__":
    import argparse
    #
    # Run Spicker on a directory of PDB files
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--executable', help="spicker executable to use")
    parser.add_argument('-m', '--models', required=True, help="Models to cluster")
    parser.add_argument('-s', '--score_type', default='rmsd', help="Use TM score")
    args = parser.parse_args()

    models = args.models
    if not os.path.exists(models): raise RuntimeError("Cannot find models: {0}".format(models))

    if os.path.isdir(models):
        models = [os.path.abspath(m) for m in glob.glob(os.path.join(models, "*.pdb"))]
    elif os.path.isfile(models):
        with open(models) as f:
            models = [l.strip() for l in f.readlines() if l.strip()]

    if not len(models):
        print "Cannot find any pdbs in: {0}".format(models)
        sys.exit(1)

    spicker_exe = os.path.abspath(args.executable) if args.executable else None

    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    spicker = Spickerer(spicker_exe=spicker_exe)
    spicker.cluster(models, score_type=args.score_type)
    print spicker.results_summary()
