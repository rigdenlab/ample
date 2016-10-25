#!/usr/bin/env python

# python imports
import glob
import logging
import os
import re
import shutil
import sys

# our imports
try:
    from ample.util import ample_util
except:
    ample_util = None
    
logger = logging.getLogger(__name__)

class SpickerResult(object):
    """
    A class to hold the result of running Spicker
    """

    def __init__(self):

        self.pdb_file = None  # Path to a list of the pdbs for this cluster
        self.cluster_size = None
        self.cluster_centroid = "N/A"
        self.pdbs = []  # ordered list of the pdbs in their results directory
        self.r_cen = []  # ordered list of the distance from the cluster centroid for each pdb
        return

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
            shutil.copy(score_matrix, os.path.join(self.run_dir,'score.matrix'))
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
                        pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)')
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
    
    def cluster(self, models, run_dir=None, score_type='rmsd', score_matrix=None, nproc=1):
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
        
        self.create_input_files(models, score_type=score_type, score_matrix=score_matrix)
        
        # We need special care if we are running with tm scores as we will be using the OPENMP
        # version of spicker which requires increasing the stack size on linux and setting the 
        # OMP_NUM_THREADS environment variable on all platforms
        # The stack size on 64-bit linux seems to be 15Mb, so I guess asking for 50 seems reasonable
        # I'm assuming that the limit is in bytes and specified by an integer so 50Mb -> 50000000
        preexec_fn=None
        env = None
        if score_type == 'tm':
            env = { 'OMP_NUM_THREADS' : str(nproc)}
            if sys.platform.lower().startswith('linux'):
                def set_stack():
                    import resource
                    stack_bytes = 50000000 # 50Mb
                    #resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY,resource.RLIM_INFINITY))
                    resource.setrlimit(resource.RLIMIT_STACK, (stack_bytes,stack_bytes))
                preexec_fn=set_stack

        logfile = os.path.abspath("spicker.log")
        rtn = ample_util.run_command([self.spicker_exe], logfile=logfile, env=env, preexec_fn=preexec_fn)
        if not rtn == 0:
            raise RuntimeError,"Error running spicker, check logfile: {0}".format(logfile)
    
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
        pdb_list = [ line.strip() for  line in open(flist , 'r') ]
        
        results = []
        # create results
        for cluster in range(len(clusterCounts)):
            result = SpickerResult()
            result.cluster_size = clusterCounts[ cluster ]
            result.pdb_file = os.path.join(self.run_dir, "spicker_cluster_{0}.list".format(cluster + 1))
            with open(result.pdb_file, "w") as f:
                for i, (idx, rcen) in enumerate(index2rcens[ cluster ]):
                    pdb = pdb_list[idx - 1]
                    if i == 0: result.cluster_centroid = pdb
                    result.pdbs.append(pdb)
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
            rstr += "* number of models: {0}\n".format(r.cluster_size)
            rstr += "* files are listed in file: {0}\n".format(r.pdb_file)
            rstr += "* centroid model is: {0}\n".format(r.cluster_centroid)
            rstr += "\n"
            
        return rstr

if __name__ == "__main__":
    import argparse, subprocess, tempfile
    
    # For running as a stand-alone script
    def run_command(cmd, logfile=None, directory=None, dolog=True, stdin=None, check=False):
        """Execute a command and return the exit_error code.
    
        We take care of outputting stuff to the logs and opening/closing logfiles
    
        Args:
        cmd - command to run as a list
        stdin - a string to use as stdin for the command
        logfile (optional) - the path to the logfile
        directory (optional) - the directory to run the job in (cwd assumed)
        dolog: bool - whether to output info to the system log
        """
    
        assert type(cmd) is list
    
        if not directory:
            directory = os.getcwd()
    
        if dolog:
            logging.debug("In directory {0}\nRunning command: {1}".format(directory, " ".join(cmd)))
    
        if logfile:
            if dolog:
                logging.debug("Logfile is: {0}".format(logfile))
            logf = open(logfile, "w")
        else:
            logf = tempfile.TemporaryFile()
            
        if stdin != None:
            stdinstr = stdin
            stdin = subprocess.PIPE
    
        # Windows needs some special treatment
        kwargs = {}
        if os.name == "nt":
            kwargs = { 'bufsize': 0, 'shell' : "False" }
        p = subprocess.Popen(cmd, stdin=stdin, stdout=logf, stderr=subprocess.STDOUT, cwd=directory, **kwargs)
    
        if stdin != None:
            p.stdin.write(stdinstr)
            p.stdin.close()
            if dolog:
                logging.debug("stdin for cmd was: {0}".format(stdinstr))
    
        p.wait()
        logf.close()
        
        return p.returncode
    
    # Mock up ample_util for when we don't have CCP4 installed
    if ample_util is None:
        class Tmp(object):pass
        ample_util = Tmp()
        ample_util.run_command = run_command
    
    #
    # Run Spicker on a directory of PDB files
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--executable',
                        help="spicker executable to use")
    parser.add_argument('-m', '--models', required=True,
                        help="Models to cluster")
    parser.add_argument('-s', '--score_type', default='rmsd',
                        help="Use TM score")
    args = parser.parse_args()
    
    models = args.models
    if not os.path.exists(models): raise RuntimeError("Cannot find models: {0}".format(models))
    
    if os.path.isdir(models):
        models = [ os.path.abspath(m) for m in glob.glob(os.path.join(models, "*.pdb")) ]
    elif os.path.isfile(models):
        with open(models) as f:
            models = [ l.strip() for l in f.readlines() if l.strip() ]
        
    if not len(models):
        print "Cannot find any pdbs in: {0}".format(models)
        sys.exit(1)

    spicker_exe=None
    if args.executable:
       spicker_exe=os.path.abspath(args.executable)
        
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    spicker = Spickerer(spicker_exe=spicker_exe)
    spicker.cluster(models, score_type=args.score_type)
    print spicker.results_summary()
