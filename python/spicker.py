#!/usr/bin/env python

# python imports
import glob
import logging
import os
import re

# our imports
try:
    import ample_util
except:
    ample_util=None

class SpickerResult( object ):
    """
    A class to hold the result of running Spicker
    """

    def __init__(self):

        self.pdb_file = None # Path to a list of the pdbs for this cluster
        self.cluster_size = None
        self.cluster_centroid = "N/A"
        self.pdb_list = [] # ordered list of the pdbs in their results directory
        self.rosetta_pdb = [] # ordered list of the pdbs in the rosetta directory
        self.r_cen = [] # ordered list of the distance from the cluster centroid for each pdb
        return

class Spickerer( object ):

    def __init__(self, spicker_exe=None,run_dir=None):
        """Initialise from a dictionary of options"""
        
        if not spicker_exe: spicker_exe=os.path.join(os.environ['CCP4'], 'bin', 'spicker')
        if not (os.path.exists(spicker_exe) and os.access(spicker_exe, os.X_OK)):
            raise RuntimeError,"Cannot find spicker executable: {0}".format(spicker_exe)
        
        self.spicker_exe =  spicker_exe
        self.run_dir = run_dir
        self.num_clusters = None
        self.results = None
        self.max_cluster_size = 200
        self.logger = logging.getLogger()
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
                    counter+=1
    
        # print counter
        return str(counter)

    def create_input_files(self,models):
    
        """
        jmht
        Create the input files required to run spicker
    
        (See notes in spicker.f FORTRAN file for a description of the required files)
    
        """
        if not len(models):
            msg = "no models provided!"
            self.logger.critical(msg)
            raise RuntimeError,msg
        
        # read_out - Input file for spicker with coordinates of the CA atoms for each of the PDB structures
        #
        # file_list - a list of the full path of all PDBs - used so we can loop through it and copy the selected
        # ones to the relevant directory after we have run spicker - the order of these must match the order
        # of the structures in the rep1.tra1 file 
        
        list_string = ''
        counter = 0
        with open( 'rep1.tra1', "w") as read_out, open( 'file_list', "w") as file_list:
            for infile in models:
                pdbname = os.path.basename(infile)
                file_list.write(infile + '\n')
                list_string = list_string + pdbname+ '\n'
                counter +=1
        
                length = self.get_length(infile)
                # 1st field is length, 2nd energy, 3rd & 4th don't seem to be used for anything
                read_out.write('\t' + length + '\t926.917       '+str(counter)+'       '+str(counter)+'\n')
                with open(infile) as read:
                    # Write out the coordinates of the CA atoms 
                    for line in read:
                        pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)')
                        result = re.match(pattern, line)
                        if result:
                            split = re.split(pattern, line)
                            if split[2] == 'CA':
                                read_out.write( '     ' + split[6] + '     ' + split[7] + '     ' +split[8] + '\n' )
    
        # from spicker.f
        #*       'rmsinp'---Mandatory, length of protein & piece for RMSD calculation;
        with open('rmsinp', "w") as rmsinp:
            rmsinp.write('1  ' + length + '\n\n')
            rmsinp.write(length + '\n')
        
        #make tra.in
        # from spicker.f
        #*       'tra.in'---Mandatory, list of trajectory names used for clustering.
        #*                  In the first line of 'tra.in', there are 3 parameters:
        #*                  par1: number of decoy files
        #*                  par2: 1, default cutoff, best for decoys from template-based
        #*                           modeling;
        #*                       -1, cutoff based on variation, best for decoys from
        #*                           ab initio modeling.
        #*                  par3: 1, closc from all decoys; -1, closc clustered decoys
        #*                  From second lines are the file names which contain coordinates
        #*                  of 3D structure decoys. All these files are mandatory
        with open('tra.in', "w") as tra:
            tra.write('1 -1 1 \nrep1.tra1\n')
    
        # Create the file with the sequence of the PDB structures
        # from spicker.f
        #*       'seq.dat'--Mandatory, sequence file, for output of PDB models.
        with open('seq.dat', "w") as seq, open(models[0], 'r') as a_pdb:
            for line in a_pdb:
                pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s*(\d*)\s')
                result = re.match(pattern, line)
                if result:
                    split = re.split(pattern, line)
                    if split[2] == 'CA':
                        seq.write('\t' +split[5] + '\t' + split[3] + '\n')
        return
    
    def cluster(self, models, run_dir=None, num_clusters=1):
        """
        Run spicker to cluster the models
        """
        
        self.num_clusters = num_clusters
        
        if run_dir: self.run_dir = run_dir
        if not self.run_dir: self.run_dir = os.path.join(os.getcwd(),'spicker')
        if not os.path.isdir(self.run_dir): os.mkdir(self.run_dir)
        os.chdir(self.run_dir)
        
        self.logger.debug("Running spicker in directory: {0}".format(self.run_dir))
        self.create_input_files(models)
        ample_util.run_command([self.spicker_exe], logfile="spicker.log")
    
        # Read the log and generate the results
        results = self.process_log()

        # Check we have enough clusters
        if len(results) < self.num_clusters:
            msg = "Only {0} clusters returned from Spicker cannot process {1} clusters!\n".format(len(results),self.num_clusters)
            self.logger.critical( msg )
            raise RuntimeError,msg
        
        # Loop through each cluster copying the files as we go
        # We only process the clusters we will be using
        for cluster in range( self.num_clusters ):
                
            result = results[ cluster ]
            result.pdb_file = os.path.join(self.run_dir, "spicker_cluster_{0}.list".format(cluster+1))
            
            with open( result.pdb_file, "w" ) as f:
                for i, pdb in enumerate(result.rosetta_pdb):
                    
                    if i > self.max_cluster_size:
                        result.cluster_size = self.max_cluster_size
                        break
                    result.pdb_list.append( pdb )
                    f.write( pdb + "\n" )
                    if i == 0: result.cluster_centroid = pdb
        self.results = results
        return
                
    def process_log( self, logfile=None ):
        """Read the spicker str.txt file and return a list of SpickerResults for each cluster.
        
        We use the R_nat value to order the files in the cluster
        """
        
        if not logfile:
            logfile = os.path.join(self.run_dir, 'str.txt')
            
        clusterCounts = []
        index2rcens = []
        
        # File with the spicker results for each cluster
        self.logger.debug("Processing spicker output file: {0}".format(logfile))
        f = open( logfile, 'r' )
        line = f.readline()
        while line:
            line = line.strip()
            
            if line.startswith("#Cluster"):
                ncluster = int( line.split()[1] )
                
                # skip 2 lines to Nstr
                f.readline()
                f.readline()
                
                line = f.readline().strip()
                if not line.startswith("Nstr="):
                    raise RuntimeError,"Problem reading file: {0}".format( logfile )
                
                ccount = int( line.split()[1] )
                clusterCounts.append( ccount )
                
                # Loop through this cluster
                i2rcen = []
                line = f.readline().strip()
                while not line.startswith("------"):
                    fields = line.split()
                    #  i_cl   i_str  R_nat   R_cen  E    #str     traj
                    # tuple of: ( index in file , distance from centroid )
                    i2rcen.append( ( int(fields[5]), float( fields[3] ) ) )
                    line = f.readline().strip()
                    
                index2rcens.append( i2rcen )
            
            line = f.readline()
                
        # Sort clusters by the R_cen - distance from cluster centroid
        for i,l in enumerate(index2rcens):
            # Sort by the distance form the centroid, so first becomes centroid
            sorted_by_rcen = sorted(l, key=lambda tup: tup[1])
            index2rcens[i] = sorted_by_rcen
    
        # Now map the indices to their files
        
        # Get ordered list of the pdb files
        flist = os.path.join( self.run_dir, 'file_list')
        pdb_list = [ line.strip() for  line in open( flist , 'r' ) ]
        
        results = []
        # create results
        for c in range( len( clusterCounts ) ):
            r = SpickerResult()
            r.cluster_size = clusterCounts[ c ]
            for i, rcen in index2rcens[ c ]:
                pdb = pdb_list[i-1]
                r.rosetta_pdb.append( pdb )
                r.r_cen.append( rcen )
            
            results.append( r )
            
        return results
        
    def results_summary(self):
        """Summarise the spicker results"""
        
        if not self.results:
            raise RuntimeError, "Could not find any results!"
        
        
        rstr = "---- Spicker Results ----\n\n"
        
        for i, r in enumerate( self.results ):
            rstr += "Cluster: {0}\n".format(i+1)
            rstr += "* number of models: {0}\n".format( r.cluster_size )
            if i <= self.num_clusters-1:
                rstr += "* files are listed in file: {0}\n".format( r.pdb_file )
                rstr += "* centroid model is: {0}\n".format( r.cluster_centroid )
            rstr += "\n"
            
        return rstr



    
if __name__ == "__main__":
    
    import subprocess,tempfile
    
    # For running as a stand-alone script
    def run_command(cmd, logfile=None, directory=None, dolog=True, stdin=None, check=False):
        """Execute a command and return the exit code.
    
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
            logging.debug("In directory {0}\nRunning command: {1}".format( directory, " ".join(cmd)  ) )
    
        if logfile:
            if dolog:
                logging.debug("Logfile is: {0}".format( logfile ) )
            logf = open( logfile, "w" )
        else:
            logf = tempfile.TemporaryFile()
            
        if stdin != None:
            stdinstr = stdin
            stdin = subprocess.PIPE
    
        # Windows needs some special treatment
        kwargs = {}
        if os.name == "nt":
            kwargs = { 'bufsize': 0, 'shell' : "False" }
        p = subprocess.Popen( cmd, stdin=stdin, stdout=logf, stderr=subprocess.STDOUT, cwd=directory, **kwargs )
    
        if stdin != None:
            p.stdin.write( stdinstr )
            p.stdin.close()
            if dolog:
                logging.debug("stdin for cmd was: {0}".format( stdinstr ) )
    
        p.wait()
        logf.close()
        
        return p.returncode
    
    class Tmp(object):pass
    ample_util=Tmp()
    ample_util.run_command=run_command
    
    #
    # Run Spicker on a directory of PDB files
    #
    import sys
    if len(sys.argv) != 2:
        print "Usage is {0} <directory_of_pdbs>".format( sys.argv[0] )
        sys.exit(1)
        
    models_dir = os.path.abspath(sys.argv[1])
    if not os.path.isdir(models_dir):
        print "Cannot find directory: {0}".format(models_dir)
        sys.exit(1)
    models=glob.glob(os.path.join(models_dir,"*.pdb"))
    if not len(models):
        print "Cannot find any pdbs in: {0}".format(models_dir)
        sys.exit(1)
        
#     spicker_exe = ample_util.find_exe("spicker")
#     if not spicker_exe:
#         print "Cannot find spicker executable in path!"
#         sys.exit(1)
        
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    spicker = Spickerer()
    spicker.cluster(models,num_clusters=10)
    print spicker.results_summary()
