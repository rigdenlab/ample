#!/usr/bin/env python

# python imports
import glob
import logging
import os
import re
import shutil
import subprocess

# our imports
import ample_util


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


class SpickerCluster( object ):

    def __init__(self, amoptd ):
        """Initialise from a dictionary of options"""
        
        self.rundir = amoptd['spicker_rundir']
        #self.clusterdir = amoptd['spicker_clusterdir']
        self.spicker_exe =  amoptd['spicker_exe']
        self.models_dir = amoptd['models_dir']
        self.num_clusters = amoptd['num_clusters']
        self.results = None

        self.logger = logging.getLogger()
        
    ###############
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

    #################
    def create_input_files(self):
    
        """
        jmht
        Create the input files required to run spicker
    
        path: PATH_TO_MODELS
        outpath: RunDir+'/spicker_run'
    
        (See notes in spicker.f FORTRAN file for a description of the required files)
    
        """

        #os.chdir( self.rundir )
        
        pdbname = ''
        list_string = ''
        counter = 0
    
        # Input file for spicker with coordinates of the CA atoms for each of the PDB structures
        read_out = open( 'rep1.tra1', "w")
    
        #jmht - a list of the full path of all PDBs - used so we can loop through it and copy the selected
        # ones to the relevant directory after we have run spicker - the order of these must match the order
        # of the structures in the rep1.tra1 file 
        file_list = open( 'file_list', "w")
    
        for infile in glob.glob( os.path.join(self.models_dir,  '*.pdb') ):
            
            pdbname = os.path.basename(infile)
            file_list.write(infile + '\n')
            list_string = list_string + pdbname+ '\n'
            counter +=1
            read = open(infile)
    
            length = self.get_length(infile)
            # 1st field is length, 2nd energy, 3rd & 4th don't seem to be used for anything
            read_out.write('\t' + length + '\t926.917       '+str(counter)+'       '+str(counter)+'\n')
            
            # Write out the coordinates of the CA atoms 
            for line in read:
                #print line
                pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)\s*(.\d*.\d*)')
                result = re.match(pattern, line)
                if result:
                    split = re.split(pattern, line)
                    #print line
                    #print split
                # print '\t' + split[5] + '\t' + split[6] + '\t' +split[7]
                    if split[2] == 'CA':
                        read_out.write( '     ' + split[6] + '     ' + split[7] + '     ' +split[8] + '\n' )
    
            read.close()
        file_list.close()
        #print list_string
        #make rmsinp
        #length = get_length(path + '/' + pdbname)
    
# from spicker.f
#*       'rmsinp'---Mandatory, length of protein & piece for RMSD calculation;
        rmsinp = open('rmsinp', "w")
        rmsinp.write('1  ' + length + '\n\n')
        rmsinp.write(length + '\n')
        rmsinp.close()
        
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
        tra = open('tra.in', "w")
        tra.write('1 -1 1 \nrep1.tra1')
    
        tra.close()
        # Create the file with the sequence of the PDB structures
# from spicker.f
#*       'seq.dat'--Mandatory, sequence file, for output of PDB models.
        seq = open('seq.dat', "w")
        a_pdb = open( os.path.join( self.models_dir, pdbname ), 'r' )
        for line in a_pdb:
            #print line
            pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s*(\d*)\s')
            result = re.match(pattern, line)
            if result:
                split = re.split(pattern, line)
                #print line
            # print split
                if split[2] == 'CA':
                    seq.write('\t' +split[5] + '\t' + split[3] + '\n')
    ################
    
    
    def run_spicker(self):
        """
        Run spicker to cluster the models
        """
        
        if not os.path.isdir( self.rundir ):
            os.mkdir( self.rundir )
        os.chdir(self.rundir)
        
        self.logger.debug("Running spicker in directory: {0}".format( self.rundir ) )
        self.create_input_files()
        ample_util.run_command([ self.spicker_exe ], logfile="spicker.log" )
    
        # Read the log and generate the results
        results = self.process_log()

	# Check we have enough clusters
        if len(results) < self.num_clusters:
            msg = "Only {0} clusters returned from Spicker cannot process {1} clusters!\n".format( len(results),self.num_clusters )
            self.logger.critical( msg )
            raise RuntimeError,msg
        
        # Loop through each cluster copying the files as we go
        # We only process the clusters we will be using
        MAXCLUSTER=200 # max number of pdbs in a cluster
        for cluster in range( self.num_clusters ):
                
            result = results[ cluster ]
            result.pdb_file = os.path.join( self.rundir, "spicker_cluster_{0}.list".format(cluster+1)  )
            
            f = open( result.pdb_file, "w" )
            for i, pdb in enumerate( result.rosetta_pdb ):
                
                if i > MAXCLUSTER:
                    result.cluster_size = MAXCLUSTER
                    break
                
                result.pdb_list.append( pdb )
                f.write( pdb + "\n" )
                
                if i == 0:
                    result.cluster_centroid = pdb
            
            f.close()
        
        self.results = results
        return
                
    def process_log( self, logfile=None ):
        """Read the spicker str.txt file and return a list of SpickerResults for each cluster.
        
        We use the R_nat value to order the files in the cluster
        """
        
        if not logfile:
            logfile = os.path.join(self.rundir, 'str.txt')
            
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
        flist = os.path.join( self.rundir, 'file_list')
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
        
#################
#def check_pid_by_kill(pid):
#    """ Check For the existence of a unix pid. """
#    try:
#        os.kill(pid, 0)
#    except OSError:
#        return False
#    else:
#        return True
###################
#def check_pid(pid):
# if os.path.exists('/proc/'+str(pid)):
#   return True
# else:
#   return False


if __name__ == "__main__":

    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)

    root = "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0"
    optd = { 'spicker_rundir' : os.path.join( root, 'spicker_run'),
            'spicker_clusterdir' : os.path.join( root, 'S_clusters'),
            'spicker_exe' : '/opt/spicker/spicker',
            'models_dir': os.path.join( root, 'models' ),
            'num_clusters' : 3
            }

    spicker = SpickerCluster( optd )
    spicker.run_spicker()
    print spicker.results_summary()
    #spicker_results = spicker.get_results()

    #run_spicker(models, spicker_run_dir, spickerexe, no_clusters_sampled, overpath)
