'''
Created on Apr 18, 2013

@author: jmht
'''

import copy
import cPickle
import logging
import os
import re
import shutil
import sys

# our imports
import ample_util
import cluster_with_MAX
import pdb_edit
import printTable
import run_spicker

class EnsembleData(object):
    """Class to hold data about an ensemble"""
    
    def __init(self):
        
        self.name = None 
        self.num_models = None
        self.num_residues = None
        self.num_atoms = None
        self.residues = []
        
        self.side_chain_treatment = None
        self.radius_threshold = None
        self.truncation_threshold = None
        
        self.pdb = None # path to the ensemble file


class Ensembler(object):
    """Class to generate ensembles from cluster of models"""
    
    def __init__(self):
        """Initialise"""
        
        # Directory where all files generated
        self.work_dir = None
        
        # Directory where the final ensembles end up
        self.ensemble_dir = None
        
        # cluster of models to work from
        self.cluster_models = []
        
        # Percent of residues to keep under each truncation level
        self.percent= None
        
        # theseus variance file
        self.variance_log = None
        
        # List of the thresholds we truncate at
        self.truncation_thresholds = []
        
        # List of directories where the truncated files are (must match truncation_thresholds )
        self.truncation_dirs = []
        
        # List of lists of the residues kept under each truncation threshold - ordered by truncation_thresholds
        self.truncation_residues = []
        
        # List of list of the truncated files under each truncation threshold - ordered by truncation_thresholds
        self.truncated_models = []
        
        # radius thresholds to cluster the models
        self.radius_thresholds = [ 1, 2, 3 ]
        
        # The list of ensembles before side-chain treatment
        self.truncated_ensembles = []
        
        # The list of final ensembles
        self.ensembles = []
        
        # Programs we use
        self.theseus_exe = None
        
        self.maxcluster_exe = None

    def pdb_list(self):
        """Return the final ensembles as a list of pdb files"""
        ensembles = []
        for ensemble in self.ensembles:
            ensembles.append( ensemble.pdb )
        return ensembles

    def edit_sidechains(self):
        """Take the ensembles and give them the 3 sidechain treatments"""
        
        pdbed = pdb_edit.PDBEdit()
        
        
        for trunc_ensemble in self.truncated_ensembles:
            
            # 3 different side-chain treatments
            
            # 1. all_atom
            ensemble = copy.deepcopy( trunc_ensemble )
            ensemble.side_chain_treatment = "all_atom"
            ensemble.side_chain_treatment = "All_atom"
            ensemble.name = "{0}_{1}".format( ensemble.side_chain_treatment, trunc_ensemble.name )
            ensemble.pdb = os.path.join( self.ensemble_dir, ensemble.name+".pdb" )
            # For all atom just copy the file
            shutil.copy2( trunc_ensemble.pdb, ensemble.pdb )
            self.ensembles.append( ensemble )
            
            # 2. Reliable side chains
            ensemble = copy.deepcopy( trunc_ensemble )
            ensemble.side_chain_treatment = "reliable_sidechains"
            ensemble.side_chain_treatment = "SCWRL_reliable_sidechains"
            ensemble.name = "{0}_{1}".format( ensemble.side_chain_treatment, trunc_ensemble.name )
            ensemble.pdb = os.path.join( self.ensemble_dir, ensemble.name+".pdb" )
            # Do the edit
            pdbed.reliable_sidechains( trunc_ensemble.pdb, ensemble.pdb )
            self.ensembles.append( ensemble )            
            
            # 3. backbone
            ensemble = copy.deepcopy( trunc_ensemble )
            ensemble.side_chain_treatment = "poly_ala"
            ensemble.name = "{0}_{1}".format( ensemble.side_chain_treatment, trunc_ensemble.name )
            ensemble.pdb = os.path.join( self.ensemble_dir, ensemble.name+".pdb" )
            # Do the edit
            pdbed.backbone( trunc_ensemble.pdb, ensemble.pdb )
            self.ensembles.append( ensemble )            
            
        return
    
    def generate_ensembles(self, cluster_models=None, root_dir=None, ensemble_id=None, percent=None ):
        """Generate the ensembles for this list of cluster models"""
        
        assert self.maxcluster_exe and self.theseus_exe, "Must set the path to maxcluster and theseus!"
        
        self.cluster_models = cluster_models
        self.percent = percent
        
        # Create the directories we'll be working in
        self.work_dir =  os.path.join( root_dir, 'fine_cluster_{0}'.format( ensemble_id ) )
        os.mkdir( self.work_dir )
        os.chdir( self.work_dir )
        
        self.ensemble_dir = os.path.join( root_dir, 'ensembles_{0}'.format( ensemble_id ) )
        os.mkdir( self.ensemble_dir )
        
        # calculate the thresholds
        self.generate_thresholds()
        
        # Truncate our models based on calculated thresholds
        self.truncate_models()
        
        # Generate the ensembles for the different truncation levels
        self.make_truncated_ensembles()
        
        # Undertake the 3 side-chain treatments
        self.edit_sidechains()


    def generate_thresholds(self):
        """
        Calculate the residue variance thresholds that will keep self.percent residues for each truncation level
        """
        #--------------------------------
        # choose threshold type
        #-------------------------------
        FIXED_INTERVALS=False
        if FIXED_INTERVALS:
            self.thresholds = [ 1, 1.5, 2 , 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8 ]
            logging.debug("Got {0} thresholds: {1}".format( len(self.thresholds), self.thresholds ))
            return
        
        #--------------------------------
        # get variations between pdbs
        #--------------------------------
        cmd = [ self.theseus_exe, "-a0" ] + self.cluster_models
        retcode = ample_util.run_command( cmd, logfile="theseus.log" )
    
        # print theseus_out
        self.variance_log = os.path.join( self.work_dir, 'theseus_variances.txt' )
        theseus_out = open( self.variance_log  )
        
        # List of variances ordered by residue index
        var_list=[]
        
        #build up a list of the variances of each residue
        for i, line in enumerate(theseus_out):
            
            # Skip header
            if i==0:
                continue
            #print line
            # Different versions of theseus may have a RES card first, so remove
            line = re.sub('RES\s*', '', line)
            pattern = re.compile('^(\d*)\s*(\w*)\s*(\d*)\s*(\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)')
            result = pattern.match(line)
            if result:
                #print line
                seq = re.split(pattern, line)
                var_list.append(float(seq[4]))
        
        theseus_out.close()
    
        length = len(var_list)
        if length == 0:
            msg = "Error generating thresholds, got len: {0}".format(length)
            logging.critical(msg)
            raise RuntimeError,msg
        
        
        percent_interval=( float(length)/100 ) *float(self.percent)
#        lowest=min(var_list)
#        highest=max(var_list)
#        print length
#        print lowest, highest
#        print int( percent_interval )
    
        ## try to find intervals for truncation
        try_list=copy.deepcopy(var_list)
        try_list.sort()
        #print try_list
    
        # print list(chunks(try_list, int(percent_interval)))
        # For chunking list
        def chunks(a_list, percent_interval):
            for i in xrange(0, len(a_list), percent_interval):
                yield a_list[i:i+percent_interval ]
    
        for x in list( chunks(try_list, int(percent_interval) ) ):
            #print x, x[-1]
            # For some cases, multiple residues share the same variance so we don't create a separate thereshold
            if x[-1] not in self.truncation_thresholds:
                self.truncation_thresholds.append(x[-1])
            
        logging.debug("Got {0} thresholds: {1}".format( len(self.truncation_thresholds), self.truncation_thresholds ))
        
        return
        
    def make_list_to_keep(self, threshold ):
        """Make a list of residues to keep under variance threshold
        INPUTS:
        threshold: the threshold variance
        
        RETURNS:
        a list of the residue indexes
        """
        add_list =[]
        
        theseus_out = open( self.variance_log )
        for line in theseus_out:
            # for alternate versions of theseus remove RES cards
            line = re.sub('RES\s*', '', line)
            pattern = re.compile('^(\d*)\s*(\w*)\s*(\d*)\s*(\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)')
            result = pattern.match(line)
            if result:
                seq = re.split(pattern, line)
                #print seq
                if (seq[6]) != 'T':
        
                    if float(seq[4]) <= float(threshold):
                        if seq[1] != '':
                    #  print 'GOT ' + str(seq[4]) + '  thresh ' + str(THESEUS_threthold) + ' KEEP ' + str(seq[3])
                            add_list.append( int(seq[3]) )
                            
        theseus_out.close()
        
        #print add_list
        return add_list
          
    def make_truncated_ensembles( self  ):
        """
        Loop through each truncation level and use maxcluster to cluster the models
        according to the three radius thresholds in self.radius_thresholds
    
        If there are more than 2 clusters returned by maxcluster use theseus to align them
        
        Generate an ensemble for 
    
        """
    
        # Loop through all truncation levels
        for tcount, truncation_threshold in enumerate( self.truncation_thresholds ):
            
            # change to the truncation directory
            truncation_dir = self.truncation_dirs[tcount]
            os.chdir( truncation_dir )
            
            # Create the list of files for maxcluster            
            fname = os.path.join( truncation_dir, "files.list" )
            f = open( fname, 'w' )
            f.write( "\n".join( self.truncated_models[tcount] )+"\n" )
            f.close()
    
            # Loop through the radius thresholds
            for radius in self.radius_thresholds:
                
                logging.debug("Clustering files under radius: {0}".format( radius ) )
                
                ensemble_dir = os.path.join( truncation_dir, 'fine_clusters_'+str(radius)+'_ensemble' )
                os.mkdir( ensemble_dir )
                os.chdir( ensemble_dir )
        
                # A list of the files that have been clustered together with maxcluster
                cluster_files = cluster_with_MAX.cluster_with_MAX_FAST( fname, radius, self.maxcluster_exe )  # use fastest method
                logging.debug("Maxcluster clustered {0} files".format ( len( cluster_files ) ) )
                if cluster_files < 2:
                    logging.info( 'Could not cluster files using radius {0}'.format( radius ) )
                    continue
                
                # Write out the files for reference
                file_list = "maxcluster_radius_{0}_files.list".format( radius )
                f = open(file_list, "w")
                for c in cluster_files:
                    f.write( c+"\n")
                f.write("\n")
                f.close()
                
                # Restrict cluster to 30
                if len( cluster_files ) > 30:
                    logging.debug("More than 30 files clustered so truncating list to first 30")
                    cluster_files = cluster_files[:30]
                    
                # For naming all files
                basename='trunc_{0}_rad_{1}'.format( truncation_threshold, radius ) 
        
                # Run theseus to generate a file containing the aligned clusters
                cmd = [ self.theseus_exe, "-r", basename, "-a0" ] + cluster_files
                retcode = ample_util.run_command( cmd, logfile=basename+"_theseus.log" )
                
                # jmht - the previous Align_rosetta_fine_clusters_with_theseus routine was running theseus twice and adding the averaged
                # ensemble from the first run to the ensemble. This seems to improve the results for the TOXD test case - maybe something to
                # look at?
                #cmd = [ self.theseus_exe, "-r", basename, "-a0" ] + cluster_files + [ basename+"run1_ave.pdb" ]
                #retcode = ample_util.run_command( cmd, logfile=basename+"_theseus.log" )
                
                # Check if theseus worked and if so rename the file with the aligned files and append the path to the ensembles
                cluster_file = os.path.join(  ensemble_dir, basename+'_sup.pdb' )
                ensemble_file = os.path.join( ensemble_dir, basename+'.pdb' )
        
                #same from here     
                if not os.path.exists( cluster_file ):
                    logging.info( 'Could not create ensemble for files (models too diverse): {0}'.format( ensemble_file ) )
                    continue

                # Rename the file                
                shutil.move( cluster_file, ensemble_file )
                
                # Generate the ensemble
                ensemble = EnsembleData()
                ensemble.name = basename
                ensemble.num_models = len( cluster_files )
                ensemble.residues = self.truncation_residues[ tcount ]
                ensemble.num_residues = len( ensemble.residues )
                ensemble.truncation_threshold = truncation_threshold
                ensemble.num_truncated_models = len( self.truncated_models[ tcount ] )
                ensemble.radius_threshold = radius
                ensemble.pdb = ensemble_file

                self.truncated_ensembles.append( ensemble )
                
        return
    
    def truncate_models(self):
        """Truncate the models according to the calculated thresholds
        
        Make a list of the residues to keep and then create directories for each truncation level
        and copy the truncated pdbs into them
        """
    
        #-------------------------------
        #truncate
        #----------------------------------
        truncate_log = open( os.path.join( self.work_dir, 'truncate.log'), "w")
        truncate_log.write('This is the number of residues kept under each truncation threshold\n\nthreshold\tnumber of residues\n')
        
        for threshold in self.truncation_thresholds:
    
            # Get a list of the indexes of the residues to keep
            add_list = self.make_list_to_keep( threshold )
            self.truncation_residues.append( add_list )
            
            truncate_log.write( str(threshold) +'\t' + str(len(add_list)) + '\n' )
            logging.info( 'Keeping: {0} residues at truncation level: {1}'.format( len(add_list), threshold ) )
            logging.debug( 'The following residues will be kept: {0}'.format( add_list ) )
            
            if len(add_list) < 1:
                #######   limit to size of pdb to keep, if files are too variable, all residues are removed
                logging.debug( 'No residues kept at this truncation level.' )
                continue
    
            trunc_out = os.path.join( self.work_dir, 'trunc_files_' + str(threshold) )
            os.mkdir(trunc_out)
            logging.info( 'truncating at: {0} in directory {1}'.format( threshold, trunc_out ) )
            self.truncation_dirs.append( trunc_out )
            
            # list of models for this truncation level
            truncated_models = []
            pdbed = pdb_edit.PDBEdit()
            for infile in self.cluster_models:
                pdbname = os.path.basename( infile )
                pdbout = os.path.join( trunc_out, pdbname )
    
                # Loop through PDB files and create new ones that only contain the residues left after truncation
                pdbed.select_residues( inpath=infile, outpath=pdbout, residues=add_list )
                
                truncated_models.append( pdbout )
            
            self.truncated_models.append( truncated_models )
                
        return
    
def create_ensembles( amoptd ):
    """Create the ensembles using the values in the amoptd dictionary"""

    #---------------------------------------
    # Generating  ensembles
    #---------------------------------------
    logger = logging.getLogger()
    logger.info('----- Clustering models --------')
    
    # Spicker Alternative for clustering then MAX
    amoptd['spicker_rundir'] = os.path.join( amoptd['work_dir'], 'spicker_run')
    spickerer = run_spicker.SpickerCluster( amoptd )
    spickerer.run_spicker()

    logger.info( spickerer.results_summary() )

    amoptd['spicker_results'] = spickerer.results
    
    # Generate list of clusters to sample
    # Done like this so it's easy to run specific clusters
    cluster_nums = [ i for i in range(1,amoptd['num_clusters']+1)]
    
    logger.info('Clustering Done. Using the following clusters: {0}'.format(cluster_nums) )
    
    amoptd['ensemble_results'] = [] 
    for cluster in cluster_nums:
        
        logger.info('----- Truncating models for cluster {0} --------'.format(cluster)   )
        
        # Get list of models from spicker results
        cluster_models = spickerer.results[ cluster-1 ].pdb_list
        
        if len( cluster_models ) < 2:
            msg = "Could not create ensemble for cluster {0} as less than 2 models!".format( cluster )
            logger.critical(msg)
            raise RuntimeError, msg
        
        ensembler = Ensembler()
        ensembler.maxcluster_exe = amoptd['maxcluster_exe']
        ensembler.theseus_exe =  amoptd['theseus_exe']
        ensembler.generate_ensembles( cluster_models=cluster_models,
                                      root_dir=amoptd['work_dir'],
                                      ensemble_id=cluster,
                                      percent=amoptd['percent'] )
        # Add results to amopt
        amoptd['ensemble_results'].append( ensembler.ensembles )
        
        # List of pdbs
        #final_ensembles = ensembler.pdb_list()
        
        # Prune down to the top model        
        if amoptd['top_model_only']:
            raise RuntimeError, "Need to fix top_model_only"
            final_ensembles = truncateedit_MAX.One_model_only( final_ensembles, amoptd['work_dir'] )

        # Add to the total list
        #ensembles.append( final_ensembles )
        
        logger.info("Truncating Done for cluster {0}".format( cluster ) )
        logger.info('Created {0} ensembles'.format( len(  amoptd['ensemble_results'][-1]  ) ) )
        
    # Write out pickle file
    f = open( amoptd['results_path'], 'w' )
    cPickle.dump( amoptd, f )
    f.close()
    logging.info("Saved results as file: {0}\n".format( amoptd['results_path'] ) )
        
    return

def ensemble_summary( amoptd ):
    """Print a summary of the ensembling process"""
    
    tableFormat = printTable.Table()
    rstr=""
    
    if amoptd.has_key('spicker_results'):
        rstr = "---- Clustering Results ----\n\n"
        tdata = [ ( "Cluster", "Num Models", "Centroid" ) ]
        for i, r in enumerate( amoptd['spicker_results'] ):
            tdata.append( ( i+1, r.cluster_size, r.cluster_centroid ) )
        rstr += tableFormat.pprint_table( tdata )
    
    if amoptd.has_key('ensemble_results') and len( amoptd['ensemble_results'] ):
        rstr += "\n---- Truncation Results ----\n\n"
        # Reconstruct trunction data
        for i, clusterEnsemble in enumerate( amoptd['ensemble_results'] ):
            
            rstr += "---- Cluster {0} ----\n".format( i+1 )
            nensembles = 0
            truncation_thresholds = {}
            
            for ensemble in clusterEnsemble:
                
                nensembles += 1
                if ensemble.truncation_threshold not in truncation_thresholds.keys():
                    truncation_thresholds[ ensemble.truncation_threshold ] = ( ensemble.num_residues,
                                                                               { ensemble.radius_threshold : ensemble.num_models } )
                else:
                    # Already got this truncation level so add radius
                    if ensemble.radius_threshold not in truncation_thresholds[ ensemble.truncation_threshold ][1].keys():
                        truncation_thresholds[ ensemble.truncation_threshold ][1][ ensemble.radius_threshold ] = ensemble.num_models
            
            
            rstr += "\n"
            tdata = [ ( "Truncation Threshold", "Num Residues", "Radius Treshold", "Num Models" ) ]
            for threshold in sorted( truncation_thresholds.keys() ):
                nresidues = truncation_thresholds[ threshold ][0] 
                for i, radius in enumerate( sorted( truncation_thresholds[ threshold ][1].keys() ) ):
                    nmodels = truncation_thresholds[ threshold ][1][radius ]
                    if i == 0:
                        tdata.append( ( threshold, nresidues, radius, nmodels  ) )
                    else:
                        tdata.append( ( "", "", radius, nmodels  ) )
                    
            rstr += tableFormat.pprint_table( tdata )
            
            rstr += "\nGenerated {0} ensembles\n\n".format( nensembles )
                
    return rstr

#    def XtestEnsemble():
#        logging.basicConfig()
#        logging.getLogger().setLevel(logging.DEBUG)
#    
#        ensembler = Ensembler()
#        ensembler.maxcluster_exe = "/opt/maxcluster/maxcluster"
#        ensembler.theseus_exe = "/opt/theseus_src/theseus"
#        ensembler.work_dir = "/opt/ample-dev1/python/TEST"
#        os.chdir( ensembler.work_dir    )
#        cf = "/home/Shared/TM/3LBW/ENSEMBLES_0/spicker_run/spicker_cluster_1.list"
#        ensembler.cluster_models = [ re.sub( "^/gpfs/home/HCEA041/djr01/jxt15-djr01", "/home/Shared", m.strip() ) for m in open( cf, "r" ) ]
#        ensembler.percent = 10
#        ensembler.generate_thresholds()       
#            
#        if False:    
#            #cf="/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/spicker_run/spicker_cluster_1.list"
#            #cluster_models = [ m.strip() for m in open( cf, "r" ) ]
#            #d = "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/S_clusters/cluster_1"
#            #import glob
#            #cluster_models = [ m for m in glob.glob( os.path.join( d, "*.pdb") ) ]
#            
#            cluster_models = [ "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/4_S_00000001.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/5_S_00000003.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/1_S_00000005.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/2_S_00000006.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/4_S_00000005.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/2_S_00000003.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/4_S_00000002.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/5_S_00000002.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/4_S_00000004.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/1_S_00000003.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/3_S_00000003.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/3_S_00000005.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/2_S_00000001.pdb",
#        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_5/models/1_S_00000004.pdb" ]
#            
#            root_dir="/opt/ample-dev1/examples/toxd-example/jtest"
#            percent=50
#            ensemble_id="FOO"
#            
#            ensembler = Ensembler()
#            ensembler.maxcluster_exe = "/opt/maxcluster/maxcluster"
#            ensembler.theseus_exe = "/opt/theseus_src/theseus"
#            ensembler.generate_ensembles( cluster_models=cluster_models, root_dir=root_dir, ensemble_id=ensemble_id, percent=percent )
#            ensembles = ensembler.as_list()
#            print ensembles
#            print len(ensembles)        
    
            
if __name__ == "__main__":
    
    # This runs the ensembling starting from a pickled file containing an amopt dictionary.
    # - used when submitting the modelling jobs to a cluster
    
    if len(sys.argv) != 2 or not os.path.isfile( sys.argv[1] ):
        print "ensemble script requires the path to a pickled amopt dictionary!"
        sys.exit(1)
    
    # Get the amopt dictionary
    fpath = sys.argv[1]
    f = open( fpath, "r" )
    amoptd = cPickle.load( f )
    f.close()
    
    #if os.path.abspath(fpath) != os.path.abspath(amoptd['results_path']):
    #    print "results_path must match the path to the pickle file"
    #    sys.exit(1)
    
    # Set up logging - could append to an existing log?
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fl = logging.FileHandler("ensemble.log")
    fl.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fl.setFormatter(formatter)
    logger.addHandler(fl)
    
    # Create the ensembles
    create_ensembles( amoptd )
