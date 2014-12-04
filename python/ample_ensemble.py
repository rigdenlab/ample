'''
Created on Apr 18, 2013

@author: jmht
'''

import copy
import glob
import logging
import math
import os
import shutil
import sys
import types
import unittest

# our imports
import ample_util
import cluster_with_MAX
import pdb_edit

class EnsembleData(object):
    """Class to hold data about an ensemble"""

    def __init(self):

        self.name = None
        self.num_models = None
        self.num_residues = None
        self.num_atoms = None
        self.residues = []
        self.centroid_model=None

        self.side_chain_treatment = None
        self.radius_threshold = None
        self.truncation_threshold = None # The variance used to truncate
        self.truncation_level = None # the index of the truncation threshold

        self.pdb = None # path to the ensemble file

        return

    def __str__(self):
        """List the data attributes of this object"""
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
              isinstance(attr, types.FunctionType) ):
                me[slot] = attr

        return "{0} : {1}".format(self.__repr__(),str(me))

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
        self.percent_interval= None

        # theseus variance file
        self.variance_log = None

        # list of tupes: (resSeq,variances)
        self.var_by_res=[]

        # List of the levels
        self.truncation_levels = []
        
        # List of the thresholds we truncate at
        self.truncation_thresholds = []

        # List of directories where the truncated files are (must match truncation_levels )
        self.truncation_dirs = []

        # List of lists of the residues kept under each truncation threshold - ordered by truncation_levels
        self.truncation_residues = []

        # List of list of the truncated files under each truncation threshold - ordered by truncation_levels
        self.truncated_models = []

        # radius thresholds to cluster the models
        self.radius_thresholds = [ 1, 2, 3 ]

        # The list of ensembles before side-chain treatment
        self.truncated_ensembles = []

        # The list of final ensembles
        self.ensembles = []

        # The maximum number of models in an ensemble
        self.max_ensemble_models = None

        # Programs we use
        self.theseus_exe = None

        self.maxcluster_exe = None
        
        return

    def pdb_list(self):
        """Return the final ensembles as a list of pdb files"""
        return [ensemble.pdb for ensemble in self.ensembles]
 
    def calculate_residues_percent(self):
        """Calculate the list of residues to keep if we are keeping self.percent residues under
        each truncation bin. The threshold is just the threshold of the most variable residue"""
        
        # Calculate variances between pdb
        var_by_res=self.calculate_variances()

        length = len(var_by_res)
        if not length > 0:
            msg = "Error reading residue variances!"
            logging.critical(msg)
            raise RuntimeError,msg
        
        # How many residues should fit in each bin
        chunk_size=int(round(float(length) * float(self.percent_interval)/100.0))
        if chunk_size < 1:
            msg = "Error generating thresholds, got < 1 AA in chunk_size"
            logging.critical(msg)
            raise RuntimeError,msg
        
        nchunks=int(round(length/chunk_size))+1
        #print "chunk_size, nchunks ",chunk_size,nchunks
        
        # Get list of residue indices sorted by variance - from most variable to least
        var_by_res.sort(key=lambda x: x[1], reverse=True)
        
        #print "var_by_res ",var_by_res
        resSeq=[ x[0] for x in var_by_res ]
        
        # Get list of residues to keep under the different intevals
        self.truncation_levels=[]
        self.truncation_thresholds=[]
        self.truncation_residues=[]
        for i in range(nchunks):
            if i==0:
                residues=copy.copy(resSeq)
                percent=100
            else:
                residues=resSeq[chunk_size*i:]
                percent=int(round(float(length-(chunk_size*i))/float(length)*100))
            
            if len(residues):
                # For the threshold we take the threshold of the most variable residue
                idx=chunk_size*(i+1)-1
                if idx > length-1: # Need to make sure we have a full final chunk
                    idx=length-1
                thresh=var_by_res[idx][1]
                self.truncation_thresholds.append(thresh)
                self.truncation_levels.append(percent)
                #print "GOT PERCENT,THRESH ",percent,thresh
                #print "residues ",residues
                residues.sort()
                self.truncation_residues.append(residues)
                
        return
    
    def calculate_residues_thresh(self):
        """Txxx
        """

        # calculate the thresholds
        self.truncation_thresholds=self.generate_thresholds()

        # We run in reverse as that's how the original code worked
        self.truncation_residues=[]
        self.truncation_levels=[]
        lt=len(self.truncation_thresholds)
        for i, truncation_threshold in enumerate(self.truncation_thresholds):
            
            truncation_level=lt-i # as going backwards
            self.truncation_levels.append(truncation_level)
            
            # Get a list of the indexes of the residues to keep
            to_keep=[resSeq for resSeq,variance in self.var_by_res if variance <= truncation_threshold]
            self.truncation_residues.append( to_keep )
        
        # We went through in reverse so put things the right way around
        self.truncation_levels.reverse()
        self.truncation_thresholds.reverse()
        self.truncation_residues.reverse()
        return
        
    def calculate_variances(self):
        
        #--------------------------------
        # get variations between pdbs
        #--------------------------------
        cmd = [ self.theseus_exe, "-a0" ] + self.cluster_models
        retcode = ample_util.run_command(cmd,
                                         logfile=os.path.join(self.work_dir,"theseus.log"),
                                         directory=self.work_dir)
        if retcode != 0:
            msg = "non-zero return code for theseus in generate_thresholds!"
            logging.critical( msg )
            raise RuntimeError, msg

        variances=[]
        variance_log=os.path.join(self.work_dir,'theseus_variances.txt')
        with open(variance_log) as f:
            for i, line in enumerate(f):
                # Skip header
                if i==0: continue

                line=line.strip()

                # Skip blank lines
                if not line: continue

                #print line
                tokens=line.split()
                # Different versions of theseus may have a RES card first, so need to check
                if tokens[0]=="RES":
                    idxResSeq=3
                    idxVariance=4
                else:
                    idxResSeq=2
                    idxVariance=3

                variances.append( ( int(tokens[idxResSeq]),float(tokens[idxVariance]) ) )

        return variances

    def edit_sidechains(self):
        """Take the ensembles and give them the 3 sidechain treatments"""

        pdbed = pdb_edit.PDBEdit()


        for trunc_ensemble in self.truncated_ensembles:

            # 3 different side-chain treatments

            # 1. all_atom
            ensemble = copy.deepcopy( trunc_ensemble )
            ensemble.side_chain_treatment = "allatom"
            ensemble.name = "{0}_{1}".format( ensemble.side_chain_treatment, trunc_ensemble.name )
            ensemble.pdb = os.path.join( self.ensemble_dir, ensemble.name+".pdb" )
            # For all atom just copy the file
            shutil.copy2( trunc_ensemble.pdb, ensemble.pdb )
            self.ensembles.append( ensemble )

            # 2. Reliable side chains
            ensemble = copy.deepcopy( trunc_ensemble )
            ensemble.side_chain_treatment = "reliable"
            #ensemble.side_chain_treatment = "SCWRL_reliable_sidechains"
            ensemble.name = "{0}_{1}".format( ensemble.side_chain_treatment, trunc_ensemble.name )
            ensemble.pdb = os.path.join( self.ensemble_dir, ensemble.name+".pdb" )
            # Do the edit
            pdbed.reliable_sidechains( trunc_ensemble.pdb, ensemble.pdb )
            self.ensembles.append( ensemble )

            # 3. backbone
            ensemble = copy.deepcopy( trunc_ensemble )
            ensemble.side_chain_treatment = "polya"
            ensemble.name = "{0}_{1}".format( ensemble.side_chain_treatment, trunc_ensemble.name )
            ensemble.pdb = os.path.join( self.ensemble_dir, ensemble.name+".pdb" )
            # Do the edit
            pdbed.backbone( trunc_ensemble.pdb, ensemble.pdb )
            self.ensembles.append( ensemble )

        return

    def generate_ensembles(self, cluster_models=None, root_dir=None, ensemble_id=None, percent=None, mode='threshold' ):
        """Generate the ensembles for this list of cluster models"""

        assert self.maxcluster_exe and self.theseus_exe, "Must set the path to maxcluster and theseus!"

        self.cluster_models = cluster_models
        self.percent_interval = percent

        # Create the directories we'll be working in
        self.work_dir =  os.path.join( root_dir, 'fine_cluster_{0}'.format( ensemble_id ) )
        os.mkdir( self.work_dir )
        os.chdir( self.work_dir )

        self.ensemble_dir = os.path.join( root_dir, 'ensembles_{0}'.format( ensemble_id ) )
        os.mkdir( self.ensemble_dir )

        # Calculate which residues to keep under the different methods
        if mode=='threshold':
            self.calculate_residues_thresh()
        elif mode=='percent':
            self.calculate_residues_percent()
        else:
            raise RuntimeError,"Unrecognised ensembling mode: {0}".format(mode)
        
        # Truncate the models
        self.truncate_models()

        # Generate the ensembles for the different truncation levels
        self.make_truncated_ensembles()

        # Undertake the 3 side-chain treatments
        self.edit_sidechains()

        return

    def generate_thresholds(self):
        """
        This is the original method developed by Jaclyn and used in all work until November 2014 (including the coiled-coil paper)
        
        Calculate the residue variance thresholds that will keep self.percent_interval residues for each truncation level
        """
        #--------------------------------
        # choose threshold type
        #-------------------------------
        FIXED_INTERVALS=False
        if FIXED_INTERVALS:
            self.thresholds = [ 1, 1.5, 2 , 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8 ]
            logging.debug("Got {0} thresholds: {1}".format( len(self.thresholds), self.thresholds ))
            return

        # Calculate variances between pdb
        self.var_by_res=self.calculate_variances()

        # List of variances ordered by residue index
        var_list=[var for (resSeq,var) in self.var_by_res]

        length = len(var_list)
        if length == 0:
            msg = "Error generating thresholds, got len: {0}".format(length)
            logging.critical(msg)
            raise RuntimeError,msg

        # How many residues should fit in each bin
        # NB - Should round up not down with int!
        chunk_size=int( ( float(length)/100 ) *float(self.percent_interval) )
        if chunk_size < 1:
            msg = "Error generating thresholds, got < 1 AA in chunk_size"
            logging.critical(msg)
            raise RuntimeError,msg

        ## try to find intervals for truncation
        truncation_thresholds=self._generate_thresholds(var_list, chunk_size)
        
        # Jens' new untested method
        #truncation_thresholds=self._generate_thresholds2(var_list, chunk_size)
        
        logging.debug("Got {0} thresholds: {1}".format( len(truncation_thresholds), truncation_thresholds ))

        return truncation_thresholds

    def _generate_thresholds(self,values,chunk_size):
        """Jaclyn's threshold method
        """
        try_list=copy.deepcopy(values)
        try_list.sort()
        #print "try_list ",try_list

        # print list(chunks(try_list, int(chunk_size)))
        # For chunking list
        def chunks(a_list, chunk_size):
            for i in xrange(0, len(a_list), chunk_size):
                yield a_list[i:i+chunk_size ]

        thresholds=[]
        for x in list( chunks(try_list, chunk_size ) ):
            #print x, x[-1]
            # For some cases, multiple residues share the same variance so we don't create a separate thereshold
            if x[-1] not in thresholds:
                thresholds.append(x[-1])
                
        return thresholds

    def _generate_thresholds2(self,values,chunk_size):
        """
        This is Jens's update to Jaclyn's method that groups the residues by variances so that we split
        them by variance, and try and fit chunk_size in each bin. Previously we tried to split by variance but didn't
        group the residues by variance, so the same variance bin could cover multiple residue groups. 
        """

        # Create tuple mapping values to counts
        data=[(i,values.count(i)) for i in sorted(set(values),reverse=True)]

        thresholds=[]
        counts=[]
        first=True
        for variance,count in data:
            if first or counts[-1] + count > chunk_size:
                thresholds.append(variance)
                counts.append(count)
                if first: first=False
            else:
                #thresholds[-1]=variance
                counts[-1]+=count

        thresholds.sort()
        return thresholds

    def make_truncated_ensembles( self  ):
        """
        Loop through each truncation level and use maxcluster to cluster the models
        according to the three radius thresholds in self.radius_thresholds

        If there are more than 2 clusters returned by maxcluster use theseus to align them

        Generate an ensemble for

        """

        # Loop through all truncation levels
        for i,truncation_level in enumerate(self.truncation_levels):
            
            # change to the truncation directory
            truncation_dir = self.truncation_dirs[i]
            os.chdir( truncation_dir )

            # Run maxcluster to generate the distance matrix
            clusterer = cluster_with_MAX.MaxClusterer( self.maxcluster_exe )
            clusterer.generate_distance_matrix( self.truncated_models[i] )

            # Loop through the radius thresholds
            num_previous_models=-1 # set to -1 so comparison always false on first pass
            for radius in self.radius_thresholds:

                logging.debug("Clustering files under radius: {0}".format( radius ) )

                # Get list of pdbs clustered according to radius threshold
                cluster_files = clusterer.cluster_by_radius( radius )
                logging.debug("Maxcluster clustered {0} files".format ( len( cluster_files ) ) )
                if cluster_files < 2:
                    logging.info( 'Could not cluster files using radius {0}'.format( radius ) )
                    continue

                # For naming all files
                #basename='trunc_{0}_rad_{1}'.format( truncation_threshold, radius )
                basename='tl{0}_r{1}'.format( truncation_level, radius )

                # Check if there are the same number of models in this ensemble as the previous one - if so
                # the ensembles will be identical and we can skip this one
                if num_previous_models == len( cluster_files ):
                    logging.debug( 'Number of decoys in cluster ({0}) is the same as under previous threshold so excluding cluster {1}'.format( len( cluster_files ), basename ) )
                    continue
                else:
                    num_previous_models = len( cluster_files )

                # Got files so create the directories
                ensemble_dir = os.path.join( truncation_dir, 'fine_clusters_{0}_ensemble'.format(radius) )
                os.mkdir( ensemble_dir )
                os.chdir( ensemble_dir )

                # Restrict cluster to self.max_ensemble_models
                assert self.max_ensemble_models
                if len( cluster_files ) > self.max_ensemble_models:
                    logging.debug("{0} files in cluster so truncating list to first {1}".format( len( cluster_files ), self.max_ensemble_models) )
                    cluster_files = cluster_files[ :self.max_ensemble_models ]

                # Write out the files for reference
                file_list = "maxcluster_radius_{0}_files.list".format( radius )
                with open(file_list, "w") as f:
                    for c in cluster_files:
                        f.write( c+"\n")
                    f.write("\n")

                # Run theseus to generate a file containing the aligned clusters
                cmd = [ self.theseus_exe, "-r", basename, "-a0" ] + cluster_files
                retcode = ample_util.run_command( cmd, logfile=basename+"_theseus.log" )
                if retcode != 0:
                    logging.debug( 'Could not create ensemble for files (models too diverse): {0}'.format( basename ) )
                    continue

                #if not os.path.exists( cluster_file ):
                #    logging.info( 'Could not create ensemble for files (models too diverse): {0}'.format( ensemble_file ) )
                #    continue

                # jmht - the previous Align_rosetta_fine_clusters_with_theseus routine was running theseus twice and adding the averaged
                # ensemble from the first run to the ensemble. This seems to improve the results for the TOXD test case - maybe something to
                # look at?
                #cmd = [ self.theseus_exe, "-r", basename, "-a0" ] + cluster_files + [ basename+"run1_ave.pdb" ]
                #retcode = ample_util.run_command( cmd, logfile=basename+"_theseus.log" )

                # Rename the file with the aligned files and append the path to the ensembles
                cluster_file = os.path.join(  ensemble_dir, basename+'_sup.pdb' )
                ensemble_file = os.path.join( ensemble_dir, basename+'.pdb' )
                shutil.move( cluster_file, ensemble_file )

                # Count the number of atoms in the ensemble-only required for benchmark mode
                natoms,nresidues=pdb_edit.PDBEdit().num_atoms_and_residues(ensemble_file,first=True)

                # Generate the ensemble
                ensemble = EnsembleData()
                ensemble.name = basename
                ensemble.num_models = len( cluster_files )
                ensemble.residues = self.truncation_residues[i]
                ensemble.num_residues = len(ensemble.residues)
                assert ensemble.num_residues==nresidues
                ensemble.num_atoms=natoms
                ensemble.truncation_level = truncation_level
                ensemble.truncation_threshold = self.truncation_thresholds[i]
                ensemble.num_truncated_models = len( self.truncated_models[i] )
                ensemble.radius_threshold = radius
                ensemble.pdb = ensemble_file

                # Get the centroid model name from the list of files given to theseus - we can't parse
                # the pdb file as theseus truncates the filename
                ensemble.centroid_model=os.path.splitext( os.path.basename(cluster_files[0]) )[0]

                self.truncated_ensembles.append( ensemble )

        return
    
    def truncate_models(self):
        self.truncation_dirs=[]
        self.truncation_models=[]
        for i in range(len(self.truncation_levels)):
            trunc_dir = os.path.join( self.work_dir, 'trunc_files_{0}'.format(self.truncation_levels[i]))
            os.mkdir(trunc_dir)
            logging.info( 'truncating at: {0} in directory {1}'.format(self.truncation_thresholds[i],trunc_dir))
            self.truncation_dirs.append(trunc_dir)

            # list of models for this truncation level
            level_models = []
            pdbed = pdb_edit.PDBEdit()
            for infile in self.cluster_models:
                pdbname = os.path.basename( infile )
                pdbout = os.path.join( trunc_dir, pdbname )
                # Loop through PDB files and create new ones that only contain the residues left after truncation
                pdbed.select_residues( inpath=infile, outpath=pdbout, residues=self.truncation_residues[i] )
                level_models.append( pdbout )

            self.truncated_models.append( level_models )
        return

class Test(unittest.TestCase):

    def setUp(self):
        """
        Get paths need to think of a sensible way to do this
        """

        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        self.ample_dir = os.sep.join( paths[ : -1 ] )
        self.theseus_exe=ample_util.find_exe('theseus')

        return

    def testThresholds(self):
        """Test we can reproduce the original thresholds"""

        ensembler=Ensembler()

        ensembler.work_dir=os.path.join(os.getcwd(),"genthresh")
        os.mkdir(ensembler.work_dir)
        ensembler.theseus_exe=self.theseus_exe
        ensembler.percent_interval=5
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        ensembler.cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        thresholds=ensembler.generate_thresholds()

        self.assertEqual(29,len(thresholds))
        self.assertEqual([0.030966, 0.070663, 0.084085, 0.099076, 0.185523, 0.723045,
                          1.79753, 3.266976, 4.886417, 6.478746, 10.378687, 11.229685,
                          16.701308, 18.823698, 22.837544, 23.741723, 25.736834, 27.761794,
                          36.255004, 37.780944, 42.318928, 49.650732, 51.481748, 56.60685,
                          58.865232, 60.570085, 70.749036, 81.15141, 87.045704],
                         thresholds,
                         "Wrong truncation thresholds"
                         )
        shutil.rmtree(ensembler.work_dir)
        return
    
    def testCalculateResiduesThresh(self):
        """Test we can calculate the original list of residues"""

        ensembler=Ensembler()
        
        # This test for percent
        ensembler.work_dir=os.path.join(os.getcwd(),"genthresh")
        os.mkdir(ensembler.work_dir)
        ensembler.theseus_exe=self.theseus_exe
        ensembler.percent_interval=5
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        ensembler.cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        ensembler.calculate_residues_thresh()

        r=[[26, 27, 31, 32],
           [26, 27, 28, 30, 31, 32],
           [26, 27, 28, 29, 30, 31, 32, 33],
           [25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
           [23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
           [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
           [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
           [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
           [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
           [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
           [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41],
           [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41],
           [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
           [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
           [9, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
           [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
           [7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
           [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
           [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45],
           [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46],
           [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47],
           [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 53],
           [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 50, 53],
           [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53],
           [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 57],
           [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 54, 56, 57],
           [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 56, 57],
           [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58],
           [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59]
           ]

        r.reverse() # was set up for the old way
        self.assertEqual(r,ensembler.truncation_residues)

        self.assertEqual(ensembler.truncation_levels,
                         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29])
        
        self.assertEqual(ensembler.truncation_thresholds,
                         [87.045704, 81.15141, 70.749036, 60.570085, 58.865232, 56.60685, 51.481748, 49.650732, 42.318928, 37.780944, 
                          36.255004, 27.761794, 25.736834, 23.741723, 22.837544, 18.823698, 16.701308, 11.229685, 10.378687, 6.478746, 
                          4.886417, 3.266976, 1.79753, 0.723045, 0.185523, 0.099076, 0.084085, 0.070663, 0.030966])

        shutil.rmtree(ensembler.work_dir)
        return
    
    def testCalculateResiduesPercent(self):

        ensembler=Ensembler()

        ensembler.work_dir=os.path.join(os.getcwd(),"genthresh")
        os.mkdir(ensembler.work_dir)
        ensembler.theseus_exe=self.theseus_exe
        ensembler.percent_interval=5
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        ensembler.cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        if not ensembler.cluster_models:
            raise RuntimeError,"Cannot find any models in: {0}".format(mdir)
        ensembler.calculate_residues_percent()
        
        r=[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59],
           [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 56, 57],
           [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 54, 57],
           [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53],
           [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 53],
           [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47],
           [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46],
           [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
           [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
           [9, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
           [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
           [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41],
           [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40],
           [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
           [21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
           [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
           [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
           [26, 27, 28, 29, 30, 31, 32, 33],
           [26, 27, 30, 31, 32],
           [31, 32]]
        
        

        self.assertEqual(ensembler.truncation_levels,
                         [100, 95, 90, 85, 80, 75, 69, 64, 59, 54, 49, 44, 39, 34, 29, 24, 19, 14, 8, 3])
        
        self.assertEqual(ensembler.truncation_thresholds,
                         [75.058226, 60.570085, 56.907929, 51.481748, 43.711297, 37.780944, 31.117855,
                          25.736834, 23.545518, 18.823698, 15.471155, 10.378687, 5.878349, 3.266976,
                          0.901114, 0.185523, 0.08613, 0.070663, 0.030966, 0.030966])

        shutil.rmtree(ensembler.work_dir)
        return


    def XtestGenThresh(self):
        l1=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
        l2=[0,5,6,7,9,1,2,2,2,3,3,8,8,8,8,8,8,8,8,8,8,8,3,3,3,3,3,4,5,9]
        l3=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5]
        
        e=Ensembler()
        
        
        for percent in [20,50]:
            
            chunk_size1=int(len(l1)*float(percent)/float(100))
            chunk_size2=int(math.ceil(len(l1)*float(percent)/float(100)))
            
            print "PERCENT ",percent
            t=e._generate_thresholds(l1,chunk_size1)
            print t
            t=e._generate_thresholds2(l1,chunk_size2)
            print t
            print
            
            chunk_size=int(math.ceil(len(l2)*float(percent)/float(100)))
            t=e._generate_thresholds(l2,chunk_size1)
            print t
            t=e._generate_thresholds2(l2,chunk_size2)
            print t
            print
            
            chunk_size=int(math.ceil(len(l3)*float(percent)/float(100)))
            t=e._generate_thresholds(l3,chunk_size1)
            print t
            t=e._generate_thresholds2(l3,chunk_size2)
            print t

        return
    
    def XtestTruncate(self):
        
        ensembler=Ensembler()
        ensembler.work_dir=os.path.join(os.getcwd(),"genthresh")
        os.mkdir(ensembler.work_dir)
        ensembler.theseus_exe=self.theseus_exe
        ensembler.percent_interval=80
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        ensembler.cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        
                #print ensembler.truncate_models2()
        return

def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testThresholds'))
    suite.addTest(Test('testCalculateResiduesThresh'))
    suite.addTest(Test('testCalculateResiduesPercent'))
    return suite
    
#
# Run unit tests
if __name__ == "__main__":

    if False:
        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        formatter = logging.Formatter('%(message)s')
        ch.setFormatter(formatter)
        root.addHandler(ch)

    unittest.TextTestRunner(verbosity=2).run(testSuite())

