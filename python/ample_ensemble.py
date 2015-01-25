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
import pdb_edit
import run_spicker
import subcluster
#from cluster_entropy import cluster - no idea where this came from...

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


class EnsembleData2(object):
    """Class to hold data about an ensemble"""

    def __init(self):

        self.name = None
        self.num_models = None
        self.num_residues = None
        self.num_atoms = None
        self.residues = None
        
        # cluster info
        self.cluster_centroid=None
        self.cluster_size=None
        
        
        self.truncation_level = None
        self.truncation_residues = None
        self.truncation_dir = None
        
        self.side_chain_treatment = None
        self.radius_threshold = None

        self.pdb = None # path to the ensemble file

        return


class Ensembler2(object):
    """Class to generate ensembles from cluster of models
    
Each method returns a list of pdbs and a corresponding list of dictionaries holding information on the returned pdbs

create_ensembles
= return a list of models

    cluster_models
    = return list of clusters - each a list of models
    
    DATA
    * cluster_method
    * cluster_number
    * models
    * centroid
    
    ######################################
    
    truncate_models
    - calculate_variances
    
    - calculate_residue_list
    
    - prune_residue_list
    
    - _truncate_models
    
    = return a list of truncation_levels - each a list of models
    
    DATA
    * truncation_level
    * variance_level
    * truncation_method
    * pruning_strategy
    
    
        ########################################
        
        - for each truncation_level:
            subcluster_models
        
        = return a list of subclusters - each a list of models
    
            ########################################
            
            for all ensembles:
               treat_side_chains


e = Ensembler(work_dir,
              cluster_method,
              cluster_exe,
              truncation_method,
              pruning_strategy,
              theseus_exe,
              subclustering_method,
              subclustering_exe,
              side_chain_treatments)
              
ensembles = e.generate_ensembles()
e.ensemble_directory
    
    """
    
    def __init__(self):
        
        # For all
        self.work_dir=None # top directory where everything gets done
        self.theseus_exe=None
        
        # clustering
        self.cluster_method="spicker" # the method for initial clustering
        self.cluster_exe=None
        self.num_clusters=1
        
        # truncation
        self.truncation_method="percent"
        self.percent_truncation=5
        self.pruning_strategy="none"
        
        # subclustering
        self.subcluster_program="maxcluster"
        self.subcluster_exe=None
        self.subclustering_method="radius"
        self.subcluster_radius_thresholds=[1,2,3]
        self.ensemble_max_models=30
        
        # Side chains
        self.side_chain_treatments=['allatom','reliable','polya']
        
        # misc
        self.logger=logging.getLogger()
        
        return
    
    def calculate_residues_percent(self,var_by_res,percent_interval):
        """Calculate the list of residues to keep if we are keeping self.percent residues under
        each truncation bin. The threshold is just the threshold of the most variable residue"""
        
        length = len(var_by_res)
        if not length > 0:
            msg = "Error reading residue variances!"
            logging.critical(msg)
            raise RuntimeError,msg
        
        # How many residues should fit in each bin
        chunk_size=int(round(float(length) * float(percent_interval)/100.0))
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
        truncation_levels=[]
        truncation_variances=[]
        truncation_residues=[]
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
                truncation_variances.append(thresh)
                truncation_levels.append(percent)
                #print "GOT PERCENT,THRESH ",percent,thresh
                #print "residues ",residues
                residues.sort()
                truncation_residues.append(residues)
                
        return truncation_levels, truncation_variances, truncation_residues
    
    def calculate_residues_thresh(self,var_by_res):
        """Txxx
        """

        # calculate the thresholds
        truncation_variances=self.generate_thresholds(var_by_res)

        # We run in reverse as that's how the original code worked
        truncation_residues=[]
        truncation_levels=[]
        lt=len(truncation_variances)
        for i, truncation_threshold in enumerate(truncation_variances):
            
            truncation_level=lt-i # as going backwards
            truncation_levels.append(truncation_level)
            
            # Get a list of the indexes of the residues to keep
            to_keep=[resSeq for resSeq,variance in var_by_res if variance <= truncation_threshold]
            truncation_residues.append( to_keep )
        
        # We went through in reverse so put things the right way around
        truncation_levels.reverse()
        truncation_variances.reverse()
        truncation_residues.reverse()
        return truncation_levels, truncation_variances, truncation_residues
        
    def calculate_variances(self,cluster_models):
        """Return a list of tuples: (resSeq,variance)"""
        
        #--------------------------------
        # get variations between pdbs
        #--------------------------------
        cmd = [ self.theseus_exe, "-a0" ] + cluster_models
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
                resSeq=int(tokens[idxResSeq])
                variance=float(tokens[idxVariance])
                assert resSeq==i,"Residue numbering doesn't match residue position in calculate_variances!"
                variances.append((resSeq,variance))
                
        return variances
    
    def cluster_models(self, models=None, cluster_method=None, num_clusters=None, cluster_exe=None):
        
        clusters=[]
        clusters_data=[]
        if cluster_method=="spicker":
            
            # Spicker Alternative for clustering
            spicker_rundir = os.path.join( self.work_dir, 'spicker')
            spickerer = run_spicker.SpickerCluster(run_dir=spicker_rundir,
                                                   spicker_exe=cluster_exe,
                                                   models=models,
                                                   num_clusters=num_clusters )
            spickerer.run_spicker()
            self.logger.info( spickerer.results_summary() )
            
            for i in range(num_clusters):
                # The models
                cluster=spickerer.results[i].pdb_list
                clusters.append(cluster)
                
                # Data on the models
                cluster_data=EnsembleData2()
                d=spickerer.results[i]
                cluster_data.cluster_centroid=d.cluster_centroid
                cluster_data.cluster_size=d.cluster_size
                clusters_data.append(cluster_data)
        else:
            raise RuntimeError,'Unrecognised clustering method: {0}'.format(cluster_method)

        return clusters, clusters_data
    
    def edit_side_chains(self,raw_ensembles,raw_ensembles_data):
        
        pdbed = pdb_edit.PDBEdit()
        ensembles=[]
        ensembles_data=[]
        for raw_ensemble, raw_ensemble_data in zip(raw_ensembles,raw_ensembles_data):
            for sct in self.side_chain_treatments:
                
                # create filename based on side chain treatment
                fpath = ample_util.filename_append(raw_ensemble,astr=sct, directory=self.ensembles_directory)
                ensemble_data=copy.copy(raw_ensemble_data)
                ensemble_data.side_chain_treatment=sct
                ensemble_data.pdb=fpath
                self.ensembles.append(fpath)
                self.ensembles_data.append(ensemble_data)
                
                # Create the files
                if sct == "allatom":
                    # For all atom just copy the file
                    shutil.copy2(raw_ensemble,fpath)
                elif sct == "reliable":
                    pdbed.reliable_sidechains(raw_ensemble,fpath)
                elif sct == "reliable":
                    pdbed.backbone(raw_ensemble,fpath)
                else:
                    raise RuntimeError,"Unrecognised side_chain_treatment: {0}\n{1}".format(sct,ensemble_data)
        
        return ensembles,ensembles_data
  
    def generate_ensembles(self, models, cluster_method, cluster_exe, work_dir):
        ensembles = []
        ensembles_data = []
        for cluster, cluster_data in self.cluster_models(models=models,
                                                         cluster_method=self.cluster_method,
                                                         num_clusters=self.num_clusters,
                                                         cluster_exe=self.cluster_exe):
            for truncated_models, truncated_models_data in self.truncate_models(cluster,
                                                                                cluster_data,
                                                                                mode=self.truncation_mode,
                                                                                percent_truncation=self.percent_truncation):
                for subcluster, subcluster_data in self.subcluster_models(truncated_models,
                                                                          truncated_models_data,
                                                                          radius_thresholds=self.subcluster_radius_thresholds,
                                                                          ensemble_max_models=self.ensemble_max_models
                                                                          ):
                    ensembles.append(subcluster)
                    ensembles_data.append(subcluster_data)
        
        # Have all-atom side-chains so now prune down the side chains
        ensembles,ensembles_data = self.edit_side_chains(ensembles,ensembles_data)
        self.ensembles=ensembles
        self.ensembles_data=ensembles_data
        
        return ensembles

    def generate_thresholds(self,var_by_res):
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

        # List of variances ordered by residue index
        var_list=[var for (resSeq,var) in var_by_res]

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
       
    def ensemble_summary(self,ensemble_data):
        
        #cluster
        #- truncation_level
        #-- subcluster
        #- variance_level
        #-- radius_level
        #--- side_chain_treatment
        
        clusters={}
        # Loop through all ensemble data objects and build up a data tree
        for e in ensemble_data:
            if e.cluster not in clusters:
                clusters[e.cluster] = {}
            if e.truncation_level not in clusters[e.cluster]:
                clusters[e.cluster] = {}
        
        return
    
    def init_from_dict(self,amoptd):
        
        self.work_dir=amoptd['work_dir'] 
        self.cluster_method=amoptd['cluster_method']
        if self.cluster_method=="spicker":
            self.cluster_exe=amoptd['spicker_exe']
            
        self.truncation_method=amoptd['truncation_method']
        self.pruning_strategy=amoptd['pruning_strategy']
        self.theseus_exe=amoptd['theseus_exe']
        self.subclustering_method="radius"
        self.subclustering_exe=None
        self.side_chain_treatments=None
        return
    
    def subcluster_models(self,
                          truncated_models,
                          truncated_models_data,
                          subcluster_program=None,
                          subcluster_exe=None,
                          radius_thresholds=None,
                          ensemble_max_models=None):
        
        ensembles=[]
        ensembles_data=[]
        
        # Use first model to get data on level
        truncation_level=truncated_models_data.truncation_level
        truncation_dir=truncated_models_data.truncation_dir
            
        # Run maxcluster to generate the distance matrix
        if subcluster_program=='maxcluster':
            clusterer = subcluster.MaxClusterer(self.subcluster_exe)
        else:
            assert False
        clusterer.generate_distance_matrix(truncated_models)

        # Loop through the radius thresholds
        num_previous_models=-1 # set to -1 so comparison always false on first pass
        for radius in radius_thresholds:

            logging.debug("subclustering files under radius: {0}".format(radius))

            # Get list of pdbs clustered according to radius threshold
            cluster_files = clusterer.cluster_by_radius(radius)
            logging.debug("Clustered {0} files".format(len(cluster_files)))
            if cluster_files < 2:
                logging.info( 'Could not cluster files using radius {0}'.format(radius))
                continue

            # For naming all files
            basename='tl{0}_r{1}'.format(truncation_level, radius)

            # Check if there are the same number of models in this ensemble as the previous one - if so
            # the ensembles will be identical and we can skip this one
            if num_previous_models == len(cluster_files):
                logging.debug( 'Number of decoys in cluster ({0}) is the same as under previous threshold so excluding cluster {1}'.format( len( cluster_files ), basename ) )
                continue
            else:
                num_previous_models = len(cluster_files)

            # Got files so create the directories
            subcluster_dir = os.path.join(truncation_dir, 'subcluster_{0}'.format(radius))
            os.mkdir(subcluster_dir)
            os.chdir(subcluster_dir)

            # Restrict cluster to max_ensemble_models
            if len(cluster_files) > ensemble_max_models:
                logging.debug("{0} files in cluster so truncating list to first {1}".format(len(cluster_files), ensemble_max_models))
                cluster_files = cluster_files[:ensemble_max_models]

            # Write out the files for reference
            file_list = "maxcluster_radius_{0}_files.list".format(radius)
            with open(file_list, "w") as f:
                for c in cluster_files:
                    f.write(c+"\n")
                f.write("\n")

            # Run theseus to generate a file containing the aligned clusters
            cmd = [ self.theseus_exe, "-r", basename, "-a0" ] + cluster_files
            retcode = ample_util.run_command( cmd, logfile=basename+"_theseus.log" )
            if retcode != 0:
                logging.debug( 'Could not create ensemble for files (models too diverse): {0}'.format( basename ) )
                continue

            # jmht - the previous Align_rosetta_fine_clusters_with_theseus routine was running theseus twice and adding the averaged
            # ensemble from the first run to the ensemble. This seems to improve the results for the TOXD test case - maybe something to
            # look at?
            #cmd = [ self.theseus_exe, "-r", basename, "-a0" ] + cluster_files + [ basename+"run1_ave.pdb" ]
            #retcode = ample_util.run_command( cmd, logfile=basename+"_theseus.log" )

            # Rename the file with the aligned files and append the path to the ensembles
            cluster_file = os.path.join(subcluster_dir, basename+'_sup.pdb')
            ensemble = os.path.join(subcluster_dir, basename+'.pdb')
            shutil.move(cluster_file, ensemble)

            # Count the number of atoms in the ensemble-only required for benchmark mode
            natoms,nresidues=pdb_edit.PDBEdit().num_atoms_and_residues(ensemble,first=True)

            # The data we've collected is the same for all pdbs in this level so just keep using the first  
            ensemble_data=copy.copy(truncated_models_data)
            ensemble_data.name = basename
            ensemble_data.num_models = len( cluster_files )
            ensemble_data.num_atoms=natoms
            ensemble_data.radius_threshold = radius
            ensemble_data.pdb = ensemble

            # Get the centroid model name from the list of files given to theseus - we can't parse
            # the pdb file as theseus truncates the filename
            ensemble_data.centroid_model=os.path.splitext( os.path.basename(cluster_files[0]) )[0]
            
            ensembles.append(ensemble)
            ensembles_data.append(ensemble_data)
        
        return ensembles,ensembles_data
    
    
    def truncate_models(self,models,models_data,truncation_method,percent_truncation):

        # Create the directories we'll be working in
        truncate_dir =  os.path.join(self.work_dir, 'truncate')
        os.mkdir( self.work_dir )
        os.chdir( self.work_dir )
        
        # Calculate variances between pdb
        var_by_res=self.calculate_variances()

        # Calculate which residues to keep under the different methods
        truncation_levels, truncation_variances, truncation_residues=None
        if truncation_method=='percent':
            truncation_levels, truncation_variances, truncation_residues=self.calculate_residues_percent(var_by_res,percent_truncation)
        elif truncation_method=='threshold':
            truncation_levels, truncation_variances, truncation_residues=self.calculate_residues_thresh(var_by_res)
        else:
            raise RuntimeError,"Unrecognised ensembling mode: {0}".format(truncation_method)

        truncated_models=[]
        truncated_models_data=[]
        for i in range(len(truncation_levels)):
            tlevel=truncation_levels[i]
            tvar=truncation_variances[i]
            tresidues=truncation_residues[i]
            trunc_dir = os.path.join( truncate_dir, 'tlevel_{0}'.format(tlevel))
            os.mkdir(trunc_dir)
            logging.info( 'truncating at: {0} in directory {1}'.format(tvar,trunc_dir))

            # list of models for this truncation level
            level_models = []
            pdbed = pdb_edit.PDBEdit()
            for infile in models:
                pdbname = os.path.basename( infile )
                pdbout = os.path.join(trunc_dir, pdbname)
                # Loop through PDB files and create new ones that only contain the residues left after truncation
                pdbed.select_residues(inpath=infile, outpath=pdbout, residues=tresidues)
                level_models.append(pdbout)
            
            # Add the model
            truncated_models.append(level_models)

            # Add the data
            model_data=copy.copy(models_data)
            model_data.truncation_level=tlevel
            model_data.truncation_variance=tvar
            model_data.truncation_residues=tresidues
            model_data.num_residues=len(tresidues)
            model_data.truncation_dir=trunc_dir
            model_data.num_truncated_models = len(models)
            
            truncated_models_data.append(model_data)
            
        return truncated_models, truncated_models_data


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
        """CONVERT TO RETURN LIST"""
        
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
            clusterer = subcluster.MaxClusterer( self.maxcluster_exe )
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


# class Test(unittest.TestCase):
# 
#     def setUp(self):
#         """
#         Get paths need to think of a sensible way to do this
#         """
# 
#         thisd =  os.path.abspath( os.path.dirname( __file__ ) )
#         paths = thisd.split( os.sep )
#         self.ample_dir = os.sep.join( paths[ : -1 ] )
#         self.theseus_exe=ample_util.find_exe('theseus')
# 
#         return
# 
#     def testThresholds(self):
#         """Test we can reproduce the original thresholds"""
# 
#         ensembler=Ensembler()
# 
#         ensembler.work_dir=os.path.join(os.getcwd(),"genthresh")
#         os.mkdir(ensembler.work_dir)
#         ensembler.theseus_exe=self.theseus_exe
#         ensembler.percent_interval=5
#         mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
#         ensembler.cluster_models=glob.glob(mdir+os.sep+"*.pdb")
#         thresholds=ensembler.generate_thresholds()
# 
#         self.assertEqual(29,len(thresholds))
#         self.assertEqual([0.030966, 0.070663, 0.084085, 0.099076, 0.185523, 0.723045,
#                           1.79753, 3.266976, 4.886417, 6.478746, 10.378687, 11.229685,
#                           16.701308, 18.823698, 22.837544, 23.741723, 25.736834, 27.761794,
#                           36.255004, 37.780944, 42.318928, 49.650732, 51.481748, 56.60685,
#                           58.865232, 60.570085, 70.749036, 81.15141, 87.045704],
#                          thresholds,
#                          "Wrong truncation thresholds"
#                          )
#         shutil.rmtree(ensembler.work_dir)
#         return
#     
#     def testCalculateResiduesThresh(self):
#         """Test we can calculate the original list of residues"""
# 
#         ensembler=Ensembler()
#         
#         # This test for percent
#         ensembler.work_dir=os.path.join(os.getcwd(),"genthresh")
#         os.mkdir(ensembler.work_dir)
#         ensembler.theseus_exe=self.theseus_exe
#         ensembler.percent_interval=5
#         mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
#         ensembler.cluster_models=glob.glob(mdir+os.sep+"*.pdb")
#         ensembler.calculate_residues_thresh()
# 
#         r=[[26, 27, 31, 32],
#            [26, 27, 28, 30, 31, 32],
#            [26, 27, 28, 29, 30, 31, 32, 33],
#            [25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
#            [23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
#            [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
#            [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
#            [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
#            [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
#            [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
#            [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41],
#            [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41],
#            [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
#            [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
#            [9, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
#            [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
#            [7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
#            [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
#            [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45],
#            [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46],
#            [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47],
#            [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 53],
#            [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 50, 53],
#            [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53],
#            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 57],
#            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 54, 56, 57],
#            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 56, 57],
#            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58],
#            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59]
#            ]
# 
#         r.reverse() # was set up for the old way
#         self.assertEqual(r,ensembler.truncation_residues)
# 
#         self.assertEqual(ensembler.truncation_levels,
#                          [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29])
#         
#         self.assertEqual(ensembler.truncation_thresholds,
#                          [87.045704, 81.15141, 70.749036, 60.570085, 58.865232, 56.60685, 51.481748, 49.650732, 42.318928, 37.780944, 
#                           36.255004, 27.761794, 25.736834, 23.741723, 22.837544, 18.823698, 16.701308, 11.229685, 10.378687, 6.478746, 
#                           4.886417, 3.266976, 1.79753, 0.723045, 0.185523, 0.099076, 0.084085, 0.070663, 0.030966])
# 
#         shutil.rmtree(ensembler.work_dir)
#         return
#     
#     def testCalculateResiduesPercent(self):
# 
#         ensembler=Ensembler()
# 
#         ensembler.work_dir=os.path.join(os.getcwd(),"genthresh")
#         os.mkdir(ensembler.work_dir)
#         ensembler.theseus_exe=self.theseus_exe
#         ensembler.percent_interval=5
#         mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
#         ensembler.cluster_models=glob.glob(mdir+os.sep+"*.pdb")
#         if not ensembler.cluster_models:
#             raise RuntimeError,"Cannot find any models in: {0}".format(mdir)
#         ensembler.calculate_residues_percent()
#         
#         r=[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59],
#            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 56, 57],
#            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 54, 57],
#            [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53],
#            [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 53],
#            [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47],
#            [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46],
#            [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
#            [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
#            [9, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
#            [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
#            [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41],
#            [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40],
#            [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
#            [21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
#            [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
#            [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
#            [26, 27, 28, 29, 30, 31, 32, 33],
#            [26, 27, 30, 31, 32],
#            [31, 32]]
#         
#         
# 
#         self.assertEqual(ensembler.truncation_levels,
#                          [100, 95, 90, 85, 80, 75, 69, 64, 59, 54, 49, 44, 39, 34, 29, 24, 19, 14, 8, 3])
#         
#         self.assertEqual(ensembler.truncation_thresholds,
#                          [75.058226, 60.570085, 56.907929, 51.481748, 43.711297, 37.780944, 31.117855,
#                           25.736834, 23.545518, 18.823698, 15.471155, 10.378687, 5.878349, 3.266976,
#                           0.901114, 0.185523, 0.08613, 0.070663, 0.030966, 0.030966])
# 
#         shutil.rmtree(ensembler.work_dir)
#         return
# 
# 
#     def XtestGenThresh(self):
#         l1=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
#         l2=[0,5,6,7,9,1,2,2,2,3,3,8,8,8,8,8,8,8,8,8,8,8,3,3,3,3,3,4,5,9]
#         l3=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5]
#         
#         e=Ensembler()
#         
#         
#         for percent in [20,50]:
#             
#             chunk_size1=int(len(l1)*float(percent)/float(100))
#             chunk_size2=int(math.ceil(len(l1)*float(percent)/float(100)))
#             
#             print "PERCENT ",percent
#             t=e._generate_thresholds(l1,chunk_size1)
#             print t
#             t=e._generate_thresholds2(l1,chunk_size2)
#             print t
#             print
#             
#             chunk_size=int(math.ceil(len(l2)*float(percent)/float(100)))
#             t=e._generate_thresholds(l2,chunk_size1)
#             print t
#             t=e._generate_thresholds2(l2,chunk_size2)
#             print t
#             print
#             
#             chunk_size=int(math.ceil(len(l3)*float(percent)/float(100)))
#             t=e._generate_thresholds(l3,chunk_size1)
#             print t
#             t=e._generate_thresholds2(l3,chunk_size2)
#             print t
# 
#         return
#     
#     def XtestTruncate(self):
#         
#         ensembler=Ensembler()
#         ensembler.work_dir=os.path.join(os.getcwd(),"genthresh")
#         os.mkdir(ensembler.work_dir)
#         ensembler.theseus_exe=self.theseus_exe
#         ensembler.percent_interval=80
#         mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
#         ensembler.cluster_models=glob.glob(mdir+os.sep+"*.pdb")
#         
#                 #print ensembler.truncate_models2()
#         return

class Test2(unittest.TestCase):

    def setUp(self):
        """
        Get paths need to think of a sensible way to do this
        """

        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        self.ample_dir = os.sep.join( paths[ : -1 ] )
        self.theseus_exe=ample_util.find_exe('theseus')
        self.spicker_exe=ample_util.find_exe('spicker')

        return

    def testThresholds(self):
        """Test we can reproduce the original thresholds"""

        ensembler=Ensembler2()

        ensembler.work_dir=os.path.join(os.getcwd(),"genthresh")
        os.mkdir(ensembler.work_dir)
        ensembler.theseus_exe=self.theseus_exe
        ensembler.percent_interval=5
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        var_by_res=ensembler.calculate_variances(cluster_models)
        thresholds=ensembler.generate_thresholds(var_by_res)
        
        self.assertEqual(30,len(thresholds),thresholds)
        self.assertEqual([0.02235, 0.041896, 0.08155, 0.085867, 0.103039, 0.185265, 0.7058, 1.772884, 3.152793,
                          4.900255, 6.206563, 10.250043, 10.772177, 16.071345, 18.292366, 22.432564, 23.265938,
                          25.348501, 27.37951, 35.877391, 37.227859, 41.759805, 49.037625, 50.992344, 56.091738,
                          58.09387, 59.917999, 70.052007, 80.235358, 86.028994],
                         thresholds,
                         "Wrong truncation thresholds"
                         )
        shutil.rmtree(ensembler.work_dir)
        return
    
    def testCalculateResiduesThresh(self):
        """Test we can calculate the original list of residues"""

        ensembler=Ensembler2()
        
        # This test for percent
        ensembler.work_dir=os.path.join(os.getcwd(),"genthresh")
        os.mkdir(ensembler.work_dir)
        ensembler.theseus_exe=self.theseus_exe
        ensembler.percent_interval=5
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        
        var_by_res=ensembler.calculate_variances(cluster_models)
        truncation_levels, truncation_variances, truncation_residues=ensembler.calculate_residues_thresh(var_by_res)
        
        self.assertEqual(truncation_levels,
                         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30])
        
        self.assertEqual(truncation_variances,
                         [86.028994, 80.235358, 70.052007, 59.917999, 58.09387, 56.091738, 50.992344, 49.037625, 41.759805, 37.227859, 35.877391,
                          27.37951, 25.348501, 23.265938, 22.432564, 18.292366, 16.071345, 10.772177, 10.250043, 6.206563, 4.900255, 3.152793,
                          1.772884, 0.7058, 0.185265, 0.103039, 0.085867, 0.08155, 0.041896, 0.02235],
                         truncation_variances)
        

        #residues=[]
        #self.assertEqual(len(r),len(truncation_residues))

        shutil.rmtree(ensembler.work_dir)
        return
    
    def testCalculateResiduesPercent(self):

        ensembler=Ensembler2()

        ensembler.work_dir=os.path.join(os.getcwd(),"genthresh")
        os.mkdir(ensembler.work_dir)
        ensembler.theseus_exe=self.theseus_exe
        ensembler.percent_interval=5
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        var_by_res=ensembler.calculate_variances(cluster_models)
        truncation_levels, truncation_variances, truncation_residues=ensembler.calculate_residues_percent(var_by_res,percent=5)


        ref=[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59],
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
            [26, 27, 28, 30, 31, 32, 33, 34],
            [26, 27, 30, 31, 32],
            [31, 32]]
        self.assertEqual(ref,truncation_residues)

        self.assertEqual(truncation_levels,
                         [100, 95, 90, 85, 80, 75, 69, 64, 59, 54, 49, 44, 39, 34, 29, 24, 19, 14, 8, 3])
        
        self.assertEqual(truncation_variances,
                         [74.272253, 59.917999, 56.260987, 50.992344, 43.268754, 37.227859, 30.767268, 25.348501,
                          23.050357, 18.292366, 15.287348, 10.250043, 5.66545, 3.152793, 0.874899, 0.185265, 0.090917,
                          0.08155, 0.039639, 0.018215])

        shutil.rmtree(ensembler.work_dir)
        return

    def testClustering(self):

        ensembler=Ensembler2()

        work_dir=os.path.join(os.getcwd(),"genthresh")
        os.mkdir(work_dir)
        ensembler.work_dir=work_dir
        ensembler.theseus_exe=self.theseus_exe
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        models=glob.glob(mdir+os.sep+"*.pdb")
        
        cluster_models,cluster_data=ensembler.cluster_models(models=models,
                                                             cluster_method='spicker',
                                                             num_clusters=3,
                                                             cluster_exe=self.spicker_exe)
        
        print cluster_models
        print cluster_data
        

        shutil.rmtree(ensembler.work_dir)
        return


def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testThresholds'))
    suite.addTest(Test('testCalculateResiduesThresh'))
    suite.addTest(Test('testCalculateResiduesPercent'))
    return suite

def testSuite2():
    suite = unittest.TestSuite()
    suite.addTest(Test2('testThresholds'))
    suite.addTest(Test2('testCalculateResiduesThresh'))
    suite.addTest(Test2('testCalculateResiduesPercent'))
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

    unittest.TextTestRunner(verbosity=2).run(testSuite2())

