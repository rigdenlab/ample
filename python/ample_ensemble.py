'''
Created on Apr 18, 2013

@author: jmht
'''

import copy
import glob
import logging
import os
import random
import shutil
import sys
import unittest

# our imports
import ample_sequence
import ample_util
import fast_protein_cluster
import pdb_edit
import spicker
import subcluster
import theseus

POLYALA='polyAla'
RELIABLE='reliable'
ALLATOM='allatom'
SIDE_CHAIN_TREATMENTS=[POLYALA,RELIABLE,ALLATOM]

def model_core_from_alignment(models,alignment_file,work_dir=None):
    
    if not work_dir: work_dir = os.path.join(os.getcwd(),'core_models')
    if not os.path.isdir(work_dir): os.mkdir(work_dir)
    
    # Read in alignment to get
    align_seq = ample_sequence.Sequence(fasta=alignment_file)
    
    # Check all alignments the same length
    
    # Get pdb names from alignment headers
    seq_names = [ h[1:].strip() for h in align_seq.headers ]
    
    GAP = '-'
    # Get array specifying which positions are core
    core = [ all([ x != GAP for x in t ]) for t in zip(*align_seq.sequences) ]
    
    # For each sequence, get a list of which positions are core
    core_positions = []
    for seq in align_seq.sequences:
        p = []
        count = 0
        for i, pos in enumerate(seq):
            if pos != GAP:
                if core[i]: p.append(count)
                count += 1
        core_positions.append(p)
        
    # Should check lengths of sequences match the length of the aa in the pdbs
        
    # Create dict mapping seq_names to core positions
    core_dict = dict( (s,core_positions[i]) for i, s in enumerate(seq_names) )
    
    # Cut the models down to core
    core_models = []
    for m in models:
        name = os.path.basename(m)
        pdbout = ample_util.filename_append(m, astr='core', directory=work_dir)
        pdb_edit.select_residues(m, pdbout, tokeep_idx=core_dict[name])
        core_models.append(pdbout)
        
    return core_models

def split_sequence(length,percent_interval):
    """split a sequence of length into chunks each separated by percent_interval"""
    
    # How many residues should fit in each bin
    chunk_size=int(round(float(length) * float(percent_interval)/100.0))
    if chunk_size < 1:
        msg = "Error splitting sequence, got chunk_size < 1: {0} : {1}".format(length,percent_interval)
        raise RuntimeError,msg
    
    nchunks=int(round(length/chunk_size))+1
    
    # Get list of residues to keep under the different intevals
    levels=[]
    indices=[]
    for i in range(nchunks):
        if i==0:
            start_stop=(0,length-1)
            percent=100
        else:
            start_stop=(chunk_size*i,length-1)
            percent=int(round(float(length-(chunk_size*i))/float(length)*100))
        levels.append(percent)
        indices.append(start_stop)
        
    return chunk_size,indices, percent

class Ensembler(object):
    """Class to generate ensembles from ab inito models (all models must have same sequence)
    
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
        self.truncation_pruning="none"
        self.percent_truncation=5
        self.pruning_strategy="none"
        
        # subclustering
        self.subcluster_method='ORIGINAL'
        self.subcluster_program="maxcluster"
        self.subcluster_exe=None
        self.subclustering_method="radius"
        self.subcluster_radius_thresholds=[1,2,3]
        self.ensemble_max_models=30
        
        # Side chains
        self.side_chain_treatments=SIDE_CHAIN_TREATMENTS
        
        # ensembles
        self.ensembles_directory=None
        self.ensembles=None
        self.ensembles_data=None
        
        # misc
        self.logger=logging.getLogger()
        
        return
    
    def align_models(self, models, basename=None, work_dir=None):
        run_theseus = theseus.Theseus(work_dir=work_dir, theseus_exe=self.theseus_exe)
        try:
            run_theseus.align_models(models, basename=basename)
        except Exception,e:
            self.logger.critical("Error running theseus: {0}".format(e))
            return False
        return run_theseus.superposed_models
    
    def _calculate_residues_percent(self,var_by_res,percent_interval):
        """Calculate the list of residues to keep if we are keeping self.percent residues under
        each truncation bin. The threshold is just the threshold of the most variable residue"""
         
        length = len(var_by_res)
        if not length > 0:
            msg = "Error reading residue variances!"
            self.logger.critical(msg)
            raise RuntimeError,msg
         
        # How many residues should fit in each bin
        chunk_size=int(round(float(length) * float(percent_interval)/100.0))
        if chunk_size < 1:
            msg = "Error generating thresholds, got < 1 AA in chunk_size"
            self.logger.critical(msg)
            raise RuntimeError,msg
         
        nchunks=int(round(length/chunk_size))+1
        #print "chunk_size, nchunks ",chunk_size,nchunks
         
        # Get list of residue indices sorted by variance - from most variable to least
        var_by_res.sort(key=lambda x: x[2], reverse=True)
         
        #print "var_by_res ",var_by_res
        idxs_all = [ x[0] for x in var_by_res ]
        resseq_all = [ x[1] for x in var_by_res ]
         
        # Get list of residues to keep under the different intevals
        truncation_levels=[]
        truncation_variances=[]
        truncation_residues=[]
        truncation_residue_idxs=[]
        MIN_RESIDUES=2 # we need at least 3 residues for theseus to work
        for i in range(nchunks):
            if i==0:
                residues=copy.copy(resseq_all)
                idxs=copy.copy(idxs_all)
                percent=100
            else:
                residues=resseq_all[chunk_size*i:]
                idxs=idxs_all[chunk_size*i:]
                percent=int(round(float(length-(chunk_size*i))/float(length)*100))
             
            if len(residues) > MIN_RESIDUES: 
                # For the threshold we take the threshold of the most variable residue
                idx=chunk_size*(i+1)-1
                if idx > length-1: # Need to make sure we have a full final chunk
                    idx=length-1
                thresh=var_by_res[idx][2]
                truncation_variances.append(thresh)
                truncation_levels.append(percent)
                #print "GOT PERCENT,THRESH ",percent,thresh
                #print "residues ",residues
                residues.sort()
                truncation_residues.append(residues)
                truncation_residue_idxs.append(idxs)
                 
        return truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs
    
    def _calculate_residues_percentX(self,var_by_res,percent_interval):
        """Calculate the list of residues to keep if we are keeping self.percent residues under
        each truncation bin. The threshold is just the threshold of the most variable residue"""
        
        length = len(var_by_res)
        if not length > 0:
            msg = "Error reading residue variances!"
            self.logger.critical(msg)
            raise RuntimeError,msg
        
        chunk_size, indices, percents = split_sequence(length,percent_interval)
        
        # Get list of residue indices sorted by variance - from most variable to least
        var_by_res.sort(key=lambda x: x[1], reverse=True)
        
        #print "var_by_res ",var_by_res
        resSeq=[ x[0] for x in var_by_res ]
        
        # Get list of residues to keep under the different intevals
        truncation_levels=[]
        truncation_variances=[]
        truncation_residues=[]
        MIN_RESIDUES=2 # we need at least 3 residues for theseus to work
        for i in range(len(indices)):
            
            start,stop = indices[i]
            percent = percents[i]
            residues=resSeq[start,stop]
            
            if len(residues) > MIN_RESIDUES: 
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
    
    def _calculate_residues_thresh(self,var_by_res,percent_interval):
        """Txxx
        """

        # calculate the thresholds
        truncation_variances=self.generate_thresholds(var_by_res,percent_interval)

        # We run in reverse as that's how the original code worked
        truncation_residues=[]
        truncation_residue_idxs=[]
        truncation_levels=[]
        lt=len(truncation_variances)
        for i, truncation_threshold in enumerate(truncation_variances):
            
            truncation_level=lt-i # as going backwards
            truncation_levels.append(truncation_level)
            
            # Get a list of the indexes of the residues to keep
            to_keep = [resSeq for idx, resSeq, variance in var_by_res if variance <= truncation_threshold]
            to_keep_idxs = [idx for idx, resSeq, variance in var_by_res if variance <= truncation_threshold]
            truncation_residues.append(to_keep)
            truncation_residue_idxs.append(to_keep_idxs)
        
        # We went through in reverse so put things the right way around
        truncation_levels.reverse()
        truncation_variances.reverse()
        truncation_residues.reverse()
        truncation_residue_idxs.reverse()
        return truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs
        
    def calculate_variances(self,cluster_models):
        """Return a list of tuples: (resSeq,variance)"""
        
        #--------------------------------
        # get variations between pdbs
        #--------------------------------
        if not os.path.exists(self.theseus_exe) and os.access(self.theseus_exe, os.X_OK):
            raise RuntimeError,"Cannot find theseus_exe: {0}".format(self.theseus_exe) 
        
        cmd = [ self.theseus_exe, "-a0" ] + cluster_models
        logfile=os.path.join(self.work_dir,"theseus.log")
        retcode = ample_util.run_command(cmd,
                                         logfile=logfile,
                                         directory=self.work_dir)
        if retcode != 0:
            msg = "non-zero return code for theseus in generate_thresholds!\n See log: {0}".format(logfile)
            self.logger.critical(msg)
            raise RuntimeError, msg

        variances=[]
        variance_log=os.path.join(self.work_dir,'theseus_variances.txt')
        with open(variance_log) as f:
            for i, line in enumerate(f):
                # Skip header
                if i==0: continue

                line=line.strip()
                if not line: continue # Skip blank lines

                #print line
                tokens=line.split()
                # Different versions of theseus may have a RES card first, so need to check
                if tokens[0]=="RES":
                    idxidx=1
                    idxResSeq=3
                    idxVariance=4
                else:
                    idxidx=0
                    idxResSeq=2
                    idxVariance=3
                idx = int(tokens[idxidx])
                assert idx == i,"Index and atom lines don't match! {0} : {1}".format(idx,i) # paranoid check
                # Theseus counts from 1, we count from 0
                idx -= 1
                resSeq = int(tokens[idxResSeq])
                variance = float(tokens[idxVariance])
                variances.append((idx,resSeq,variance))
                
        return variances
    
    def cluster_models(self,
                       models=None,
                       cluster_method=None,
                       num_clusters=None,
                       cluster_exe=None,
                       nproc=1,
                       max_cluster_size=200
                       ):
        clusters=[]
        clusters_data=[]
        if cluster_method=="spicker":
            # Spicker Alternative for clustering
            self.logger.info('* Running SPICKER to cluster models *')
            spicker_rundir = os.path.join(self.work_dir, 'spicker')
            spickerer = spicker.Spickerer(spicker_exe=cluster_exe)
            spickerer.cluster(models,
                              num_clusters=num_clusters,
                              max_cluster_size=max_cluster_size,
                              run_dir=spicker_rundir)
            self.logger.debug(spickerer.results_summary())
            
            for i in range(num_clusters):
                # The models
                cluster=spickerer.results[i].pdb_list
                clusters.append(cluster)
                # Data on the models
                cluster_data=self.create_dict()
                d=spickerer.results[i]
                cluster_data['cluster_num']=i+1
                cluster_data['cluster_centroid']=d.cluster_centroid
                cluster_data['cluster_num_models']=d.cluster_size
                cluster_data['cluster_method']=cluster_method
                cluster_data['num_clusters']=num_clusters
                clusters_data.append(cluster_data)
        
        elif cluster_method=="fast_protein_cluster":
            fpc=fast_protein_cluster.FPC()
            SCORE_TYPE='rmsd'
            CLUSTER_METHOD='kmeans'
            self.logger.info('Running fast_protein_cluster with: score_type: {0} cluster_method: {1}'.format(SCORE_TYPE,
                                                                                                             CLUSTER_METHOD))
            fpc_rundir = os.path.join(self.work_dir, 'fast_protein_cluster')
            self.logger.info('fast_protein_cluster running in directory: {0}'.format(fpc_rundir))
            clusters, clusters_data = fpc.cluster(models=models,
                                                num_clusters=num_clusters,
                                                score_type=SCORE_TYPE,
                                                cluster_method=CLUSTER_METHOD,
                                                work_dir=fpc_rundir,
                                                fpc_exe=cluster_exe,
                                                nproc=nproc,
                                                max_cluster_size=max_cluster_size)
        else:
            raise RuntimeError,'Unrecognised clustering method: {0}'.format(cluster_method)

        return clusters, clusters_data
    
    def create_dict(self):
        """Create an empty dictionary
        Not strictly necessary but it's a place to remember what we capture
        """
        d={}
        d['cluster_method'] = None
        d['num_clusters'] = None
        d['cluster_num'] = None
        d['cluster_centroid'] = None
        d['cluster_num_models'] = None
        
        # truncation info
        d['truncation_level'] = None
        d['percent_truncation'] = None
        d['truncation_method'] = None
        d['truncation_residues'] = None
        d['truncation_dir'] = None
        d['truncation_variance'] = None
        d['truncation_num_residues'] = None

        # subclustering info
        d['subcluster_num_models'] = None
        d['subcluster_radius_threshold'] = None
        d['subcluster_centroid_model'] = None
    
        # ensemble info
        d['name'] = None
        d['side_chain_treatment'] = None
        d['ensemble_num_atoms'] = None
        d['ensemble_pdb'] = None # path to the ensemble file
        
        return d
    
    def edit_side_chains(self,raw_ensemble,raw_ensemble_data,ensembles_directory):
        
        assert os.path.isdir(ensembles_directory),"Cannot find ensembles directory: {0}".format(ensembles_directory)
        ensembles=[]
        ensembles_data=[]
        for sct in self.side_chain_treatments:
            
            # create filename based on side chain treatment
            fpath = ample_util.filename_append(raw_ensemble,astr=sct, directory=ensembles_directory)
            
            # Create the files
            if sct == ALLATOM:
                # For all atom just copy the file
                shutil.copy2(raw_ensemble,fpath)
            elif sct == RELIABLE:
                pdb_edit.reliable_sidechains(raw_ensemble,fpath)
            elif sct == POLYALA:
                pdb_edit.backbone(raw_ensemble,fpath)
            else:
                raise RuntimeError,"Unrecognised side_chain_treatment: {0}".format(sct)
            
            # Count the number of atoms in the ensemble-only required for benchmark mode
            natoms,nresidues=pdb_edit.num_atoms_and_residues(fpath,first=True)
            
            # Process ensemble data
            ensemble_data=copy.copy(raw_ensemble_data)
            ensemble_data['side_chain_treatment']=sct
            ensemble_data['name']='c{0}_tl{1}_r{2}_{3}'.format(ensemble_data['cluster_num'],
                                                               ensemble_data['truncation_level'],
                                                               ensemble_data['subcluster_radius_threshold'],
                                                               sct)
            ensemble_data['ensemble_pdb']=fpath
            ensemble_data['ensemble_num_atoms']=natoms
            # check
            assert ensemble_data['truncation_num_residues']==nresidues,"Unmatching number of residues!"
            
            ensembles.append(fpath)
            ensembles_data.append(ensemble_data)
                
        return ensembles,ensembles_data
  
    def generate_ensembles(self,models,
                           cluster_method=None,
                           cluster_exe=None,
                           num_clusters=None,
                           percent_truncation=None,
                           truncation_method=None,
                           truncation_pruning=None,
                           ensembles_directory=None,
                           work_dir=None,
                           nproc=None):
        
        # Work dir set each time
        if not work_dir:
            raise RuntimeError,"Need to set work_dir!"
        self.work_dir=work_dir
        
        if not cluster_method:
            cluster_method=self.cluster_method
        if not cluster_exe:
            cluster_exe=self.cluster_exe
        if not num_clusters:
            num_clusters=self.num_clusters
        if not percent_truncation:
            percent_truncation=self.percent_truncation
        if not truncation_method:
            truncation_method=self.truncation_method
        if not truncation_pruning:
            truncation_pruning=self.truncation_pruning
        if not ensembles_directory:
            self.ensembles_directory=os.path.join(work_dir,"ensembles")
        else:
            self.ensembles_directory=ensembles_directory
        
        if not len(models):
            raise RuntimeError,"Cannot find any models for ensembling!" 
        if not all([os.path.isfile(m) for m in models]):
            raise RuntimeError,"Problem reading models given to Ensembler: {0}".format(models) 
        
        self.logger.info('Ensembling models in directory: {0}'.format(self.work_dir))
    
        # Create final ensembles directory
        if not os.path.isdir(self.ensembles_directory):
            os.mkdir(self.ensembles_directory)

        self.ensembles = []
        self.ensembles_data = []
        for cluster, cluster_data in zip(*self.cluster_models(models=models,
                                                              cluster_method=cluster_method,
                                                              num_clusters=num_clusters,
                                                              cluster_exe=cluster_exe,
                                                              nproc=nproc)):
            if len(cluster) < 2:
                self.logger.info("Cannot truncate cluster {0} as < 2 models!".format(cluster_data['cluster_num']))
                continue
            for truncated_models, truncated_models_data in zip(*self.truncate_models(cluster,
                                                                                     cluster_data,
                                                                                     truncation_method=truncation_method,
                                                                                     truncation_pruning=truncation_pruning,
                                                                                     percent_truncation=percent_truncation)):
                for subcluster, subcluster_data in zip(*self.subcluster_models(truncated_models,
                                                                               truncated_models_data,
                                                                               subcluster_program=self.subcluster_program,
                                                                               subcluster_exe=self.subcluster_program,
                                                                               ensemble_max_models=self.ensemble_max_models)):
                    for ensemble, ensemble_data in zip(*self.edit_side_chains(subcluster, subcluster_data, self.ensembles_directory)):
                        self.ensembles.append(ensemble)
                        self.ensembles_data.append(ensemble_data)
        
        return self.ensembles

    def generate_ensembles_homologs(self,
                                    models,
                                    alignment_file = None,
                                    percent_truncation = None,
                                    truncation_method = None,
                                    ensembles_directory = None,
                                    work_dir = None,
                                    nproc = None):
        
        # Work dir set each time
        if not work_dir: raise RuntimeError,"Need to set work_dir!"
        self.work_dir=work_dir
        
        if not percent_truncation:
            percent_truncation=self.percent_truncation
        if not truncation_method:
            truncation_method=self.truncation_method
        if not ensembles_directory:
            self.ensembles_directory=os.path.join(work_dir,"ensembles")
        else:
            self.ensembles_directory=ensembles_directory
        
        if not len(models):
            raise RuntimeError,"Cannot find any models for ensembling!" 
        if not all([os.path.isfile(m) for m in models]):
            raise RuntimeError,"Problem reading models given to Ensembler: {0}".format(models) 
        
        self.logger.info('Ensembling models in directory: {0}'.format(self.work_dir))
    
        # Create final ensembles directory
        if not os.path.isdir(self.ensembles_directory): os.mkdir(self.ensembles_directory)
        
        # If we've been given an alignment file, trim the models down to a core
        if alignment_file:
            core_models_dir = os.path.join(work_dir,'core_models')
            core_models = model_core_from_alignment(models, alignment_file=alignment_file, work_dir=core_models_dir)
        else: core_models = models
        
        self.ensembles = []
        self.ensembles_data = []
        for truncated_models, truncated_models_data in zip(*self.truncate_models(core_models,
                                                                                 models_data={'cluster_num': 1},
                                                                                 truncation_method=truncation_method,
                                                                                 truncation_pruning=None,
                                                                                 percent_truncation=percent_truncation,
                                                                                 homologs=True
                                                                                 )):
            tlevel = truncated_models_data['truncation_level']
            ensemble_dir = os.path.join(truncated_models_data['truncation_dir'],
                                        "ensemble_{0}".format(tlevel))
            os.mkdir(ensemble_dir)
            os.chdir(ensemble_dir)
             
            # Need to create an alignment file for theseus
            basename = "e{0}".format(tlevel)
            pre_ensemble = self.align_models(truncated_models, basename=basename, work_dir=ensemble_dir)
            if not pre_ensemble:
                self.logger.critical("Skipping ensemble {0} due to error with Theseus".format(basename))
                continue
            pre_ensemble_data = copy.copy(truncated_models_data)
            pre_ensemble_data['subcluster_radius_threshold'] = 0
             
            for ensemble, ensemble_data in zip(*self.edit_side_chains(pre_ensemble, pre_ensemble_data, self.ensembles_directory)):
                self.ensembles.append(ensemble)
                self.ensembles_data.append(ensemble_data)
        
        return self.ensembles

    def generate_thresholds(self,var_by_res,percent_interval):
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
            self.logger.debug("Got {0} thresholds: {1}".format( len(self.thresholds), self.thresholds ))
            return

        # List of variances ordered by residue index
        var_list=[var for (_,_,var) in var_by_res]

        length = len(var_list)
        if length == 0:
            msg = "Error generating thresholds, got len: {0}".format(length)
            self.logger.critical(msg)
            raise RuntimeError,msg

        # How many residues should fit in each bin
        # NB - Should round up not down with int!
        chunk_size=int( ( float(length)/100 ) *float(percent_interval) )
        if chunk_size < 1:
            msg = "Error generating thresholds, got < 1 AA in chunk_size"
            self.logger.critical(msg)
            raise RuntimeError,msg

        ## try to find intervals for truncation
        truncation_thresholds=self._generate_thresholds(var_list, chunk_size)
        
        # Jens' new untested method
        #truncation_thresholds=self._generate_thresholds2(var_list, chunk_size)
        
        self.logger.debug("Got {0} thresholds: {1}".format( len(truncation_thresholds), truncation_thresholds ))

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
    
    def prune_residues(self,residues,chunk_size=1,allowed_gap=2):
        """Remove any residues that are < chunk_size where the gap before and after is > allowed_gap"""
        
        assert chunk_size > 0 and allowed_gap > 0, "chunk_size and allowed_gap must be > 0!: {0} {1}".format(chunk_size,allowed_gap)
        
        if not len(residues): return residues,None
        lenr=len(residues)
        if lenr <= chunk_size:
            return [],residues
        
        # Build up a list of residues to remove
        to_remove=[]
        start=residues[0]
        last=residues[0]
        this_residue=None
        last_chunk_end=residues[0]-(allowed_gap+1) # make sure starting gap is bigger than allowed 
        i=1
        idxLast=lenr-1
        while i <= idxLast:
            this_residue=residues[i]
            if i==idxLast or this_residue != last+1:
                if i==idxLast: # Need to fiddle things at the end
                    if this_residue != last + 1: # handle last if a singleton
                        start=this_residue
                        last_chunk_end=last
                    last=this_residue
                    postgap=allowed_gap+1 # at end so just larger then gap
                else:
                    postgap=(this_residue-last)-1
                
                pregap=(start-last_chunk_end)-1
                this_chunk_size=(last-start)+1
                
                # remove if it satisfies the requirements
                if (this_chunk_size <= chunk_size and pregap >= allowed_gap and postgap >= allowed_gap):
                    chunk=[x for x in range(start,last+1)]
                    to_remove += chunk
                
                # reset start and last_chunk_end
                start=this_residue
                last_chunk_end=last
                
            # update loop
            last=this_residue
            i+=1
        
        # Remove the chunks and return
        if len(to_remove):
            return [r for r in residues if r not in to_remove],to_remove
        else:
            return residues,None

    def subcluster_models(self,
                           truncated_models,
                           truncated_models_data,
                           subcluster_program=None,
                           subcluster_exe=None,
                           ensemble_max_models=None):
        if self.subcluster_method=="ORIGINAL":
            f=self.subcluster_models_fixed_radii
        elif self.subcluster_method=="FIXED_ENSEMBLES":
            f=self.subcluster_models_floating_radii
        else:
            assert False
        return f(truncated_models,
                 truncated_models_data,
                 subcluster_program,
                 subcluster_exe,
                 ensemble_max_models)
        
    def subcluster_models_fixed_radii(self,
                                      truncated_models,
                                      truncated_models_data,
                                      subcluster_program=None,
                                      subcluster_exe=None,
                                      ensemble_max_models=None):
        
        # Theseus only works with > 3 residues
        if truncated_models_data['truncation_num_residues'] <= 2: return [],[]
        
        radius_thresholds=self.subcluster_radius_thresholds
        ensembles=[]
        ensembles_data=[]
        
        # Use first model to get data on level
        cluster_num=truncated_models_data['cluster_num']
        truncation_level=truncated_models_data['truncation_level']
        truncation_dir=truncated_models_data['truncation_dir']
            
        # Run maxcluster to generate the distance matrix
        if subcluster_program=='maxcluster':
            clusterer = subcluster.MaxClusterer(self.subcluster_exe)
        else:
            assert False
        clusterer.generate_distance_matrix(truncated_models)
        #clusterer.dump_matrix(os.path.join(truncation_dir,"subcluster_distance.matrix")) # for debugging

        # Loop through the radius thresholds
        num_previous_models=-1 # set to -1 so comparison always false on first pass
        for radius in radius_thresholds:

            self.logger.debug("subclustering files under radius: {0}".format(radius))

            # Get list of pdbs clustered according to radius threshold
            cluster_files = clusterer.cluster_by_radius(radius)
            self.logger.debug("Clustered {0} files".format(len(cluster_files)))
            if len(cluster_files) < 2:
                self.logger.debug('Clustered fewer than 2 files using radius {0} -  in truncation dir {1} SKIPPING'.format(radius,truncation_dir))
                continue

            # For naming all files
            basename='c{0}_tl{1}_r{2}'.format(cluster_num, truncation_level, radius)

            # Check if there are the same number of models in this ensemble as the previous one - if so
            # the ensembles will be identical and we can skip this one
            if num_previous_models == len(cluster_files):
                self.logger.debug( 'Number of decoys in cluster ({0}) is the same as under previous threshold so excluding cluster {1}'.format( len( cluster_files ), basename ) )
                continue
            else:
                num_previous_models = len(cluster_files)

            # Got files so create the directories
            subcluster_dir = os.path.join(truncation_dir, 'subcluster_{0}'.format(radius))
            os.mkdir(subcluster_dir)
            os.chdir(subcluster_dir)

            # Restrict cluster to max_ensemble_models
            if len(cluster_files) > ensemble_max_models:
                self.logger.debug("{0} files in cluster so truncating list to first {1}".format(len(cluster_files), ensemble_max_models))
                cluster_files = cluster_files[:ensemble_max_models]
                
            cluster_file = self.align_models(cluster_files)
            if not cluster_file:
                msg="Error running theseus on ensemble {0} in directory: {1}\nSkipping subcluster: {0}".format(basename,
                                                                                                               subcluster_dir)
                self.logger.critical(msg)
                continue
             
            ensemble = os.path.join(subcluster_dir, basename+'.pdb')
            shutil.move(cluster_file, ensemble)

            # The data we've collected is the same for all pdbs in this level so just keep using the first  
            ensemble_data=copy.copy(truncated_models_data)
            ensemble_data['subcluster_num_models'] = len( cluster_files )
            ensemble_data['subcluster_radius_threshold'] = radius
            ensemble_data['ensemble_pdb'] = ensemble

            # Get the centroid model name from the list of files given to theseus - we can't parse
            # the pdb file as theseus truncates the filename
            ensemble_data['subcluster_centroid_model']=os.path.abspath(cluster_files[0])
            
            ensembles.append(ensemble)
            ensembles_data.append(ensemble_data)
        
        return ensembles,ensembles_data

    def subcluster_models_floating_radii(self,
                                         truncated_models,
                                         truncated_models_data,
                                         subcluster_program=None,
                                         subcluster_exe=None,
                                         ensemble_max_models=None):
        self.logger.info("subclustering with floaing radii")

        # Run maxcluster to generate the distance matrix
        if subcluster_program=='maxcluster':
            clusterer = subcluster.MaxClusterer(self.subcluster_exe)
        else:
            assert False
        clusterer.generate_distance_matrix(truncated_models)
        #clusterer.dump_matrix(os.path.join(truncation_dir,"subcluster_distance.matrix")) # for debugging
        
        subclusters=[]
        subclusters_data=[]
        clusters=[]
        radii=[]
        len_truncated_models=len(truncated_models)
        for i in range(len(self.subcluster_radius_thresholds)):
            radius=None
            nmodels=None
            if i > 0 and radii[i-1] > self.subcluster_radius_thresholds[i]:
                radius = radii[i-1]
                nmodels = len(clusters[i-1])
                cluster_files, radius = self._subcluster_nmodels(nmodels, radius, clusterer, direction='up',increment=1)
            else:
                radius = self.subcluster_radius_thresholds[i]
                cluster_files = clusterer.cluster_by_radius(radius)
                
            cluster_files=tuple(sorted(cluster_files)) # Need to sort so that we can check if we've had this lot before
            cluster_size=len(cluster_files)

            if radius in radii or cluster_size==1:
                # Increase radius till we have one more than the last one
                if cluster_size==1:
                    nmodels=2
                else:
                    radius=radii[i-1]
                    nmodels=len(clusters[i-1])+1
                cluster_files, radius = self._subcluster_nmodels(nmodels, radius, clusterer, direction='up',increment=1)
                cluster_files=sorted(cluster_files)
            elif cluster_size >= ensemble_max_models or cluster_files in clusters:
                # Randomly pick ensemble_max_models
                cluster_files = self._pick_nmodels(cluster_files, clusters, ensemble_max_models)
                if not cluster_files:
                    self.logger.debug('Could not cluster files under radius: {0} - could not find different models'.format(radius))
                    
            
            # Need to check in case we couldn't cluster under this radius
            cluster_size=len(cluster_files)
            if cluster_size==1 or radius in radii:
                self.logger.debug('Could not cluster files under radius: {0} - got {1} files'.format(radius,len(cluster_files)))
                break
            
            self.logger.debug('Subclustering {0} files under radius {1}'.format(cluster_size,radius))
            try:
                cluster_ensemble, data = self._subcluster_radius(list(cluster_files), radius, truncated_models_data)
            except RuntimeError:
                self.logger.debug('Could not cluster files under radius: {0} as theseus failed'.format(radius,len(cluster_files)))
                # If theseus fails, we just move
                break
            
            subclusters.append(cluster_ensemble)
            subclusters_data.append(data)
            clusters.append(tuple(cluster_files)) # append as tuple so it is hashable
            radii.append(radius)
            if cluster_size==len_truncated_models: break
            
        return subclusters, subclusters_data
    
    def _pick_nmodels(self, models, clusters, ensemble_max_models):
        MAXTRIES=50
        tries = 0
        clusters=set(clusters)
        nmodels=min(len(models),ensemble_max_models)
        while True:
            subcluster = random.sample(models,nmodels)
            subcluster = tuple(sorted(subcluster))
            if subcluster not in clusters: break
            tries += 1
            if tries >= MAXTRIES: return None
        return subcluster
    
    def _subcluster_nmodels(self,nmodels,radius,clusterer,direction,increment):
        """
        """
        MINRADIUS=0.0001
        MAXRADIUS=100
        subcluster_models=clusterer.cluster_by_radius(radius)
        len_models=len(subcluster_models)
        self.logger.debug("_subcluster_nmodels: {0} {1} {2} {3} {4}".format(len_models,nmodels,radius,direction,increment))
        if len_models == nmodels or radius < MINRADIUS or radius > MAXRADIUS:
            self.logger.debug("_subcluster_nmodels returning: nmodels: {0} radius: {1}".format(len_models,radius ))
            return subcluster_models, radius
        
        def lower_increment(increment):
            increment = increment / float(10)
            if increment <= 0.00001: raise RuntimeError,"increment out of bounds"
            return increment
        
        # Am sure the logic could be improved here, but it seems to  work
        try:
            if len_models > nmodels:
                # If we have more models than we want and we are increasing the radius, we've overshot, so we need to
                # decrease the radius but by a smaller increment
                # If the radius is the same as the increment, we need to decrease the incrememnt before we subtract it
                # as both of the above require decreasing the increment we have one test and just change the direction
                # for the overshoot
                if direction == 'up' or abs(radius-increment) < 0.0000001:
                    if direction == 'up': direction = 'down'
                    increment = lower_increment(increment)
                radius -= increment
            elif len_models < nmodels:
                if direction == 'down' :
                    direction = 'up'
                    increment = lower_increment(increment)
                radius += increment
        except RuntimeError:
            # Can't get a match so just return what we have
            self.logger.debug("_subcluster_nmodels exceeded increment. Returning: nmodels: {0} radius: {1}".format(len(subcluster_models),radius ))
            return subcluster_models, radius
            
        return self._subcluster_nmodels(nmodels, radius, clusterer, direction, increment)
    
    def _subcluster_radius(self,models,radius,truncated_models_data):
        # Extract data from dictionary
        cluster_num=truncated_models_data['cluster_num']
        truncation_level=truncated_models_data['truncation_level']
        truncation_dir=truncated_models_data['truncation_dir']

        # Got files so create the directories
        subcluster_dir = os.path.join(truncation_dir, 'subcluster_{0}'.format(radius))
        os.mkdir(subcluster_dir)
        os.chdir(subcluster_dir)

        basename='c{0}_tl{1}_r{2}'.format(cluster_num, truncation_level, radius)
        cluster_file = self.align_models(models)
        if not cluster_file:
            msg="Error running theseus on ensemble {0} in directory: {1}\nSkipping subcluster: {0}".format(basename,
                                                                                                subcluster_dir)
            self.logger.critical(msg)
            raise RuntimeError,msg
        
        ensemble = os.path.join(subcluster_dir, basename+'.pdb')
        shutil.move(cluster_file, ensemble)

        # The data we've collected is the same for all pdbs in this level so just keep using the first  
        subcluster_data=copy.copy(truncated_models_data)
        subcluster_data['subcluster_num_models'] = len(models)
        subcluster_data['subcluster_radius_threshold'] = radius
        subcluster_data['ensemble_pdb'] = ensemble

        # Get the centroid model name from the list of files given to theseus - we can't parse
        # the pdb file as theseus truncates the filename
        subcluster_data['subcluster_centroid_model'] = os.path.abspath(models[0])
        return ensemble, subcluster_data
  
    def truncate_models(self,
                        models,
                        models_data=None,
                        truncation_method=None,
                        percent_truncation=None,
                        truncation_pruning='none',
                        homologs=False
                        ):
        
        assert len(models) > 1,"Cannot truncate as < 2 models!"
        assert truncation_method and percent_truncation,"Missing arguments: {0} {1}".format(truncation_method,percent_truncation)

        # Create the directories we'll be working in
        truncate_dir =  os.path.join(self.work_dir, 'truncate_{0}'.format(models_data['cluster_num']))
        os.mkdir(truncate_dir)
        os.chdir(truncate_dir)
        
        # Calculate variances between pdb - and align them if necessary
        run_theseus = theseus.Theseus(work_dir=truncate_dir, theseus_exe=self.theseus_exe)
        run_theseus.align_models(models, homologs=homologs)
        #if homologs: models = run_theseus.aligned_models
        var_by_res = run_theseus.var_by_res()
        
        # Calculate which residues to keep under the different methods
        truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs=None,None,None,None
        if truncation_method=='percent':
            truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = self._calculate_residues_percent(var_by_res,percent_truncation)
        elif truncation_method=='thresh':
            truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs = self._calculate_residues_thresh(var_by_res,percent_truncation)
        else:
            raise RuntimeError,"Unrecognised ensembling mode: {0}".format(truncation_method)

        truncated_models=[]
        truncated_models_data=[]
        pruned_residues=None
        for tlevel,tvar,tresidues,tresidue_idxs in zip(truncation_levels, truncation_variances, truncation_residues, truncation_residue_idxs):
            # Prune singletone/doubletone etc. residues if required
# This code is currently disabled as it relies on resseqs and we now use residue indexes to allow us to work with pdbs with different reseqs
#             self.logger.debug("truncation_pruning: {0}".format(truncation_pruning))
#             if truncation_pruning=='single':
#                 tresidues,pruned_residues=self.prune_residues(tresidues, chunk_size=1, allowed_gap=2)
#                 if pruned_residues: self.logger.debug("prune_residues removing: {0}".format(pruned_residues))
#             elif truncation_pruning=='none':
#                 pass
#             else:
#                 raise RuntimeError,"Unrecognised truncation_pruning: {0}".format(truncation_pruning)
            
            # Skip if there are no residues
            if not tresidue_idxs:
                self.logger.debug("Skipping truncation level {0} with variance {1} as no residues".format(tlevel,tvar))
                continue
            
            trunc_dir = os.path.join(truncate_dir, 'tlevel_{0}'.format(tlevel))
            os.mkdir(trunc_dir)
            self.logger.info( 'truncating at: {0} in directory {1}'.format(tvar,trunc_dir))

            # list of models for this truncation level
            level_models = []
            for infile in models:
                pdbname = os.path.basename(infile)
                pdbout = os.path.join(trunc_dir, pdbname)
                # Loop through PDB files and create new ones that only contain the residues left after truncation
                pdb_edit.select_residues(pdbin=infile, pdbout=pdbout, tokeep_idx=tresidue_idxs)
                level_models.append(pdbout)
            
            # Add the model
            truncated_models.append(level_models)

            # Add the data
            model_data=copy.copy(models_data)
            model_data['truncation_level']        = tlevel
            model_data['truncation_variance']     = tvar
            model_data['truncation_residues']     = tresidues
            model_data['truncation_num_residues'] = len(tresidues)
            model_data['truncation_dir']          = trunc_dir
            model_data['percent_truncation']      = percent_truncation
            model_data['truncation_method']       = truncation_method
            model_data['truncation_pruning']      = truncation_pruning
            model_data['pruned_residues']         = pruned_residues
            
            truncated_models_data.append(model_data)
            
        return truncated_models, truncated_models_data

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
        
        cls.theseus_exe=ample_util.find_exe('theseus')
        cls.spicker_exe=ample_util.find_exe('spicker')
        cls.maxcluster_exe=ample_util.find_exe('maxcluster')


        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        formatter = logging.Formatter('%(message)s')
        ch.setFormatter(formatter)
        root.addHandler(ch)

        return

    def testPruneResidues(self):
        """Test we can reproduce the original thresholds"""
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()

        ensembler=Ensembler()
        
        residues=[]
        pres,pruned=ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres,[])
        self.assertEqual(pruned,None)
          
        residues=[1]
        pres,pruned=ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres,[])
        self.assertEqual(pruned,[1])
         
        residues=[1,2]
        pres,pruned=ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres,[])
        self.assertEqual(pruned,[1,2])
         
        residues=[1,2,4]
        pres,pruned=ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres,[1,2,4])
        self.assertEqual(pruned,None)
         
        # Big enough gap 
        residues=[1,2,3,4,8,13,14,15]
        pres,pruned=ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres,[1,2,3,4,13,14,15])
        self.assertEqual(pruned,[8])
           
        residues=[1,2,3,4,8,9,13,14,15]
        pres,pruned=ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres,[1,2,3,4,13,14,15])
        self.assertEqual(pruned,[8,9])
           
        residues=[1,2,3,4,8,9,10,13,14,15]
        pres,pruned=ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres,residues)
        self.assertEqual(pruned,None)
           
        # end gap not big enough
        residues=[1,2,3,4,8,9,11,12,13]
        pres,pruned=ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres,residues)
        self.assertEqual(pruned,None)
           
        # Lone residue at start
        residues=[1,11,12,13]
        pres,pruned=ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres,[11,12,13])
        self.assertEqual(pruned,[1])
        
        # Lone residue at end
        residues=[11,12,13,19]
        pres,pruned=ensembler.prune_residues(residues, chunk_size=2, allowed_gap=2)
        self.assertEqual(pres,[11,12,13])
        self.assertEqual(pruned,[19])
          
        # Mixed
        residues=[1,3,4,7,10,11,13,15,16,19]
        pres,pruned=ensembler.prune_residues(residues, chunk_size=1, allowed_gap=2)
        self.assertEqual(pres,[1,3,4,10,11,13,15,16])
        self.assertEqual(pruned,[7,19])

        return

    def testThresholds(self):
        """Test we can reproduce the original thresholds"""
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()

        ensembler=Ensembler()
        
        ensembler.work_dir=os.path.join(self.tests_dir,"genthresh1")
        if os.path.isdir(ensembler.work_dir):
            shutil.rmtree(ensembler.work_dir)
        os.mkdir(ensembler.work_dir)
        
        ensembler.theseus_exe=self.theseus_exe
        percent_interval=5
        mdir=os.path.join(self.testfiles_dir,"models")
        cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        var_by_res=ensembler.calculate_variances(cluster_models)
        thresholds=ensembler.generate_thresholds(var_by_res,percent_interval)
        
        self.assertEqual(30,len(thresholds),thresholds)
        reft=[0.02235, 0.041896, 0.08155, 0.085867, 0.103039, 0.185265, 0.7058, 1.772884, 3.152793,
              4.900255, 6.206563, 10.250043, 10.772177, 16.071345, 18.292366, 22.432564, 23.265938,
              25.348501, 27.37951, 35.877391, 37.227859, 41.759805, 49.037625, 50.992344, 56.091738,
              58.09387, 59.917999, 70.052007, 80.235358, 86.028994]
        
        self.assertTrue(all([ abs(r-c) < 0.0001 for r,c in zip(reft,thresholds)]),
                         "Wrong truncation thresholds: {0}".format(thresholds))
        
        shutil.rmtree(ensembler.work_dir)
        return
    
    def testResiduesThresh(self):
        """Test we can calculate the original list of residues"""
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler=Ensembler()
        
        # This test for percent
        ensembler.work_dir=os.path.join(self.tests_dir,"genthresh2")
        if os.path.isdir(ensembler.work_dir):
            shutil.rmtree(ensembler.work_dir)
        os.mkdir(ensembler.work_dir)
        os.chdir(ensembler.work_dir)
        
        ensembler.theseus_exe=self.theseus_exe
        percent_interval=5
        mdir=os.path.join(self.testfiles_dir,"models")
        cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        
        var_by_res=ensembler.calculate_variances(cluster_models)
        truncation_levels, truncation_variances, truncation_residues,  truncation_residues_idxs = ensembler._calculate_residues_thresh(var_by_res,percent_interval)
        
        self.assertEqual(truncation_levels,
                         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30])
        
        
        refv = [86.028994, 80.235358, 70.052007, 59.917999, 58.09387, 56.091738, 50.992344, 49.037625, 41.759805, 37.227859, 35.877391,
                27.37951, 25.348501, 23.265938, 22.432564, 18.292366, 16.071345, 10.772177, 10.250043, 6.206563, 4.900255, 3.152793,
                1.772884, 0.7058, 0.185265, 0.103039, 0.085867, 0.08155, 0.041896, 0.02235]
        self.assertTrue(all([ abs(r-c) < 0.0001 for r,c in zip(refv,truncation_variances)]),
                         "Wrong truncation variances: {0}".format(truncation_variances))

        residues=[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59],
                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58],
                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 56, 57],
                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 54, 56, 57],
                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 57],
                  [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53],
                  [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 50, 53],
                  [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 53],
                  [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47],
                  [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46],
                  [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45],
                  [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
                  [7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
                  [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
                  [9, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
                  [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
                  [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
                  [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41],
                  [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41],
                  [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
                  [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
                  [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
                  [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
                  [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
                  [23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
                  [25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
                  [26, 27, 28, 30, 31, 32, 33, 34],
                  [26, 27, 30, 31, 32, 33],
                  [26, 27, 31, 32],
                  [31, 32]]
        self.assertEqual(residues,truncation_residues)

        shutil.rmtree(ensembler.work_dir)
        return
    
    def testResiduesPercent(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler=Ensembler()

        ensembler.work_dir=os.path.join(self.tests_dir,"genthresh3")
        if os.path.isdir(ensembler.work_dir):
            shutil.rmtree(ensembler.work_dir)
        os.mkdir(ensembler.work_dir)
        os.chdir(ensembler.work_dir)
        ensembler.theseus_exe=self.theseus_exe
        ensembler.percent_interval=5
        mdir=os.path.join(self.testfiles_dir,"models")
        cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        var_by_res=ensembler.calculate_variances(cluster_models)
        truncation_levels, truncation_variances, truncation_residues, truncation_residues_idxs = ensembler._calculate_residues_percent(var_by_res,percent_interval=5)


        residues=[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59],
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
            [26, 27, 30, 31, 32]]
        
        for i in range(len(residues)):
            self.assertEqual(residues[i],truncation_residues[i],"Mismatching residues for level {0}\n{1}\n{2}".format(i,residues[i],truncation_residues[i]))

        self.assertEqual(truncation_levels,
                         [100, 95, 90, 85, 80, 75, 69, 64, 59, 54, 49, 44, 39, 34, 29, 24, 19, 14, 8])
        
        refv=[74.272253, 59.917999, 56.260987, 50.992344, 43.268754, 37.227859, 30.767268, 25.348501,
              23.050357, 18.292366, 15.287348, 10.250043, 5.66545, 3.152793, 0.874899, 0.185265, 0.090917,
              0.08155, 0.039639]
        self.assertTrue(all([ abs(r-c) < 0.0001 for r,c in zip(refv,truncation_variances)]),
                         "Wrong truncation variances: {0}".format(truncation_variances))
        shutil.rmtree(ensembler.work_dir)
        return

    def testClustering(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler=Ensembler()

        ensembler.work_dir=os.path.join(self.tests_dir,"genthresh4")
        if os.path.isdir(ensembler.work_dir):
            shutil.rmtree(ensembler.work_dir)
        os.mkdir(ensembler.work_dir)
        
        ensembler.theseus_exe=self.theseus_exe
        
        mdir=os.path.join(self.testfiles_dir,"models")
        models=glob.glob(mdir+os.sep+"*.pdb")
        
        cluster_models,cluster_data=ensembler.cluster_models(models=models,
                                                             cluster_method='spicker',
                                                             num_clusters=3,
                                                             cluster_exe=self.spicker_exe)
        
        names=sorted([os.path.basename(m) for m in cluster_models[0]])
        self.assertEqual(names,
                         sorted(['5_S_00000005.pdb', '4_S_00000005.pdb', '5_S_00000004.pdb', '4_S_00000002.pdb',
                          '4_S_00000003.pdb', '3_S_00000006.pdb', '3_S_00000004.pdb', '2_S_00000005.pdb',
                          '2_S_00000001.pdb', '3_S_00000003.pdb', '1_S_00000005.pdb', '1_S_00000002.pdb', 
                          '1_S_00000004.pdb']) )
        
        d = cluster_data[2]
        self.assertEqual(os.path.basename(d['cluster_centroid']),'2_S_00000003.pdb')
        self.assertEqual(d['cluster_num_models'],1)
        shutil.rmtree(ensembler.work_dir)
        return
    
    def testEnsemblingPercent(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler=Ensembler()

        work_dir=os.path.join(self.tests_dir,"genthresh5")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe=self.theseus_exe
        ensembler.cluster_exe=self.spicker_exe
        ensembler.subcluster_exe=self.maxcluster_exe
        
        mdir=os.path.join(self.testfiles_dir,"models")
        models=glob.glob(mdir+os.sep+"*.pdb")

        num_clusters=1
        cluster_method='spicker'
        percent_truncation=5
        truncation_method="percent"
        ensembles=ensembler.generate_ensembles(models,
                                                 cluster_method=cluster_method,
                                                 cluster_exe=self.spicker_exe,
                                                 num_clusters=num_clusters,
                                                 percent_truncation=percent_truncation,
                                                 truncation_method=truncation_method,
                                                 work_dir=work_dir)

        eref=sorted(['c1_tl100_r2_allatom.pdb', 'c1_tl100_r2_reliable.pdb', 'c1_tl100_r2_polyAla.pdb', 'c1_tl100_r3_allatom.pdb',
              'c1_tl100_r3_reliable.pdb', 'c1_tl100_r3_polyAla.pdb', 'c1_tl95_r2_allatom.pdb', 'c1_tl95_r2_reliable.pdb',
              'c1_tl95_r2_polyAla.pdb', 'c1_tl95_r3_allatom.pdb', 'c1_tl95_r3_reliable.pdb', 'c1_tl95_r3_polyAla.pdb', 'c1_tl90_r1_allatom.pdb',
              'c1_tl90_r1_reliable.pdb', 'c1_tl90_r1_polyAla.pdb', 'c1_tl90_r2_allatom.pdb', 'c1_tl90_r2_reliable.pdb', 'c1_tl90_r2_polyAla.pdb',
              'c1_tl90_r3_allatom.pdb', 'c1_tl90_r3_reliable.pdb', 'c1_tl90_r3_polyAla.pdb', 'c1_tl85_r1_allatom.pdb', 'c1_tl85_r1_reliable.pdb',
              'c1_tl85_r1_polyAla.pdb', 'c1_tl85_r2_allatom.pdb', 'c1_tl85_r2_reliable.pdb', 'c1_tl85_r2_polyAla.pdb', 'c1_tl85_r3_allatom.pdb',
              'c1_tl85_r3_reliable.pdb', 'c1_tl85_r3_polyAla.pdb', 'c1_tl80_r1_allatom.pdb', 'c1_tl80_r1_reliable.pdb', 'c1_tl80_r1_polyAla.pdb',
              'c1_tl80_r2_allatom.pdb', 'c1_tl80_r2_reliable.pdb', 'c1_tl80_r2_polyAla.pdb', 'c1_tl80_r3_allatom.pdb', 'c1_tl80_r3_reliable.pdb',
              'c1_tl80_r3_polyAla.pdb', 'c1_tl75_r1_allatom.pdb', 'c1_tl75_r1_reliable.pdb', 'c1_tl75_r1_polyAla.pdb', 'c1_tl75_r2_allatom.pdb',
              'c1_tl75_r2_reliable.pdb', 'c1_tl75_r2_polyAla.pdb', 'c1_tl75_r3_allatom.pdb', 'c1_tl75_r3_reliable.pdb', 'c1_tl75_r3_polyAla.pdb',
              'c1_tl69_r1_allatom.pdb', 'c1_tl69_r1_reliable.pdb', 'c1_tl69_r1_polyAla.pdb', 'c1_tl69_r2_allatom.pdb', 'c1_tl69_r2_reliable.pdb',
              'c1_tl69_r2_polyAla.pdb', 'c1_tl64_r1_allatom.pdb', 'c1_tl64_r1_reliable.pdb', 'c1_tl64_r1_polyAla.pdb', 'c1_tl64_r2_allatom.pdb',
              'c1_tl64_r2_reliable.pdb', 'c1_tl64_r2_polyAla.pdb', 'c1_tl59_r1_allatom.pdb', 'c1_tl59_r1_reliable.pdb', 'c1_tl59_r1_polyAla.pdb',
              'c1_tl59_r2_allatom.pdb', 'c1_tl59_r2_reliable.pdb', 'c1_tl59_r2_polyAla.pdb', 'c1_tl54_r1_allatom.pdb', 'c1_tl54_r1_reliable.pdb',
              'c1_tl54_r1_polyAla.pdb', 'c1_tl54_r2_allatom.pdb', 'c1_tl54_r2_reliable.pdb', 'c1_tl54_r2_polyAla.pdb', 'c1_tl49_r1_allatom.pdb',
              'c1_tl49_r1_reliable.pdb', 'c1_tl49_r1_polyAla.pdb', 'c1_tl49_r2_allatom.pdb', 'c1_tl49_r2_reliable.pdb', 'c1_tl49_r2_polyAla.pdb',
              'c1_tl44_r1_allatom.pdb', 'c1_tl44_r1_reliable.pdb', 'c1_tl44_r1_polyAla.pdb', 'c1_tl44_r2_allatom.pdb', 'c1_tl44_r2_reliable.pdb',
              'c1_tl44_r2_polyAla.pdb', 'c1_tl39_r1_allatom.pdb', 'c1_tl39_r1_reliable.pdb', 'c1_tl39_r1_polyAla.pdb', 'c1_tl34_r1_allatom.pdb',
              'c1_tl34_r1_reliable.pdb', 'c1_tl34_r1_polyAla.pdb', 'c1_tl29_r1_allatom.pdb', 'c1_tl29_r1_reliable.pdb', 'c1_tl29_r1_polyAla.pdb', 
              'c1_tl24_r1_allatom.pdb', 'c1_tl24_r1_reliable.pdb', 'c1_tl24_r1_polyAla.pdb', 'c1_tl19_r1_allatom.pdb', 'c1_tl19_r1_reliable.pdb',
              'c1_tl19_r1_polyAla.pdb', 'c1_tl14_r1_allatom.pdb', 'c1_tl14_r1_reliable.pdb', 'c1_tl14_r1_polyAla.pdb', 'c1_tl8_r1_allatom.pdb',
              'c1_tl8_r1_reliable.pdb', 'c1_tl8_r1_polyAla.pdb'])
        self.assertEqual(sorted([os.path.basename(m) for m in ensembles]),eref)
        d = ensembler.ensembles_data[5]

        self.assertEqual(d['percent_truncation'],percent_truncation)
        self.assertEqual(d['truncation_method'],truncation_method)
        self.assertEqual(d['cluster_method'],cluster_method)
        self.assertEqual(d['num_clusters'],num_clusters)
        self.assertTrue(abs(d['truncation_variance']-13.035172) < 0001)
        self.assertEqual(d['ensemble_num_atoms'],984)
        self.assertEqual(os.path.basename(d['subcluster_centroid_model']),'5_S_00000005.pdb')
        
        shutil.rmtree(ensembler.work_dir)
        return
    
    def testEnsemblingThresh(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler=Ensembler()

        work_dir=os.path.join(self.tests_dir,"genthresh6")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe=self.theseus_exe
        ensembler.cluster_exe=self.spicker_exe
        ensembler.subcluster_exe=self.maxcluster_exe
        
        mdir=os.path.join(self.testfiles_dir,"models")
        models=glob.glob(mdir+os.sep+"*.pdb")
        
        num_clusters=1
        cluster_method='spicker'
        percent_truncation=5
        truncation_method="thresh"
        ensembles=ensembler.generate_ensembles(models,
                                               cluster_method=cluster_method,
                                               cluster_exe=self.spicker_exe,
                                               num_clusters=num_clusters,
                                               percent_truncation=percent_truncation,
                                               truncation_method=truncation_method,
                                               work_dir=work_dir)
        
        self.assertEqual(len(ensembles),162,len(ensembles))
        d = ensembler.ensembles_data[5]
        
        self.assertTrue(abs(d['truncation_variance']-27.389253) < 0001)
        self.assertEqual(d['percent_truncation'],percent_truncation)
        self.assertEqual(d['truncation_method'],truncation_method)
        self.assertEqual(d['cluster_method'],cluster_method)
        self.assertEqual(d['num_clusters'],num_clusters)
        self.assertEqual(d['subcluster_radius_threshold'],3)
        self.assertEqual(d['side_chain_treatment'],ALLATOM)
        self.assertEqual(d['ensemble_num_atoms'],984)
        self.assertEqual(os.path.basename(d['subcluster_centroid_model']),'5_S_00000005.pdb')
        
        shutil.rmtree(ensembler.work_dir)
        return
    
    def testSubclusterNew1(self):
        """Divergent models"""
        
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler=Ensembler()

        work_dir=os.path.join(self.tests_dir,"genthresh7")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe=self.theseus_exe
        ensembler.cluster_exe=self.spicker_exe
        ensembler.subcluster_exe=self.maxcluster_exe
        ensembler.subcluster_method="FIXED_ENSEMBLES"
        
        mdir=os.path.join(self.testfiles_dir,"2qsk_models")
        truncated_models=glob.glob(mdir+os.sep+"*.pdb")

        truncated_models_data = { 'cluster_num'      : 1,
                                  'truncation_level' : 1,
                                  'truncation_dir'   : work_dir } 
        
        subcluster, data = ensembler.subcluster_models(truncated_models,
                                                       truncated_models_data,
                                                       subcluster_program='maxcluster',
                                                       subcluster_exe=self.maxcluster_exe,
                                                       ensemble_max_models=30)
        
        self.assertEqual(data[0]['subcluster_num_models'],2)
        self.assertTrue(abs(data[0]['subcluster_radius_threshold']-5.82) < 0.0001)
        self.assertEqual(data[1]['subcluster_num_models'],3)
        self.assertTrue(abs(data[1]['subcluster_radius_threshold']-6.82) < 0.0001)
        self.assertEqual(data[2]['subcluster_num_models'],4)
        self.assertTrue(abs(data[2]['subcluster_radius_threshold']-6.92) < 0.0001)
        shutil.rmtree(work_dir)
        return
    
    def testSubclusterNew2(self):
        """Similar models"""
        
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler=Ensembler()

        work_dir=os.path.join(self.tests_dir,"genthresh8")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe=self.theseus_exe
        ensembler.cluster_exe=self.spicker_exe
        ensembler.subcluster_exe=self.maxcluster_exe
        ensembler.subcluster_method="FIXED_ENSEMBLES"
        
        mdir=os.path.join(self.testfiles_dir,"1p9g_models")
        truncated_models=glob.glob(mdir+os.sep+"*.pdb")

        truncated_models_data = { 'cluster_num'      : 1,
                                  'truncation_level' : 1,
                                  'truncation_dir'   : work_dir } 
        
        subcluster, data = ensembler.subcluster_models(truncated_models,
                                                       truncated_models_data,
                                                       subcluster_program='maxcluster',
                                                       subcluster_exe=self.maxcluster_exe,
                                                       ensemble_max_models=30)

        self.assertEqual(data[0]['subcluster_num_models'],30,)
        self.assertTrue(abs(data[0]['subcluster_radius_threshold']-1) < 0.0001)
        self.assertEqual(data[1]['subcluster_num_models'],30)
        self.assertTrue(abs(data[1]['subcluster_radius_threshold']-2) < 0.0001)
        self.assertEqual(data[2]['subcluster_num_models'],30)
        self.assertTrue(abs(data[2]['subcluster_radius_threshold']-3) < 0.0001)
        shutil.rmtree(work_dir)
        
        return
    
    def testSubclusterNew3(self):
        """standard models"""
        
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler=Ensembler()

        work_dir=os.path.join(self.tests_dir,"genthresh9")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        ensembler.theseus_exe=self.theseus_exe
        ensembler.cluster_exe=self.spicker_exe
        ensembler.subcluster_exe=self.maxcluster_exe
        ensembler.subcluster_method="FIXED_ENSEMBLES"
        
        mdir=os.path.join(self.testfiles_dir,"models")
        truncated_models=glob.glob(mdir+os.sep+"*.pdb")

        truncated_models_data = { 'cluster_num'      : 1,
                                  'truncation_level' : 1,
                                  'truncation_dir'   : work_dir } 
        
        subcluster, data = ensembler.subcluster_models(truncated_models,
                                                       truncated_models_data,
                                                       subcluster_program='maxcluster',
                                                       subcluster_exe=self.maxcluster_exe,
                                                       ensemble_max_models=30)
        self.assertEqual(data[0]['subcluster_num_models'],2)
        self.assertTrue(abs(data[0]['subcluster_radius_threshold']-1.9) < 0.0001,"GOT {0}".format(data[0]['subcluster_radius_threshold']))
        self.assertEqual(data[1]['subcluster_num_models'],3)
        self.assertTrue(abs(data[1]['subcluster_radius_threshold']-2) < 0.0001,"GOT {0}".format(data[1]['subcluster_radius_threshold']))
        self.assertEqual(data[2]['subcluster_num_models'],8)
        self.assertTrue(abs(data[2]['subcluster_radius_threshold']-3) < 0.0001,"GOT {0}".format(data[2]['subcluster_radius_threshold']))
        shutil.rmtree(work_dir)
        return
    
    
    def XtestSubclusterNew(self):
        
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler=Ensembler()

        work_dir=os.path.join(self.tests_dir,"genthresh9")
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
        os.mkdir(work_dir)
        
        mdir="1P9G"
        truncated_models=glob.glob(mdir+os.sep+"*.pdb")

        subcluster_exe=self.maxcluster_exe
        clusterer = subcluster.MaxClusterer(subcluster_exe)
        clusterer.generate_distance_matrix(truncated_models)

        max_models=30
        radius=1
        direction="down"
        increment=0.1
        models, new_radius = ensembler._subcluster_nmodels(max_models,
                                                           radius,
                                                           clusterer,
                                                           direction,
                                                           increment)
        self.assertEqual(len(models),36)
        self.assertTrue(abs(new_radius-0.005) < 0.0001)
        return
    
    def testHomologs(self):
            
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        ensembler=Ensembler()
        ensembler.theseus_exe=self.theseus_exe
        
        work_dir=os.path.join(self.tests_dir,"homologs_test")
        if os.path.isdir(work_dir): shutil.rmtree(work_dir)
        os.mkdir(work_dir)

        pdb_list = [ '1ujb.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir,'examples','homologs',pdb) for pdb in pdb_list ]
        alignment_file = os.path.join(self.ample_dir,'examples','homologs','testthree.afasta')
        
        ensembles = ensembler.generate_ensembles_homologs(models, alignment_file=alignment_file, work_dir=work_dir)
         
        self.assertEqual(len(ensembles),57)
        
        shutil.rmtree(work_dir)
        return
    
    def testCoreFromAlignment(self):
        
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()
        work_dir=os.path.join(self.tests_dir,"homologs_core")
        if os.path.isdir(work_dir): shutil.rmtree(work_dir)
        os.mkdir(work_dir)

        models = glob.glob(os.path.join("/opt/ample-dev1/examples/homologs","*.pdb"))
        alignment_file = os.path.join("/opt/ample-dev1/examples/homologs","testthree.afasta")
        
        core_models = model_core_from_alignment(models, alignment_file, work_dir=work_dir)
        
        got = {}
        for m in core_models:
            name = os.path.splitext(os.path.basename(m))[0]
            resseqd = pdb_edit.resseq(m)
            resseq = resseqd[resseqd.keys()[0]]
            got[name] = resseq

        ref = { '1ujb_core':
               [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 
                32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 
                56, 57, 58, 59, 60, 61, 62, 63, 64, 67, 69, 70, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 84, 85, 
                86, 87, 88, 89, 90, 92, 93, 94, 95, 96, 97, 98, 99, 101, 102, 103, 104, 105, 106, 107, 108, 109, 
                110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 124, 126, 127, 128, 129, 130, 
                131, 132, 133, 134, 135, 136, 137, 138, 139, 141, 142, 143, 144, 145, 146, 147, 149, 150, 151],
               '3c7tA_core' :
                [73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 
                 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 
                 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 172, 174, 176, 180, 181, 182, 
                 183, 184, 185, 186, 187, 188, 189, 233, 234, 235, 236, 237, 238, 239, 241, 242, 243, 244, 245, 246, 247, 
                 248, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 
                 276, 277, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 
                 309, 310, 311, 314, 315, 316],
               '2a6pA_core' :
                [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 
                 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 
                 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 122, 123, 124, 125, 126, 127, 
                 128, 129, 130, 131, 132, 133, 134, 135, 136, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 
                 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 167, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 
                 179, 180, 181, 184, 185, 186, 187, 188, 189, 190, 191, 193, 194, 195] }
        
        self.assertEqual(got,ref)
        shutil.rmtree(work_dir)
        return

#
# Run unit tests
if __name__ == "__main__":

    if True:
        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        formatter = logging.Formatter('%(message)s')
        ch.setFormatter(formatter)
        root.addHandler(ch)

    unittest.main(verbosity=2)

