'''
Created on Apr 18, 2013

@author: jmht
'''

import copy
import glob
import logging
import os
import shutil
import sys
import unittest

# our imports
import ample_util
import pdb_edit
import run_spicker
import subcluster

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
        
        # ensembles
        self.ensembles_directory=None
        self.ensembles=None
        self.ensembles_data=None
        
        # misc
        self.logger=logging.getLogger()
        
        return
    
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
    
    def _calculate_residues_thresh(self,var_by_res,percent_interval):
        """Txxx
        """

        # calculate the thresholds
        truncation_variances=self.generate_thresholds(var_by_res,percent_interval)

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
            self.logger.info( '* Running SPICKER to cluster models *' )
            spicker_rundir = os.path.join( self.work_dir, 'spicker')
            spickerer = run_spicker.SpickerCluster(run_dir=spicker_rundir,
                                                   spicker_exe=cluster_exe,
                                                   models=models,
                                                   num_clusters=num_clusters )
            spickerer.run_spicker()
            self.logger.debug( spickerer.results_summary() )
            
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
        d['num_residues'] = None

        # subclustering info
        d['num_models'] = None
        d['subcluster_radius_threshold'] = None
        d['subcluster_centroid_model'] = None
    
        # ensemble info
        d['name'] = None
        d['side_chain_treatment'] = None
        d['num_atoms'] = None
        d['pdb'] = None # path to the ensemble file
        
        return d
    
    def edit_side_chains(self,raw_ensemble,raw_ensemble_data,ensembles_directory):
        
        assert os.path.isdir(ensembles_directory),"Cannot find ensembles directory: {0}".format(ensembles_directory)
        ensembles=[]
        ensembles_data=[]
        for sct in self.side_chain_treatments:
            
            # create filename based on side chain treatment
            fpath = ample_util.filename_append(raw_ensemble,astr=sct, directory=ensembles_directory)
            
            # Create the files
            if sct == "allatom":
                # For all atom just copy the file
                shutil.copy2(raw_ensemble,fpath)
            elif sct == "reliable":
                pdb_edit.reliable_sidechains(raw_ensemble,fpath)
            elif sct == "polya":
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
            ensemble_data['pdb']=fpath
            ensemble_data['num_atoms']=natoms
            # check
            assert ensemble_data['num_residues']==nresidues,"Unmatching number of residues!"
            
            ensembles.append(fpath)
            ensembles_data.append(ensemble_data)
                
        return ensembles,ensembles_data
  
    def generate_ensembles(self,models,
                           cluster_method=None,
                           cluster_exe=None,
                           num_clusters=None,
                           percent_truncation=None,
                           truncation_method=None,
                           ensembles_directory=None,
                           work_dir=None):
        
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
                                                              cluster_exe=cluster_exe)):
            if len(cluster) < 2:
                self.logger.info("Cannot truncate cluster {0} as < 2 models!".format(cluster_data.cluster_num))
                continue
            for truncated_models, truncated_models_data in zip(*self.truncate_models(cluster,
                                                                                     cluster_data,
                                                                                     truncation_method=truncation_method,
                                                                                     percent_truncation=percent_truncation)):
                for subcluster, subcluster_data in zip(*self.subcluster_models(truncated_models,
                                                                               truncated_models_data,
                                                                               radius_thresholds=self.subcluster_radius_thresholds,
                                                                               subcluster_program=self.subcluster_program,
                                                                               subcluster_exe=self.subcluster_program,
                                                                               ensemble_max_models=self.ensemble_max_models)):
                    for ensemble, ensemble_data in zip(*self.edit_side_chains(subcluster, subcluster_data, self.ensembles_directory)):
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
        var_list=[var for (_,var) in var_by_res]

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
        cluster_num=truncated_models_data['cluster_num']
        truncation_level=truncated_models_data['truncation_level']
        truncation_dir=truncated_models_data['truncation_dir']
            
        # Run maxcluster to generate the distance matrix
        if subcluster_program=='maxcluster':
            clusterer = subcluster.MaxClusterer(self.subcluster_exe)
        else:
            assert False
        clusterer.generate_distance_matrix(truncated_models)

        # Loop through the radius thresholds
        num_previous_models=-1 # set to -1 so comparison always false on first pass
        for radius in radius_thresholds:

            self.logger.debug("subclustering files under radius: {0}".format(radius))

            # Get list of pdbs clustered according to radius threshold
            cluster_files = clusterer.cluster_by_radius(radius)
            self.logger.debug("Clustered {0} files".format(len(cluster_files)))
            if cluster_files < 2:
                self.logger.info( 'Could not cluster files using radius {0}'.format(radius))
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
                self.logger.debug( 'Could not create ensemble for files (models too diverse): {0}'.format( basename ) )
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

            # The data we've collected is the same for all pdbs in this level so just keep using the first  
            ensemble_data=copy.copy(truncated_models_data)
            ensemble_data['num_models'] = len( cluster_files )
            ensemble_data['subcluster_radius_threshold'] = radius
            ensemble_data['pdb'] = ensemble

            # Get the centroid model name from the list of files given to theseus - we can't parse
            # the pdb file as theseus truncates the filename
            ensemble_data['subcluster_centroid_model']=os.path.basename(cluster_files[0])
            
            ensembles.append(ensemble)
            ensembles_data.append(ensemble_data)
        
        return ensembles,ensembles_data
    
    def truncate_models(self,models,models_data,truncation_method,percent_truncation):
        
        assert len(models) > 1,"Cannot truncate as < 2 models!"

        # Create the directories we'll be working in
        truncate_dir =  os.path.join(self.work_dir, 'truncate_{0}'.format(models_data['cluster_num']))
        os.mkdir(truncate_dir)
        os.chdir(truncate_dir)
        
        # Calculate variances between pdb
        var_by_res=self.calculate_variances(models)

        # Calculate which residues to keep under the different methods
        truncation_levels, truncation_variances, truncation_residues=None,None,None
        if truncation_method=='percent':
            truncation_levels, truncation_variances, truncation_residues=self._calculate_residues_percent(var_by_res,percent_truncation)
        elif truncation_method=='thresh':
            truncation_levels, truncation_variances, truncation_residues=self._calculate_residues_thresh(var_by_res,percent_truncation)
        else:
            raise RuntimeError,"Unrecognised ensembling mode: {0}".format(truncation_method)

        truncated_models=[]
        truncated_models_data=[]
        for tlevel,tvar,tresidues in zip(truncation_levels, truncation_variances, truncation_residues):
            trunc_dir = os.path.join(truncate_dir, 'tlevel_{0}'.format(tlevel))
            os.mkdir(trunc_dir)
            self.logger.info( 'truncating at: {0} in directory {1}'.format(tvar,trunc_dir))

            # list of models for this truncation level
            level_models = []
            for infile in models:
                pdbname = os.path.basename( infile )
                pdbout = os.path.join(trunc_dir, pdbname)
                # Loop through PDB files and create new ones that only contain the residues left after truncation
                pdb_edit.select_residues(inpath=infile, outpath=pdbout, residues=tresidues)
                level_models.append(pdbout)
            
            # Add the model
            truncated_models.append(level_models)

            # Add the data
            model_data=copy.copy(models_data)
            model_data['truncation_level']=tlevel
            model_data['truncation_variance']=tvar
            model_data['truncation_residues']=tresidues
            model_data['num_residues']=len(tresidues)
            model_data['truncation_dir']=trunc_dir
            model_data['percent_truncation'] = percent_truncation
            model_data['truncation_method'] = truncation_method
            
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
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
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
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        
        var_by_res=ensembler.calculate_variances(cluster_models)
        truncation_levels, truncation_variances, truncation_residues=ensembler._calculate_residues_thresh(var_by_res,percent_interval)
        
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
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        var_by_res=ensembler.calculate_variances(cluster_models)
        truncation_levels, truncation_variances, truncation_residues=ensembler._calculate_residues_percent(var_by_res,percent_interval=5)


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
            [26, 27, 30, 31, 32],
            [31, 32]]
        self.assertEqual(residues,truncation_residues)

        self.assertEqual(truncation_levels,
                         [100, 95, 90, 85, 80, 75, 69, 64, 59, 54, 49, 44, 39, 34, 29, 24, 19, 14, 8, 3])
        
        refv=[74.272253, 59.917999, 56.260987, 50.992344, 43.268754, 37.227859, 30.767268, 25.348501,
              23.050357, 18.292366, 15.287348, 10.250043, 5.66545, 3.152793, 0.874899, 0.185265, 0.090917,
            0.08155, 0.039639, 0.018215]
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
        
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        models=glob.glob(mdir+os.sep+"*.pdb")
        
        cluster_models,cluster_data=ensembler.cluster_models(models=models,
                                                             cluster_method='spicker',
                                                             num_clusters=3,
                                                             cluster_exe=self.spicker_exe)
        
        names=[os.path.basename(m) for m in cluster_models[0]]
        self.assertEqual(names,
                         ['4_S_00000002.pdb', '4_S_00000005.pdb', '5_S_00000005.pdb', '4_S_00000003.pdb',
                          '2_S_00000005.pdb', '5_S_00000004.pdb', '1_S_00000005.pdb', '3_S_00000003.pdb',
                          '3_S_00000006.pdb', '1_S_00000002.pdb', '2_S_00000001.pdb', '1_S_00000004.pdb',
                          '3_S_00000004.pdb']
                         )
        
        d = cluster_data[2]
        self.assertEqual(os.path.basename(d['cluster_centroid']),'1_S_00000001.pdb')
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
        
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
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

        eref=['c1_tl100_r2_allatom.pdb', 'c1_tl100_r2_reliable.pdb', 'c1_tl100_r2_polya.pdb', 'c1_tl100_r3_allatom.pdb',
              'c1_tl100_r3_reliable.pdb', 'c1_tl100_r3_polya.pdb', 'c1_tl95_r2_allatom.pdb', 'c1_tl95_r2_reliable.pdb',
              'c1_tl95_r2_polya.pdb', 'c1_tl95_r3_allatom.pdb', 'c1_tl95_r3_reliable.pdb', 'c1_tl95_r3_polya.pdb', 'c1_tl90_r1_allatom.pdb',
              'c1_tl90_r1_reliable.pdb', 'c1_tl90_r1_polya.pdb', 'c1_tl90_r2_allatom.pdb', 'c1_tl90_r2_reliable.pdb', 'c1_tl90_r2_polya.pdb',
              'c1_tl90_r3_allatom.pdb', 'c1_tl90_r3_reliable.pdb', 'c1_tl90_r3_polya.pdb', 'c1_tl85_r1_allatom.pdb', 'c1_tl85_r1_reliable.pdb',
              'c1_tl85_r1_polya.pdb', 'c1_tl85_r2_allatom.pdb', 'c1_tl85_r2_reliable.pdb', 'c1_tl85_r2_polya.pdb', 'c1_tl85_r3_allatom.pdb',
              'c1_tl85_r3_reliable.pdb', 'c1_tl85_r3_polya.pdb', 'c1_tl80_r1_allatom.pdb', 'c1_tl80_r1_reliable.pdb', 'c1_tl80_r1_polya.pdb',
              'c1_tl80_r2_allatom.pdb', 'c1_tl80_r2_reliable.pdb', 'c1_tl80_r2_polya.pdb', 'c1_tl80_r3_allatom.pdb', 'c1_tl80_r3_reliable.pdb',
              'c1_tl80_r3_polya.pdb', 'c1_tl75_r1_allatom.pdb', 'c1_tl75_r1_reliable.pdb', 'c1_tl75_r1_polya.pdb', 'c1_tl75_r2_allatom.pdb',
              'c1_tl75_r2_reliable.pdb', 'c1_tl75_r2_polya.pdb', 'c1_tl75_r3_allatom.pdb', 'c1_tl75_r3_reliable.pdb', 'c1_tl75_r3_polya.pdb',
              'c1_tl69_r1_allatom.pdb', 'c1_tl69_r1_reliable.pdb', 'c1_tl69_r1_polya.pdb', 'c1_tl69_r2_allatom.pdb', 'c1_tl69_r2_reliable.pdb',
              'c1_tl69_r2_polya.pdb', 'c1_tl64_r1_allatom.pdb', 'c1_tl64_r1_reliable.pdb', 'c1_tl64_r1_polya.pdb', 'c1_tl64_r2_allatom.pdb',
              'c1_tl64_r2_reliable.pdb', 'c1_tl64_r2_polya.pdb', 'c1_tl59_r1_allatom.pdb', 'c1_tl59_r1_reliable.pdb', 'c1_tl59_r1_polya.pdb',
              'c1_tl59_r2_allatom.pdb', 'c1_tl59_r2_reliable.pdb', 'c1_tl59_r2_polya.pdb', 'c1_tl54_r1_allatom.pdb', 'c1_tl54_r1_reliable.pdb',
              'c1_tl54_r1_polya.pdb', 'c1_tl54_r2_allatom.pdb', 'c1_tl54_r2_reliable.pdb', 'c1_tl54_r2_polya.pdb', 'c1_tl49_r1_allatom.pdb',
              'c1_tl49_r1_reliable.pdb', 'c1_tl49_r1_polya.pdb', 'c1_tl49_r2_allatom.pdb', 'c1_tl49_r2_reliable.pdb', 'c1_tl49_r2_polya.pdb',
              'c1_tl44_r1_allatom.pdb', 'c1_tl44_r1_reliable.pdb', 'c1_tl44_r1_polya.pdb', 'c1_tl44_r2_allatom.pdb', 'c1_tl44_r2_reliable.pdb',
              'c1_tl44_r2_polya.pdb', 'c1_tl39_r1_allatom.pdb', 'c1_tl39_r1_reliable.pdb', 'c1_tl39_r1_polya.pdb', 'c1_tl34_r1_allatom.pdb',
              'c1_tl34_r1_reliable.pdb', 'c1_tl34_r1_polya.pdb', 'c1_tl29_r1_allatom.pdb', 'c1_tl29_r1_reliable.pdb', 'c1_tl29_r1_polya.pdb', 
              'c1_tl24_r1_allatom.pdb', 'c1_tl24_r1_reliable.pdb', 'c1_tl24_r1_polya.pdb', 'c1_tl19_r1_allatom.pdb', 'c1_tl19_r1_reliable.pdb',
              'c1_tl19_r1_polya.pdb', 'c1_tl14_r1_allatom.pdb', 'c1_tl14_r1_reliable.pdb', 'c1_tl14_r1_polya.pdb', 'c1_tl8_r1_allatom.pdb',
              'c1_tl8_r1_reliable.pdb', 'c1_tl8_r1_polya.pdb']
        self.assertEqual([os.path.basename(m) for m in ensembles],eref)
        d = ensembler.ensembles_data[5]

        self.assertEqual(d['percent_truncation'],percent_truncation)
        self.assertEqual(d['truncation_method'],truncation_method)
        self.assertEqual(d['cluster_method'],cluster_method)
        self.assertEqual(d['num_clusters'],num_clusters)
        self.assertEqual(d['truncation_variance'],13.035172)
        self.assertEqual(d['num_atoms'],290)
        self.assertEqual(d['subcluster_centroid_model'],'4_S_00000002.pdb')
        
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
        
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
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
        
        self.assertEqual(d['truncation_variance'],27.389253)
        self.assertEqual(d['percent_truncation'],percent_truncation)
        self.assertEqual(d['truncation_method'],truncation_method)
        self.assertEqual(d['cluster_method'],cluster_method)
        self.assertEqual(d['num_clusters'],num_clusters)
        self.assertEqual(d['num_atoms'],290)
        self.assertEqual(d['side_chain_treatment'],'polya')
        self.assertEqual(d['subcluster_centroid_model'],'4_S_00000002.pdb')
        
        shutil.rmtree(ensembler.work_dir)
        return

def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testThresholds'))
#     suite.addTest(Test('testResiduesThresh'))
#     suite.addTest(Test('testResiduesPercent'))
#     suite.addTest(Test('testClustering'))
#     suite.addTest(Test('testEnsemblingThresh'))
#     suite.addTest(Test('testEnsemblingPercent'))
    return suite
    
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

    unittest.TextTestRunner(verbosity=2).run(testSuite())

