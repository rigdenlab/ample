'''
Created on Apr 18, 2013

@author: jmht
'''

import copy
import glob
import logging
import os
import re
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
        self.truncation_threshold = None

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
        self.percent= None

        # theseus variance file
        self.variance_log = None

        # list of tupes: (resSeq,variances)
        self.var_by_res=[]

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

        # The maximum number of models in an ensemble
        self.max_ensemble_models = None

        # Programs we use
        self.theseus_exe = None

        self.maxcluster_exe = None

        return

    def pdb_list(self):
        """Return the final ensembles as a list of pdb files"""
        return [ensemble.pdb for ensemble in self.ensembles]

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

        return

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
        retcode = ample_util.run_command(cmd,
                                         logfile=os.path.join(self.work_dir,"theseus.log"),
                                         directory=self.work_dir)
        if retcode != 0:
            msg = "non-zero return code for theseus in generate_thresholds!"
            logging.critical( msg )
            raise RuntimeError, msg

        #build up a list of the variances of each residue
        self.variance_log=os.path.join(self.work_dir,'theseus_variances.txt')
        self.var_by_res=self.read_variances(self.variance_log)

        # List of variances ordered by residue index
        var_list=[var for (resSeq,var) in self.var_by_res]

        length = len(var_list)
        if length == 0:
            msg = "Error generating thresholds, got len: {0}".format(length)
            logging.critical(msg)
            raise RuntimeError,msg

        # How many residues should fit in each bin
        chunk_size=int( ( float(length)/100 ) *float(self.percent) )
        if chunk_size < 1:
            msg = "Error generating thresholds, got < 1 AA in chunk_size"
            logging.critical(msg)
            raise RuntimeError,msg
#        lowest=min(var_list)
#        highest=max(var_list)
#        print length
#        print lowest, highest
#        print int( chunk_size )

        ## try to find intervals for truncation
        self.truncation_thresholds=self._generate_thresholds(var_list, chunk_size)

        logging.debug("Got {0} thresholds: {1}".format( len(self.truncation_thresholds), self.truncation_thresholds ))

        return

    def _generate_thresholds(self,values,chunk_size):
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
                thresholds[-1]=variance
                counts[-1]+=count

        return thresholds

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

            # Run maxcluster to generate the distance matrix
            clusterer = cluster_with_MAX.MaxClusterer( self.maxcluster_exe )
            clusterer.generate_distance_matrix( self.truncated_models[tcount] )

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
                basename='trunc_{0}_rad_{1}'.format( truncation_threshold, radius )

                # Check if there are the same number of models in this ensemble as the previous one - if so
                # the ensembles will be identical and we can skip this one
                if num_previous_models == len( cluster_files ):
                    logging.debug( 'Number of decoys in cluster ({0}) is the same as under previous threshold so excluding cluster {1}'.format( len( cluster_files ), basename ) )
                    continue
                else:
                    num_previous_models = len( cluster_files )

                # Got files so create the directories
                ensemble_dir = os.path.join( truncation_dir, 'fine_clusters_'+str(radius)+'_ensemble' )
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
                ensemble.residues = self.truncation_residues[ tcount ]
                ensemble.num_residues = len(ensemble.residues)
                assert ensemble.num_residues==nresidues
                ensemble.num_atoms=natoms
                ensemble.truncation_threshold = truncation_threshold
                ensemble.num_truncated_models = len( self.truncated_models[ tcount ] )
                ensemble.radius_threshold = radius
                ensemble.pdb = ensemble_file

                # Get the centroid model name from the list of files given to theseus - we can't parse
                # the pdb file as theseus truncates the filename
                ensemble.centroid_model=os.path.splitext( os.path.basename(cluster_files[0]) )[0]

                self.truncated_ensembles.append( ensemble )

        return

    def read_variances(self,variance_log):
        variances=[]
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

        assert len(variances)>0,"Failed to read variances from: {0}".format(variance_log)
        return variances

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
            to_keep=[resSeq for resSeq,variance in self.var_by_res if variance <= threshold]
            self.truncation_residues.append(to_keep)

            truncate_log.write( str(threshold) +'\t' + str(len(to_keep)) + '\n' )
            logging.info( 'Keeping: {0} residues at truncation level: {1}'.format( len(to_keep), threshold ) )
            logging.debug( 'The following residues will be kept: {0}'.format( to_keep ) )

            if len(to_keep) < 1:
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
                pdbed.select_residues( inpath=infile, outpath=pdbout, residues=to_keep )

                truncated_models.append( pdbout )

            self.truncated_models.append( truncated_models )

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
        """xxx"""

        ensembler=Ensembler()

        ensembler.work_dir=os.path.join(os.getcwd(),"genthresh")
        os.mkdir(ensembler.work_dir)
        ensembler.theseus_exe=self.theseus_exe
        ensembler.percent=5
        mdir=os.path.join(self.ample_dir,"examples","toxd-example","models")
        ensembler.cluster_models=glob.glob(mdir+os.sep+"*.pdb")
        ensembler.generate_thresholds()

        self.assertEqual(29,len(ensembler.truncation_thresholds))
        self.assertEqual([0.030966, 0.070663, 0.084085, 0.099076, 0.185523, 0.723045,
                          1.79753, 3.266976, 4.886417, 6.478746, 10.378687, 11.229685,
                          16.701308, 18.823698, 22.837544, 23.741723, 25.736834, 27.761794,
                          36.255004, 37.780944, 42.318928, 49.650732, 51.481748, 56.60685,
                          58.865232, 60.570085, 70.749036, 81.15141, 87.045704],
                         ensembler.truncation_thresholds,
                         "Wrong truncation thresholds"
                         )

        truncation_residues=[]
        for threshold in ensembler.truncation_thresholds:
            to_keep=[resSeq for resSeq,variance in ensembler.var_by_res if variance <= threshold]
            truncation_residues.append(to_keep)

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

        self.assertEqual(r,truncation_residues)

        shutil.rmtree(ensembler.work_dir)
        return


    def testGenThresh(self):
        l=[0,5,6,7,9,1,2,2,2,3,3,8,8,8,8,8,8,8,8,8,8,8,3,3,3,3,3,4,5,9]
        l=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]

        return

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
    unittest.main()

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
