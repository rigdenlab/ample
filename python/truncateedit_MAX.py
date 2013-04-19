#!/usr/bin/python2.6

# python imports
import copy
import glob
import logging
import os
import re
import shutil
import subprocess
import sys
import time

# Our imports
import ample_util
import cluster_with_MAX
import pdb_edit

def One_model_only(list_of_ensembles, rundir):
    if not os.path.exists(rundir + '/Top_model_ensembles'):
        os.mkdir(rundir + '/Top_model_ensembles')

    for pdb in list_of_ensembles:
        if os.path.exists(pdb):
            name = os.path.split(pdb)
            #print name
            #print pdb
            outpdb = open(rundir + '/Top_model_ensembles/'+name[-1], "w")

            for line in open(pdb):
                if re.search('ENDMDL', line):
                    break
                outpdb.write(line)
            outpdb.close()

    outlist = []
    for outpdb in os.listdir(rundir + '/Top_model_ensembles'):
        outlist.append( rundir + '/Top_model_ensembles/'+ outpdb)
    return  outlist

#############cluster_with_MAX####
def make_list_to_keep( theseus_out, THESEUS_threthold ):
    """
    Make a list of residues to keep under variance threshold
    INPUTS:
    theseus_out: theseus_variances.txt output file
    THESEUS_threthold: the threshold variance

    RETURNS:
    a list of the residue indexes
    """
    add_list =[]

    theseus_out = open(theseus_out)
    for line in theseus_out:
        # for alternate versions of theseus remove RES cards
        line = re.sub('RES\s*', '', line)
        pattern = re.compile('^(\d*)\s*(\w*)\s*(\d*)\s*(\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)')
        result = pattern.match(line)
        if result:
            seq = re.split(pattern, line)
            #print seq
            if (seq[6]) != 'T':

                if float(seq[4]) <= float(THESEUS_threthold):
                    if seq[1] != '':
                #  print 'GOT ' + str(seq[4]) + '  thresh ' + str(THESEUS_threthold) + ' KEEP ' + str(seq[3])
                        add_list.append( int(seq[3]) )

    #print add_list
    return add_list
#END make_list_to_keep

def chunks(a_list, percent_interval):


    for i in xrange(0, len(a_list), percent_interval):
        yield a_list[i:i+percent_interval ]
##############################
def fly_threshold(theseus_out, percent):
    """
    Make a list of residues to keep under variance threshold
    Threshold is chosen based on number of residues kept

    Args:
    theseus_out: theseus_variances.txt output file
    percent: % of resdiues to have under each truncation level
    
    Return:
    list of the thresholds to truncate at
    """
    # List of variances ordered by residue index
    var_list=[]

    # print theseus_out
    theseus_out = open(theseus_out)
    
    #jmht - build up a list of the variances of each residue
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

    length = len(var_list)
    #print length
    #try 10 percent
    percent_interval=(float(length)/100) *float(percent)
    #print int(percent_interval)

    #lowest=min(var_list)
    #highest=max(var_list)
    # print lowest, highest

    ## try to find intervals for truncation
    try_list=copy.deepcopy(var_list)
    try_list.sort()

    Thresholds = []
    # print list(chunks(try_list, int(percent_interval)))

    for x in list(chunks(try_list, int(percent_interval))):
        # print x[-1]
        Thresholds.append(x[-1])

    return  Thresholds
###End fly_threshold

################################
def Run_Rosetta(string, no_of_files, radius, Rosetta_cluster, RDB):  # rosetta can fail so repeat job
    rosettaout = os.getcwd()
    condition = 0

    while condition == 0:
        os.system (Rosetta_cluster+ ' -database '+RDB+' -cluster:radius ' +str(radius) + ' -in:file:s '+string+' > collect_out')

        my_sdtout = open (rosettaout + '/collect_out')
        for line in my_sdtout:
            file_pattern = re.compile('^Clustering\s*(\d*)\s*structures')
            file_result = file_pattern.match(line)
            if file_result:
        #print line

                file_result2 = re.split(file_pattern, line )
                #print file_result2
                if int(file_result2[1]) !=  no_of_files:
                    print 'fail there are ' + str(no_of_files) + ' files, got ' + file_result2[1]
                if int(file_result2[1]) ==  no_of_files:
                    print ' SUCCESS there are ' + str(no_of_files) + ' files, got ' + file_result2[1]
                    condition = 1
    return

def make_ensembles(trunc_out, threshold, theseus_exe, maxcluster_exe ):
    """
    Given a directory of truncated PDB files, use maxcluster to cluster them
    according to the three radius thresholds in RADS

    If there are more than 2 clusters returned by maxcluster used theseus to align them

    INPUTS:
    trunc_out: directory with truncated PDB files (e.g. fine_cluster_2/trunc_files_2)
    threshold: a float with the threshold the files were truncated with - used for file naming
    THESUS: path to theseus executable
    maxcluster_exe: path to maxcluster executable

    """

    ensembles_made = []

    RADS = [ 1, 2, 3 ] # radius thresholds
        
    # We make a list of all the files truncated at this level and use
    # this as input to all the programs
    file_list = glob.glob( os.path.join( trunc_out, '*.pdb' ) )
    
    # Create file holdings a list  of all files (needed by maxcluster)
    fname = os.path.join( trunc_out, "files.list" )
    f = open( fname, 'w' )
    f.write( "\n".join( file_list )+"\n" )
    f.close()

    for RAD in RADS:
        
        logging.debug("Clustering files under radius: {0}".format(RAD) )
        
        ensemble_dir = os.path.join( trunc_out, 'fine_clusters_'+str(RAD)+'_ensemble' )
        os.mkdir( ensemble_dir )
        os.chdir( ensemble_dir )

        # A list of the files that have been clustered together with maxcluster
        cluster_files = cluster_with_MAX.cluster_with_MAX_FAST( fname, RAD, maxcluster_exe )  # use fastest method
        logging.debug("Maxcluster clustered {0} files".format ( len( cluster_files ) ) )
        if cluster_files < 2:
            logging.info( 'Could not create ensemble for radius {0} (models too diverse)'.format( RAD ) )
            continue
        
        # Write out the files
        file_list = "maxcluster_radius_{0}_files.list".format( RAD )
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
        basename='trunc_{0}_rad_{1}'.format( threshold, RAD ) 

        # Run theseus to generate a file containing the aligned clusters
        cmd = [ theseus_exe, "-r", basename, "-a0" ] + cluster_files
        retcode = ample_util.run_command( cmd, logfile=basename+"_theseus.log" )
        
        # jmht - the previous Align_rosetta_fine_clusters_with_theseus routine was running theseus twice and adding the averaged
        # ensemble from the first run to the ensemble. This seems to improve the results for the TOXD test case - maybe something to
        # look at?
        
        # Check if theseus worked and if so rename the file with the aligned files and append the path to the ensembles
        cluster_file = os.path.join(  ensemble_dir, basename+'_sup.pdb' )
        ensemble_file = os.path.join( ensemble_dir, basename+'.pdb' )

        #same from here     
        if os.path.exists( cluster_file ):
            shutil.move( cluster_file, ensemble_file )
            ensembles_made.append( ensemble_file )
        else:
            logging.info( 'Could not create ensemble for files (models too diverse): {0}'.format( ensemble_file ) )
            
    #print ensembles_made

    return ensembles_made
###END make_ensembles

def truncate( theseus_exe, models_list, work_dir, percent, FIXED_INTERVALS=False ):
    """
    Truncate the models in one folder
    * Run theseus to find the variances
    * For each truncation level, create a directory containing the truncated PDB files

    INPUTS:
    theseus_exe: path to theseus
    models_list: list of paths to the PDBs to be truncated
    work_dir: RunDir+'/fine_cluster_'+str(samples) [ root directory where truncated clusters will go ]
    percent: (var_args['percent'] or 5 ) - mutually exclusive with FIXED_INTERVALS. Truncate models so that "percent"
             residues are kept each cycle
    FIXED_INTERVALS: use the hard-coded intervals for truncation: [1, 1.5, 2 , 2.5, 3, 3.5 ,4, 4.5, 5, 5.5, 6, 7 , 8]

    OUTPUTS:
    truncate_log: file listing how many residues kept under each truncation threshold
    List of PDB files in trunc_out (e.g. fine_cluster_2/trunc_files_2.5) that contain the PDBs truncated to this level
    """

    work_dir = os.path.abspath( work_dir )
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
    os.chdir(work_dir)

    #--------------------------------
    # get variations between pdbs
    #--------------------------------
    cmd = [ theseus_exe, "-a0" ] + models_list
    retcode = ample_util.run_command( cmd, logfile="theseus.log" )

    #--------------------------------
    # choose threshold type
    #-------------------------------
    T_data = os.path.join( work_dir, 'theseus_variances.txt' )
    if FIXED_INTERVALS:
        thresholds = [ 1, 1.5, 2 , 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8 ]
    else:
        thresholds = fly_threshold(T_data, percent)
    
    logging.debug("Got {0} thresholds: {1}".format( len(thresholds), thresholds ))

    #-------------------------------
    #truncate
    #----------------------------------
    truncate_log = open(work_dir + '/truncate.log', "w")
    truncate_log.write('This is the number of residues kept under each truncation threshold\n\nthreshold\tnumber of residues\n')
    truncation_result = []
    for threshold in thresholds:

        # Get a list of the indexes of the residues to keep
        add_list = make_list_to_keep(T_data, threshold)
        
        truncate_log.write( str(threshold) +'\t' + str(len(add_list)) + '\n' )
        logging.info( 'Keeping: {0} residues at truncation level: {1}'.format( len(add_list), threshold ) )
        logging.debug( 'The following residues will be kept: {0}'.format( add_list ) )
        
        if len(add_list) < 1:
            #######   limit to size of pdb to keep, if files are too variable, all residues are removed
            logging.debug( 'No residues kept at this truncation level.' )
            continue

        trunc_out = os.path.join( work_dir, 'trunc_files_' + str(threshold) )
        os.mkdir(trunc_out)
        logging.info( 'truncating at: {0} in directory {1}'.format( threshold, trunc_out ) )
        
        pdbed = pdb_edit.PDBEdit()
        for infile in models_list:
            pdbname = os.path.split( infile )[1]
            pdbout = os.path.join( trunc_out, pdbname )

            # Loop through PDB files and create new ones that only contain the residues left after truncation
            pdbed.select_residues( inpath=infile, outpath=pdbout, residues=add_list )
        
        truncation_result.append( ( threshold, trunc_out ) )
            
    return truncation_result

###END truncate

def truncate_Phenix(PHENIX, models_list, work_dir ):
    print 'asembling with Phenix'
    print PHENIX
    phenix_name = PHENIX

    if not re.search('ensembler', PHENIX):
        phenix_name = PHENIX+'.ensembler'
    print phenix_name

    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
        
    os.chdir(work_dir)

    number_of_models = 0
    list_of_pdbs = []
    string = ''
    for infile in models_list:
        number_of_models +=1
        string  = string + infile + ' '
        list_of_pdbs.append(infile)

    cmd = phenix_name+' '+string
    #print cmd
    os.system(cmd)

    return[ os.path.join(work_dir,'ensemble_merged.pdb') ]

###END truncate_Phenix

###################
#/home/jaclyn/LARGE_RUN/1EN2/fine/THESEUS_polyala_cluster0_trunc_4/ALIGNED_rad_3/theseus
#/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0

if __name__ =="__main__":

#    #            /home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/cluster_60/sorted_cluster_0
#    THESEUS =  '/home/jaclyn/LARGE_RUN/1EN2/fine/THESEUS_polyala_cluster0_trunc_4/ALIGNED_rad_3/theseus'
#    models_path='/home/jaclyn/BEE/new_case/olga/models'
#    work_dir ='/home/jaclyn/DOMAINS/testing/truncted'
#    MAX = '/home/jaclyn/programs/maxcluster/maxcluster'
#    percent = 5
#
#    #list_of_ensembles = truncate_Phenix(THESEUS, models_path, work_dir, MAX, percent, True )
#
#    rundir = '/home/jaclyn/BEE/ample_test/RO'
#    list_of_ensembles = ['/home/jaclyn/BEE/ample_test/MODEL.pdb', '/home/jaclyn/BEE/ample_test/MODEL (copy).pdb']
#
#    One_model_only(list_of_ensembles, rundir)
    
    
        
    logging.basicConfig()
    logging.getLogger().setLevel(logging.DEBUG)
    
    trunc_dir = "fine_cluster_X"
    spicker_cluster = ["/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/models/1_S_00000003.pdb",
                        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/models/1_S_00000005.pdb",
                        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/models/2_S_00000002.pdb",
                        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/models/2_S_00000005.pdb",
                        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/models/2_S_00000006.pdb",
                        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/models/3_S_00000001.pdb",
                        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/models/3_S_00000002.pdb",
                        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/models/3_S_00000004.pdb",
                        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/models/3_S_00000005.pdb",
                        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/models/4_S_00000004.pdb",
                        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/models/5_S_00000002.pdb",
                        "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_0/models/5_S_00000003.pdb" ]
    
    theseus_exe = "/opt/theseus_src/theseus"
    work_dir = trunc_dir
    maxcluster_exe = "/opt/maxcluster/maxcluster"
    percent = 50
    
    #print " ".join( spicker_cluster )
    trunc_list = truncate( theseus_exe, spicker_cluster, work_dir, percent )
    
    list_of_ensembles = []
    for threshold, tdir in trunc_list:
        list_of_ensembles += make_ensembles( tdir,
                                            threshold,
                                            theseus_exe,
                                            maxcluster_exe )
    
    