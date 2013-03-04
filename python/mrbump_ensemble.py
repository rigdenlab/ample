'''
Created on Feb 28, 2013

@author: jmht
'''

# python imports
import logging
import os
import re
import shutil
import sys
import unittest

# our imports
import clusterize
import Final_display_results
import printTable
import run_mr_bump_shelx_parallel


def mrbump_ensemble_cluster( ensembles, amopt, clusterID="X" ):
    """
    Process the list of ensembles using MrBump on a cluster
    INPUTS:
    ensembles - list of pdb files, each file containing a group of ensembles
    amopt - amopt object
    clusterId: number/id of the cluster that is being processed
    """
    logger = logging.getLogger()
    logger.info("Running MR and model building on a cluster\n\n")

    mrBuild = clusterize.ClusterRun()
    mrBuild.QTYPE = "SGE"

    mrBuildClusterDir = os.path.join(amopt.d['mrbump_dir'], "cluster_run" + str(clusterID))
    os.mkdir(mrBuildClusterDir)
    #mrBuild.getMTZInfo(amopt.d['mtz'], mrBuildClusterDir)

    # Reset the queue list
    mrBuild.qList = []
    jobID = 0
    for pdbfile in ensembles:
        mrBuild.mrBuildOnCluster(mrBuildClusterDir, pdbfile, jobID, amopt)
        jobID = jobID + 1
    mrBuild.monitorQueue()
    
    work_dir = amopt.d['work_dir']
    shutil.rmtree(work_dir + '/fine_cluster_' + str(clusterID))
    # shutil.rmtree(work_dir+'/pre_models')
    for l in os.listdir(work_dir + '/spicker_run'):
        if os.path.splitext(l)[1] == 'pdb':
            os.remove(work_dir + '/spicker_run/' + l)
    os.remove(work_dir + '/spicker_run/rep1.tra1')
    T = printTable.Table()

    T.bumppath = mrBuildClusterDir
    T.cluster = True
    table = T.maketable()
    out = sys.stdout
    T.pprint_table(out, table)

        # cleanup
    # for each_run in os.listdir(mrBuildClusterDir ):
        #   if os.path.isdir(  os.path.join(mrBuildClusterDir, each_run)):
        #      name=re.split('_', each_run)
        #      mrBuildOutputDir=os.path.join(bump_dir, "cluster_run"+str(cluster)+"result"+name[1])
        #      os.mkdir(mrBuildOutputDir)
        #      shutil.move (os.path.join(mrBuildClusterDir, each_run, "search_"+name[1]+"_mrbump","phaser_shelx" ),mrBuildOutputDir  )
        #      shutil.move (os.path.join(mrBuildClusterDir, each_run, "search_"+name[1]+"_mrbump","molrep_shelx" ),mrBuildOutputDir  )
        #      shutil.move (os.path.join(mrBuildClusterDir, each_run, "search_"+name[1]+"_mrbump","data" ),mrBuildOutputDir  )
        #      shutil.move (os.path.join(mrBuildClusterDir, each_run, "logs" ),mrBuildOutputDir  )
        # shutil.rmtree(mrBuildClusterDir)
    
    # Monitor the cluster queue to see when all jobs have finished
##End mrbump_ensemble_cluster

def mrbump_ensemble_local( ensembles, amopt, clusterID="X" ):
    """
    Process the list of ensembles using MrBump on a local machine
    INPUTS:
    ensembles - list of pdb files, each file containing a group of ensembles
    amopt - amopt object
    clusterId: number/id of the cluster that is being processed
    """
    # Split ensembles and rum mrbump
    split_ensembles = run_mr_bump_shelx_parallel.split(ensembles, amopt.d['nproc'])
    run_mr_bump_shelx_parallel.split_into_runs(split_ensembles, amopt)

    bump_dir = amopt.d['mrbump_dir']
    # Process results
    if amopt.d['use_shelxe']:
        Final_display_results.make_log(bump_dir, os.path.join(amopt.d['work_dir'], 'Final_results.log'))
        # print '\n\nFinal Results:\n\n'
        # T=printTable.Table()
        # T.bumppath = work_dir +'/MRBUMP_cluster'+str(cluster)
        # T.cluster = False
        # table = T.maketable()
        # out = sys.stdout
        # T.pprint_table(out, table)
    else:
        resultslog = open(amopt.d['work_dir'] + os.sep + 'Results.log', "w")
        print 'getting results from: {}',format( bump_dir )
        for mrbumplog in os.listdir(bump_dir):
            if re.search('.log', mrbumplog):
                # print mrbumplog
                for line in open(mrbumplog):
                    if re.search('^(\d)\s*loc0_', line):
                        if not re.search('method', line):
                            print line
                            resultslog.write(line)
        resultslog.close()
        
##End mrbump_ensemble_local

class Test(unittest.TestCase):


    def testName(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()