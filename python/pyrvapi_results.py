'''
Created on 3 Mar 2015

@author: jmht

'''

import cPickle
import os
import subprocess
import sys

sys.path.insert(0,"/opt/ample-dev1/python")
import ensemble
import mrbump_results

import pyrvapi
# try: import pyrvapi
# except: pyrvapi=None

class MyClass(object):
    '''
    classdocs
    '''


    def __init__(self, params):
        '''
        Constructor
        '''


def create_table(pyrvapi, table_id, tdata):
    # Make column headers
    for i in range(len(tdata[0])): # Skip name as it's the row header
        pyrvapi.rvapi_put_horz_theader(table_id, tdata[0][i], "", i) # Add table data
    
    for i in range(1, len(tdata) - 1):
        for j in range(len(tdata[i])):
            pyrvapi.rvapi_put_table_string(table_id, str(tdata[i][j]), i - 1, j)
    
    # Now colour the ensemble name cells
    for i in range(len(tdata) - 1):
        pyrvapi.rvapi_shape_table_cell(table_id, #tableId
            i, #row
            0, #column
            "", #tooltip
            "", # cell_css
            "table-blue-vh", #cell_style
            1, # rowSpan
            1) # colSpan
    
    pyrvapi.rvapi_flush()
    return

# Infrastructure to run
ccp4 = os.environ["CCP4"]
share_jsrview = os.path.join(ccp4, "share", "jsrview")
jsrview = os.path.join(ccp4, "libexec", "jsrview")

pklfile="/opt/ample-dev1/tests/testfiles/resultsd1.pkl"
with open(pklfile) as f: results_dict=cPickle.load(f)
dir_test=os.getcwd()

pyrvapi.rvapi_init_document ( "TestRun",dir_test,"RVAPI Demo 1",1,7,share_jsrview,None,None,None)
subprocess.Popen([jsrview, os.path.join("index.html")])

pyrvapi.rvapi_add_header("AMPLE Results")
pyrvapi.rvapi_add_tab("tab1","Summary",True) # Last arg is "open" - i.e. show or hide

#
# Ensemble Results
#
pyrvapi.rvapi_add_section("ensembles","Ensembles","tab1",0,0,1,1,False)

# Get the ensembling data
ensembles_data=results_dict['ensembles_data']
clusters, cluster_method, truncation_method, percent_truncation = ensemble.collate_cluster_data(ensembles_data)

rstr = ""
rstr += "Ensemble Results<br/>"
rstr += "----------------<br/><br/>"
rstr += "Cluster method: {0}<br/>".format(cluster_method)
rstr += "Truncation method: {0}<br/>".format(truncation_method)
rstr += "Percent truncation: {0}<br/>".format(percent_truncation)
rstr += "Number of clusters: {0}<br/><br/>".format(len(clusters.keys()))
rstr += "Generated {0} ensembles<br/><br/>".format(len(ensembles_data))
pyrvapi.rvapi_add_text(rstr,"ensembles",0,0,1,1 )
#

pyrvapi.rvapi_add_table1("ensembles/ensemble_table","Ensembling Results",1,0,1,1,True)
# for cluster_num in sorted(clusters.keys()):
#     rstr += "\n"
#     rstr += "Cluster {0}\n".format(cluster_num)
#     rstr += "Number of models: {0}\n".format(clusters[cluster_num]['cluster_num_models'])
#     rstr += "Cluster centroid: {0}\n".format(clusters[cluster_num]['cluster_centroid'])
#     rstr += "\n"
#     tdata = cluster_table_data(clusters, cluster_num)
#     rstr += tableFormat.pprint_table(tdata)        
# 
cluster_num=1
tdata = ensemble.cluster_table_data(clusters, cluster_num)
create_table(pyrvapi,"ensemble_table", tdata)

#
# MRBUMP Results
#
pyrvapi.rvapi_add_section("mrbump","MRBUMP","tab1",0,0,1,1,True)
pyrvapi.rvapi_add_table1("mrbump/mrbump_table","MRBUMP Results",1,0,1,1,True)
tdata = mrbump_results.ResultsSummary().results_table(results_dict['mrbump_results'])
create_table(pyrvapi,"mrbump_table", tdata)
