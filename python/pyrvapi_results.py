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

try: import pyrvapi
except: pyrvapi=None

class MyClass(object):
    '''
    classdocs
    '''


    def __init__(self, params):
        '''
        Constructor
        '''


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
pyrvapi.rvapi_add_tab("tab1","Summary",False)
pyrvapi.rvapi_add_section("sec1","Results","tab1",0,0,1,1,True)
pyrvapi.rvapi_add_table1 ("sec1/table1","Ensembling Results",1,0,1,1,True)

# Get the ensembling data
ensembles_data=results_dict['ensembles_data']
clusters, cluster_method, truncation_method, percent_truncation = ensemble.collate_cluster_data(ensembles_data)
cluster_num=1
tdata = ensemble.cluster_table_data(clusters, cluster_num)

# Make column headers
for i in range(1,len(tdata[0])-1): # Skip name as it's the row header
    #print pyrvapi.rvapi_put_horz_theader( "table1","Row 1","Tooltip 1",0 ),
    print pyrvapi.rvapi_put_horz_theader("table1",tdata[0][i],"",i-1),

# Make ensemble name row headers
for i in range(1,len(tdata)-1):
    print pyrvapi.rvapi_put_vert_theader("table1",tdata[i][0],"",i-1),
    
for i in range(1,len(tdata)-1):
    for j in range(1,len(tdata[i])-1):
        print pyrvapi.rvapi_put_table_string ("table1",str(tdata[i][j]),i-1,j-1),

print pyrvapi.rvapi_flush(),
