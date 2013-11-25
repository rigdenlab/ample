#!/usr/bin/env python


import cPickle
import os
import sys

sys.path.append("/Users/jmht/Documents/AMPLE/ample-dev1/python")
sys.path.append("/Users/jmht/Documents/AMPLE/ample-dev1/scripts")
#sys.path.append("/opt/ample-dev1/python")
#sys.path.append("/opt/ample-dev1/scripts")

rundir = "/Users/jmht/Documents/AMPLE"
os.chdir( rundir )


from analyse_run import AmpleResult
from contacts import Contacts
import dssp


CLUSTERNUM=0
dataDir = "/Users/jmht/Documents/AMPLE/data/coiled-coils/data"

pfile = "/Users/jmht/Documents/AMPLE/data/coiled-coils/ar_results.pkl"
pfile = "/Users/jmht/Documents/AMPLE/data/coiled-coils/test.pkl"
f = open( pfile )
resultsDict = cPickle.load( f  )
f.close()

# cPickle.dump( obj, f) 


#print resultsDict
#print len( resultsDict )

# all = []
# for result in resultsDict:
#     if result.goodContacts > 40:
#         all.append( result )
# 
# pfile = "/Users/jmht/Documents/AMPLE/data/coiled-coils/test.pkl"
# f = open( pfile, 'w' )
# cPickle.dump( all, f )
# f.close()


for result in [ resultsDict[ 0] ]: 
    
    c = Contacts()
    c.best = result.contactData
    
    dssp_file = os.path.join( dataDir, result.pdbCode+".dssp")
    
    dsspP = dssp.DsspParser( dssp_file )
    
    print c.helixFromContacts( dsspP )