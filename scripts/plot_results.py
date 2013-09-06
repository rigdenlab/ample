'''
Created on 4 Sep 2013

@author: jmht
'''

import sys
sys.path.append("/opt/ample-dev1/scripts")

import cPickle
import csv
import os
from analyse_run import AmpleResult

f = "/home/jmht/Documents/test/ar_results.pkl"
fh = open(f)
allResults = cPickle.load( fh )
fh.close()


TMdir = "/media/data/shared/TM"
rundir="/home/jmht/Dropbox/MRes/Project/Thesis/Analysis"



pdata = [ ['phaserLLG','phaserTFZ','CC','success'] ]
count=0
for r in allResults:
    success=0
    if r.shelxeCC > 30 and r.shelxeAvgChainLength > 10:
        jsuccess=1
    
    if r.mrProgram == 'phaser' and r.phaserLLG != None and r.phaserTFZ != None:
        count+=1
        pdata.append( [ r.phaserLLG, r.phaserTFZ, r.shelxeCC, success  ]  )
        
print "GOT COUNT",count



# molrepSuccess=0
# phaserSuccess=0
# molrepFail=0
# phaserFail=0
# 
# for r in allResults:
#     if r.shelxeCC > 30 and r.shelxeAvgChainLength > 10:
#         if r.mrProgram == 'phaser':
#             phaserSuccess += 1
#         else:
#             molrepSuccess += 1
#     else:
#         if r.mrProgram == 'phaser':
#             phaserFail += 1
#         else:
#             molrepFail += 1
#             
# print "Molrep: {0}:{1}".format( molrepSuccess, molrepFail)
# print "Phaser: {0}:{1}".format( phaserSuccess, phaserFail)
# print "TOTAL ",molrepSuccess+molrepFail+phaserSuccess+phaserFail
# sys.exit()

# CC = [ r.shelxeCC for r in allResults ]
# cpath = os.path.join( rundir, 'plot.csv' )
# csvfile =  open( cpath, 'wb')
# csvwriter = csv.writer(csvfile, delimiter=',',
#                         quotechar='"', quoting=csv.QUOTE_MINIMAL)
# 
# csvwriter.writerow( ["CC"] )
# for p in CC:
#     if not p:
#         p = -1
#     csvwriter.writerow( [p] )
#     
# csvfile.close()

# duff=99999
# failed = []
# success = []
# pdata = [ ['reforiginRmsd','ensembleNativeMaxsub','success'] ]
# for r in allResults:
#     jsuccess=None
#     if r.shelxeCC > 30 and r.shelxeAvgChainLength > 10:
#         jsuccess=1
#         success.append( r )
#     else:
#         jsuccess=0
#         failed.append( r )
#      
#     rmsd = r.reforiginRmsd
#     if r.reforiginRmsd == None or r.reforiginRmsd == duff:
#         rmsd = -1
#          
#     pdata.append( [ rmsd, r.ensembleNativeMaxsub, jsuccess  ]  )
  
    
# ensembleData = {}
# for r in allResults:
#     if r.pdbCode not in ensembleData.keys():
#         ensembleData[ r.pdbCode ] = {}
#     
#         pfile = os.path.join( TMdir, r.pdbCode, "ROSETTA_MR_0/resultsd.pkl")
#         f = open( pfile )
#         ad = cPickle.load( f  )
#         f.close()
#         
#         ensembleData[ r.pdbCode ][ 'clustersize'] = ad['spicker_results'][0].cluster_size
#         ensembleData[ r.pdbCode ][ 'ensemblesize'] = len( ad['ensemble_results'][0] )
#




cpath = os.path.join( rundir, 'plot.csv' )
csvfile =  open( cpath, 'wb')
csvwriter = csv.writer(csvfile, delimiter=',',
                        quotechar='"', quoting=csv.QUOTE_MINIMAL)

for p in pdata:
    csvwriter.writerow( p )
    
csvfile.close()

sys.exit(0)
# sfail=0
# for s in success:
#     if s.reforiginRmsd == duff:
#         sfail += 1
#  
# ffail=0
# for f in failed:
#     if f.reforiginRmsd == duff:
#         ffail += 1
#  
# print sfail
# print ffail





# print len(success)
# print len(failed)
# 
# clusters = []
# ensembles = []
# for k in sorted( ensembleData.keys() ):
#     c = ensembleData[ k ][ 'clustersize']
#     e = ensembleData[ k ][ 'ensemblesize']
#     clusters.append( c )
#     ensembles.append( e )
#     print "{0}\nClustersize: {1}\nEnsemblesize: {2}\n\n".format( k, c, e )
# 
# 
# print sorted(clusters)
# print sorted(ensembles)
# 
# a=0
# for e in ensembles:
#     a+=e
#     
# print float(a)/len(ensembles)