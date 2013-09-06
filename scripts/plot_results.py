'''
Created on 4 Sep 2013

@author: jmht

r1% = 42.1568627451
r2% = 31.3725490196
r3% = 26.4705882353

polya% = 41.1764705882
reliable% = 27.4509803922
allatom% = 31.3725490196
'''

import sys
sys.path.append("/opt/ample-dev1/scripts")

import cPickle
import csv
import os
from analyse_run import AmpleResult
import math


TMdir = "/media/data/shared/TM"
rundir="/Users/jmht/Dropbox/MRes/Project/Thesis/Analysis"
rundir="/home/jmht/Dropbox/MRes/Project/Thesis/Analysis"

f = "/Users/jmht/Dropbox/MRes/Project/Thesis/Analysis/ar_results.pkl"
f = "/home/jmht/Documents/test/ar_results.pkl"
fh = open(f)
allResults = cPickle.load( fh )
fh.close()


sspred = {}
solvent = {}
sspredC = {}
sspredE = {}
sspredH = {}
success = {}

for r in allResults:
    
    if not success.has_key( r.pdbCode ):
        success[ r.pdbCode ] = 0
        
    if not solvent.has_key( r.pdbCode ):
        solvent[ r.pdbCode ] = r.solventContent

    if r.shelxeCC > 30 and r.shelxeAvgChainLength > 10:
        success[ r.pdbCode ] = 1
     
    if not sspred.has_key( r.pdbCode ):
        dC = math.fabs( r.ss_pred['percentC'] - r.ss_dssp['percentC'] )
        dE = math.fabs( r.ss_pred['percentE'] - r.ss_dssp['percentE'] )
        dH = math.fabs( r.ss_pred['percentH'] - r.ss_dssp['percentH'] )
        sspred[r.pdbCode ] = 100 - (dC + dE + dH) /3
        sspredC[r.pdbCode ] = 100 - dC
        sspredE[r.pdbCode ] = 100 - dE
        sspredH[r.pdbCode ] = 100 - dH

print sspred
print success

pdata = [ ['pdb', 'sspred','sspredC', 'sspredE', 'sspredH', 'solvent', 'success' ] ]
for k in sorted(sspred.keys()):
    pdata.append( [ k, sspred[k], sspredC[k], sspredE[k], sspredH[k], solvent[k], success[k] ] ) 


# pdata = [ ['rfree','rfact','dRF', 'CC','success'] ]
# for r in allResults:
#     success=0
#     if r.shelxeCC > 30 and r.shelxeAvgChainLength > 10:
#         success=1
#     
#     if r.rfree and r.rfact:
#         dRF = float(r.rfree)-float(r.rfact)
#       
#         pdata.append( [ r.rfree, r.rfact, dRF, r.shelxeCC, success  ]  )


# r1=0
# r2=0
# r3=0
# polya=0
# reliable=0
# allatom=0
# nensembles=[]
# nsuccess=0
# for r in allResults:
#     if r.shelxeCC > 30 and r.shelxeAvgChainLength > 10:
#         nsuccess += 1
#         if r.ensembleRadiusThreshold == 1:
#              r1+=1
#         elif r.ensembleRadiusThreshold == 2:
#              r2+=1
#         elif r.ensembleRadiusThreshold == 3:
#              r3+=1
#         else:
#             raise RuntimeError,"BAD RADIUS"
# 
#         if r.ensembleSideChainTreatment == 'SCWRL_reliable_sidechains':
#              reliable+=1
#         elif r.ensembleSideChainTreatment == 'poly_ala':
#              polya+=1
#         elif r.ensembleSideChainTreatment == 'All_atom':
#              allatom+=1
#         else:
#             raise RuntimeError,"BAD sidechain"
#         
#         nensembles.append( r.ensembleNumModels )
#         
# 
# print "r1% = {0}".format( (float(r1)/nsuccess)*100   )
# print "r2% = {0}".format( (float(r2)/nsuccess)*100   )
# print "r3% = {0}".format( (float(r3)/nsuccess)*100   )
# print r1+r2+r3
# print "polya% = {0}".format( (float(polya)/nsuccess)*100   )
# print "reliable% = {0}".format( (float(reliable)/nsuccess)*100   )
# print "allatom% = {0}".format( (float(allatom)/nsuccess)*100   )
# 
# cpath = os.path.join( rundir, 'plot.csv' )
# csvfile =  open( cpath, 'wb')
# csvwriter = csv.writer(csvfile, delimiter=',',
#                         quotechar='"', quoting=csv.QUOTE_MINIMAL)
#  
# csvwriter.writerow( ["numModels"] )
# for n in nensembles:
#     csvwriter.writerow( [n] )
# csvfile.close()
# sys.exit()



# pdata = [ ['phaserLLG','phaserTFZ','CC','success'] ]
# sn=0
# sp=0
# nsuccess=0
# for r in allResults:
#     success=0
#     if r.shelxeCC > 30 and r.shelxeAvgChainLength > 10:
#         nsuccess+=1
#         success=1
#      
#     if r.mrProgram == 'phaser' and r.phaserLLG != None and r.phaserTFZ != None:
#         llg = r.phaserLLG
#         if success==1 and llg < 0:
#             sn+=1
#         elif success==1 and llg >= 0:
#             sp+=1
#         if llg < -1000:
#             llg = -1000
#         pdata.append( [ llg, r.phaserTFZ, r.shelxeCC, success  ]  )
#         
# print sp
# print sn
# print nsuccess

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