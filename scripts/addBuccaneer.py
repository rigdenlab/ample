"""

"""


import sys
sys.path.append( "/opt/ample-dev1/python" )
sys.path.append( "/opt/ample-dev1/scripts" )

import cPickle
import csv
import os

import ample_util
from analyse_run import AmpleResult

import parse_buccaneer

runDir = "/media/data/shared/coiled-coils/ensemble/ensemble_redo_failures1"
workDir = os.getcwd()

pfile = os.path.join( workDir, "ar_results.pkl" )
with open( pfile ) as f:
    allResults = cPickle.load( f )

# Hack to add extra attributes
#a = AmpleResult()
p = parse_buccaneer.BuccaneerLogParser()
for r in allResults:
    
    print "processing ",r.pdbCode, r.ensembleName
    #r.orderedAttrs = a.orderedAttrs
    #r.orderedTitles = a.orderedTitles
    buccaneerLog = os.path.join( runDir,
                                 r.pdbCode,
                                 "ROSETTA_MR_0",
                                 "MRBUMP",
                                 "cluster_1",
                                 "search_{0}_mrbump".format(r.ensembleName),
                                 "data",
                                 "loc0_ALL_{0}".format(r.ensembleName),
                                 "unmod",
                                 "mr",
                                 "{0}".format(r.mrProgram),
                                 "build",
                                 "shelxe",
                                 "rebuild",
                                 "build",
                                 "buccaneer.log" )
    
    if os.path.isfile(buccaneerLog):
        p.parse(buccaneerLog)
        r.buccFinalRfact = p.finalRfact
        r.buccFinalRfree = p.finalRfree
    
    
pfile = ample_util.filename_append( pfile, astr="bucc")
f = open( pfile, 'w' )
ampleDict = cPickle.dump( allResults, f  )

cpath = os.path.join( workDir, 'results_bucc.csv' )
csvfile =  open( cpath, 'wb')
csvwriter = csv.writer(csvfile, delimiter=',',
                        quotechar='"', quoting=csv.QUOTE_MINIMAL)

header=False
for r in allResults:
    if not header:
        #csvwriter.writerow( r.titlesAsList() )
        csvwriter.writerow( r.valueAttrAsList() )
        header=True
    csvwriter.writerow( r.valuesAsList() )
    
csvfile.close()


