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
import buccaneer

ensembleDir = "/media/data/shared/coiled-coils/ensemble/ensemble_buccaneer/"
#runDir = "/home/jmht/Documents/test/CC/run5"
runDir = os.getcwd()
os.chdir( runDir )

pfile = os.path.join( runDir, "ar_results.pkl" )
with open( pfile ) as f:
    ensembleResults = cPickle.load( f )

# Hack to add extra attributes
a = AmpleResult()
bp = buccaneer.BuccaneerLogParser()
for r in ensembleResults:
    print "processing ",r.pdbCode, r.ensembleName
    r.orderedAttrs = a.orderedAttrs
    r.orderedTitles = a.orderedTitles
    
    blog = os.path.join( ensembleDir,
                         r.pdbCode,
                         "{0}/build/buccaneer.log".format( r.ensembleName ) )

    print "processing log ",blog
    
    if os.path.isfile( blog ):
        bp.parse( blog )
        r.buccFinalRfree = bp.finalRfree
        r.buccFinalRfact = bp.finalRfact
    else:
        r.buccFinalRfree = "nolog"
        r.buccFinalRfact = "nolog"
        print "MISSING LOGFILE ",r.pdbCode,r.ensembleName, blog


pfile = ample_util.filename_append( pfile, astr="bucc")
f = open( pfile, 'w' )
ampleDict = cPickle.dump( ensembleResults, f  )

cpath = os.path.join( runDir, 'results_bucc.csv' )
csvfile =  open( cpath, 'wb')
csvwriter = csv.writer(csvfile, delimiter=',',
                        quotechar='"', quoting=csv.QUOTE_MINIMAL)

header=False
for r in ensembleResults:
    if not header:
        #csvwriter.writerow( r.titlesAsList() )
        csvwriter.writerow( r.valueAttrAsList() )
        header=True
    csvwriter.writerow( r.valuesAsList() )
    
csvfile.close()


# sys.path.append( "/Users/jmht/Documents/AMPLE/ample-dev1/scripts" )
# import pdb_edit
# pe = pdb_edit.PDBEdit()
# i = pe.get_info("/Users/jmht/Documents/AMPLE/data/coiled-coils/pdbs/1G1J.pdb")
# print i.solventContent
