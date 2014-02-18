import os
import sys
sys.path.append("/Users/jmht/Documents/AMPLE/ample-dev1/scripts")
sys.path.append("/Users/jmht/Documents/AMPLE/ample-dev1/python")

import ample_util
from analyse_run import AmpleResult

ocsv = sys.argv[1] 
ncsv = ample_util.filename_append( filename=ocsv, astr="obj" )

a = AmpleResult()
with open( ocsv ) as o, open( ncsv, 'w' ) as n:
    for i, line in enumerate( o ):
        if i == 0:
            fields = line.strip().split(",")
            nf = []
            for f in fields:
                i = a.orderedTitles.index( f )
                nf.append( a.orderedAttrs[ i ] )
            line = ",".join( nf ) + "\n"
        n.write( line )	
