#!/usr/bin/env python


import csv
import cPickle
import logging
import os
import sys
#sys.path.append("/gpfs/home/HCEA041/djr01/jxt15-djr01/ample-dev1/python")
sys.path.append("/opt/ample-dev1/python")

import mrbump_results
import run_spicker
import ensemble
import printTable


def len_ensemble( epdb ):

    f = open( epdb, 'r' )

    nmodels = 0
    residues = 0
    # Find where first model starts
    line = f.readline().strip()
    while line:
        line = f.readline().strip()
        
        if line.startswith("MODEL "):
            nmodels+=1
            # Count the residues & atoms
            atoms=0
            residues=0
            last=None
        
            aline = f.readline().strip()
            while aline:
        
                if not aline.startswith("ATOM"):
                    if not f.readline().strip().startswith("ENDMDL"):
                        raise RuntimeError,"Error looking for end of atoms"
                    break
        
                atoms+=1
                # See if the residue is different from thelast
                fields = aline.split()
                nresidue = fields[5]
                if not last or nresidue != last:
                    last = nresidue
                    residues+=1
        
                aline = f.readline().strip()
            #end while

    f.close()

    return nmodels, residues


#print len_ensemble( "/home/Shared/TM/1GU8/ENSEMBLES_0/ensembles_1/poly_ala_trunc_38.476482_rad_2.pdb" )


#sys.exit()

root="/home/Shared/TM"
root="/gpfs/home/HCEA041/djr01/jxt15-djr01/TM"

dirs = [ "1GU8", "2BHW", "2BL2", "2EVU", "2O9G", "2UUI", "2WIE", "2X2V", "2XOV", "3GD8", "3HAP", "3LBW", "3LDC", "3OUF", "3PCV", "3RLB", "3TX3", "3U2F", "4DVE" ]
dirs = [ "1GU8", "2BHW", "2BL2", "2EVU", "2O9G", "2UUI", "2WIE", "2X2V", "2XOV", "3GD8", "3HAP", "3LBW", "3LDC", "3OUF", "3PCV", "3RLB", "3U2F", "4DVE" ]
# Need to think about 3LBW
#dirs = [ "1GU8", "2BHW", "2BL2", "2EVU", "2O9G", "2UUI", "2WIE", "2X2V", "2XOV", "3GD8", "3HAP", "3LDC", "3OUF", "3PCV", "3RLB", "3U2F", "4DVE" ]
#dirs = [ "2BL2" ]
#dirs = [ "1GU8", "2BHW", "2BL2", "2EVU", "2O9G", "2UUI", "2WIE", "2X2V", "2XOV", "3GD8", "3HAP",  "3LDC", "3OUF", "3PCV", "3RLB", "3U2F", "4DVE" ]
#dirs = [ "1GU8", "2BHW", "2EVU", "2O9G", "2UUI", "2XOV", "3GD8", "3HAP",  "3RLB"]
#dirs = [ "2WIE"]


# Get the  data on each TM
# fields are: title, resolution, length
TMdict = {}
with open( os.path.join(root, "misc/TM_Data.csv"), "r" ) as inputf:
    
    reader = csv.reader(inputf, delimiter=',', quotechar='"')
    
    for i,fields in enumerate(reader):
        if i ==0:
            continue
        TMdict[ fields[0] ] = fields[1:]

epkl = "ENSEMBLES_0/ensemble_results.pkl"
bpkl = "ENSEMBLE_2_0/ample_results.pkl"

logging.basicConfig()
logging.getLogger().setLevel(logging.CRITICAL)
#logging.getLogger().setLevel(logging.DEBUG)

# Want
# sizes of all spicker clusters (centroid model?)
# truncation level, number of models & number residues in models
#  for rads 1, 2 & 3 -
#     number of models
#
# mrbump - top result - number residues, number of models in search model
# for all successful results, min & max number of models, size


def qensemble(name, ensemble_results):
    """Get nmodels & nresidues for the given model form the ensemble results"""
    for r in ensemble_results:
        #print "Checking {0} against {1}".format(name, r.name)
        #print "Checking {0} against {1}".format(name[9:-6], r.name)
        if r.name == name[9:-6]:
            # if only...
            #return r.num_models, r.num_residues
            # Horrible, but I got the nmodels wrong...            
            pdb = r.pdb
            s = pdb.split(os.sep)
            path =  os.path.join( root, os.sep.join( s[7:] ) )
            #print "checking path ",path
            return len_ensemble( path )


table = printTable.Table()

for tdir in dirs:
    
    print "\n\nResults for: {0}".format(tdir)
    print "Title: {0}\nResolution: {1}\nLength: {2}".format( TMdict[tdir][0],TMdict[tdir][1],TMdict[tdir][2] )
    dpath = os.path.join( root, tdir )
    
    pklf = open( os.path.join( dpath, epkl ) )
    rdict = cPickle.load( pklf  )
    pklf.close()
    
    # spicker
    spicker_results = rdict['spicker_results']
    clusters=""
    for sr in  spicker_results:
        clusters+="{0}:".format(sr.cluster_size)
    clusters=clusters[:-1]
    print "Spicker clusters: "+clusters
    
    # Take second cluster
    eresults = rdict['ensemble_results'][1]
    
    # MRBUMP molrep
    #if False and os.path.exists( os.path.join( dpath, bpkl ) ):
    if os.path.exists( os.path.join( dpath, bpkl ) ):
        pklf = open( os.path.join( dpath, bpkl ) )
        rdict = cPickle.load( pklf  )
        pklf.close()
        mrb_results = rdict['mrbump_results'][0]
        # We changed things while doing this run
        for i,res in enumerate(mrb_results):
            #mrb_results[i].name = res.name[9:-6]
            pass
    else:
        # unfinished run
        r = mrbump_results.ResultsSummary( os.path.join(dpath,"ENSEMBLE_2_0/MRBUMP/cluster_1") )
        #print r.summariseResults()
        r.extractResults()
        mrb_results = r.results
    
    
    # Combine and resort
    #mrb_results = mrb_phaser + mrb_molrep
    #sortf = lambda x: float( x.shelxCC )
    #mrb_results.sort(key=sortf)
    #mrb_results.reverse()
    
#    resultsTable = []
#    resultsTable.append( ("name", "MR_program", "Solution", "final_Rfact", "final_Rfree", "SHELXE_CC", "#Models", "#Residues") )
#        
#    for result in mrb_results:
#        #name = result.name[9:-6]
#        nmodels, nresidues = qensemble( result.name, eresults)
#        rl = [ result.name,
#              result.program,
#              result.solution,
#              result.rfact,
#              result.rfree,
#              result.shelxCC,
#              nmodels,
#              nresidues
#              ]
#        resultsTable.append( rl )


    resultsTable = []
    resultsTable.append( ("name", "MR_program", "Solution", "final_Rfact", "final_Rfree", "SHELXE_CC", "#Models", "#Residues") )
    
    for result in mrb_results:
        #name = result.name[9:-6]
        (nmodels, nresidues) = qensemble( result.name, eresults)
        rl = [ result.name,
              result.program,
              result.solution,
              result.rfact,
              result.rfree,
              result.shelxCC,
              nmodels,
              nresidues
              ]
        resultsTable.append( rl )

    
    print table.pprint_table( resultsTable )
