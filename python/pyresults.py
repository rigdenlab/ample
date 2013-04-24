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

root="/gpfs/home/HCEA041/djr01/jxt15-djr01/TM"
root="/home/Shared/TM"

dirs = [ "1GU8", "2BHW", "2BL2", "2EVU", "2O9G", "2UUI", "2WIE", "2X2V", "2XOV", "3GD8", "3HAP", "3LBW", "3LDC", "3OUF", "3PCV", "3RLB", "3TX3", "3U2F", "4DVE" ]
dirs = [ "1GU8", "2BHW", "2BL2", "2EVU", "2O9G", "2UUI", "2WIE", "2X2V", "2XOV", "3GD8", "3HAP", "3LBW", "3LDC", "3OUF", "3PCV", "3RLB", "3U2F", "4DVE" ]
dirs = [ "1GU8", "2BHW", "2BL2", "2EVU", "2O9G", "2UUI", "2WIE", "2X2V", "2XOV", "3GD8", "3HAP",  "3LDC", "3OUF", "3PCV", "3RLB", "3U2F", "4DVE" ]
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
mpkl = "ENSEMBLE_1_MOLREP_0/ensemble_results.pkl"
ppkl = "ENSEMBLE_1_PHASER0/ensemble_results.pkl"

logging.basicConfig()
logging.getLogger().setLevel(logging.CRITICAL)

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
        if r.name == name:
            # if only...
            #return r.num_models, r.num_residues
            # Horrible, but I got the nmodels wrong...            
            pdb = r.pdb
            s = pdb.split(os.sep)
            path =  os.path.join( root, os.sep.join( s[7:] ) )
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
    
    # Take first cluster
    eresults = rdict['ensemble_results'][0]
    
    #for e in eresults:
    #    print e.num_models
    
    # MRBUMP molrep
    if os.path.exists( os.path.join( dpath, mpkl ) ):
        pklf = open( os.path.join( dpath, mpkl ) )
        rdict = cPickle.load( pklf  )
        pklf.close()
        mrb_molrep = rdict['mrbump_results'][0]
        # We changed things while doing this run
        for i,res in enumerate(mrb_molrep):
            mrb_molrep[i].name = res.name[9:-6]
    else:
        # unfinished run
        r = mrbump_results.ResultsSummary( os.path.join(dpath,"ENSEMBLE_1_MOLREP_0/MRBUMP/cluster_1"), cluster=False )
        #print r.summariseResults()
        r.extractResults()
        mrb_molrep = r.results
    
    # MRBUMP phaser
    if os.path.exists( os.path.join( dpath, ppkl ) ):
        pklf = open( os.path.join( dpath, ppkl ) )
        rdict = cPickle.load( pklf  )
        pklf.close()
        mrb_phaser = rdict['mrbump_results'][0]
        # We changed things while doing this run
        for i,res in enumerate(mrb_phaser):
            mrb_phaser[i].name = res.name[9:-6]
    else:
        # unfinished run
        r = mrbump_results.ResultsSummary( os.path.join(dpath,"ENSEMBLE_1_PHASER0/MRBUMP/cluster_1"), cluster=False )
        r.extractResults()
        mrb_phaser = r.results
    
    
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
    
    print "RESULTS FOR MOLREP"
    for result in mrb_molrep:
        #name = result.name[9:-6]
        nmodels, nresidues = qensemble( result.name, eresults)
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
    
    resultsTable = []
    resultsTable.append( ("name", "MR_program", "Solution", "final_Rfact", "final_Rfree", "SHELXE_CC", "#Models", "#Residues") )
    
    print "RESULTS FOR PHASER"
    for result in mrb_phaser:
        #name = result.name[9:-6]
        nmodels, nresidues = qensemble( result.name, eresults)
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