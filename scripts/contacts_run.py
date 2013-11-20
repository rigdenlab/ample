#!/usr/bin/env python


import cPickle
import os
import shutil
import sys

#sys.path.append("/Users/jmht/Documents/AMPLE/ample-dev1/python")
sys.path.append("/opt/ample-dev1/python")
sys.path.append("/opt/ample-dev1/scripts")


import ample_util
import analyse_run
import contacts
import dssp
import mrbump_results
import pdb_edit
import reforigin
import residue_map

if __name__ == "__main__":
    
    CLUSTERNUM=0
    rundir = "/home/jmht/Documents/test/ncont/new"
    rootDir = "/media/data/shared/coiled-coils"
    os.chdir( rundir )
    
#    pfile = os.path.join( CCdir,"results.pkl" )
#    f = open( pfile )
#    resultsDict = cPickle.load( f  )
#    f.close()
    
    allResults = []
    
    #for pdbcode in [ "1GU8", "2BHW", "2BL2", "2EVU", "2O9G", "2UUI", "2WIE", "2X2V", "2XOV", "3GD8", "3HAP", "3LBW", "3LDC", "3OUF", "3PCV", "3RLB", "3U2F", "4DVE" ]:
    #for pdbcode in sorted( resultsDict.keys() ):
    for pdbcode in [ l.strip() for l in open( os.path.join( rootDir, "dirs.list") ) if not l.startswith("#") ]:
    #for pdbcode in [ "1DEB" ]:
        
        workdir = os.path.join( rundir, pdbcode )
        if not os.path.isdir( workdir ):
            os.mkdir( workdir )
        os.chdir( workdir )
            
        print "\nResults for ",pdbcode
        
        # Directory where all the data for this run live
        dataDir = os.path.join( rootDir, pdbcode )
        
        # Get the path to the original pickle file
        pfile = os.path.join( dataDir, "ROSETTA_MR_0/resultsd.pkl")
        f = open( pfile )
        ampleDict = cPickle.load( f  )
        f.close()    
    
        # First process all stuff that's the same for each structure
        
        # Get path to native Extract all the nativeInfo from it
        nativePdb = os.path.join( dataDir, "{0}.pdb".format( pdbcode ) )
        pdbedit = pdb_edit.PDBEdit()
        nativeInfo = pdbedit.get_info( nativePdb )
        
        # First check if the native has > 1 model and extract the first if so
        if len( nativeInfo.models ) > 1:
            print "nativePdb has > 1 model - using first"
            nativePdb1 = ample_util.filename_append( filename=nativePdb, astr="model1", directory=workdir )
            pdbedit.extract_model( nativePdb, nativePdb1, modelID=nativeInfo.models[0].serial )
            nativePdb = nativePdb1
            
        # Standardise the PDB to rename any non-standard AA, remove solvent etc
        nativePdbStd = ample_util.filename_append( filename=nativePdb, astr="std", directory=workdir )
        pdbedit.standardise( nativePdb, nativePdbStd )
        nativePdb = nativePdbStd
        
        # Secondary Structure assignments
        dssp_file = os.path.join( dataDir, "{0}.dssp".format( pdbcode  )  )
        #dssp_file = os.path.join( dataDir, "{0}.dssp".format( pdbcode.lower()  )  )
        dsspP = dssp.DsspParser( dssp_file )
        
        # See if the residue numbers match, and if not create a residue map
        
        # Get hold of a full model so we can do the mapping of residues
        refModelPdb = os.path.join( dataDir, "models/S_00000001.pdb".format( pdbcode ) )
        resSeqMap = residue_map.residueSequenceMap()
        
        modelInfo = pdbedit.get_info( refModelPdb )
        
        # Need to update nativeInfo as it might have changed
        nativeInfo = pdbedit.get_info( nativePdb )
        
        #try:
        resSeqMap.fromInfo( refInfo=modelInfo,
                                refChainID=modelInfo.models[0].chains[0],
                                targetInfo=nativeInfo,
                                targetChainID=nativeInfo.models[0].chains[0]
                                )
        
        #except Exception, e:
        #    print "GOT EXCEPTION {0}".format( e )
        
        #continue
    
        # Loop over each result
        r = mrbump_results.ResultsSummary( os.path.join( dataDir, "ROSETTA_MR_0/MRBUMP/cluster_1") )
        r.extractResults()
        for mrbumpResult in r.results:
        #for mrbumpResult in resultsDict[ pdbcode ]:
            
            #print "processing result ",mrbumpResult
            
            # yuck...
            if mrbumpResult.solution == 'unfinished':
                s = mrbumpResult.name.split("_")
                #ensembleName = "_".join( s[2:-2] )
                ensembleName = "_".join( s[2:-1] )
            else:
                # MRBUMP Results have loc0_ALL_ prepended and  _UNMOD appended
                ensembleName = mrbumpResult.name[9:-6]
                
            #if ensembleName != "poly_ala_trunc_0.466257_rad_2":
            #    continue
            
            # Extract information on the models and ensembles
            eresults = ampleDict['ensemble_results']
            got=False
            for e in ampleDict[ 'ensemble_results' ][ CLUSTERNUM ]:
                if e.name == ensembleName:
                    got=True
                    break 
        
            if not got:
                raise RuntimeError,"Failed to get ensemble results"
        
            # No results so move on
            if mrbumpResult.solution == 'unfinished':
                continue
            
            # Need to remove last component as we recorded the refmac directory
            resultDir = os.sep.join( mrbumpResult.resultDir.split(os.sep)[:-1] )

            if mrbumpResult.program != "phaser":
                continue
                
            print "Checking ensemble ",ensembleName
            
            phaserPdb = os.path.join( resultDir,"refine","{0}_loc0_ALL_{1}_UNMOD.1.pdb".format(mrbumpResult.program, ensembleName) )
            if not os.path.isfile( phaserPdb ):
                continue
                
            placedPdb = phaserPdb
                
            # Get the reforigin RMSD of the phaser placed model as refined with refmac
            #refinedPdb = os.path.join( resultDir, "refine", "refmac_{0}_loc0_ALL_{1}_UNMOD.pdb".format( mrbumpResult.program, ensembleName ) )
            
            # debug - copy into work directory as reforigin struggles with long pathnames
            shutil.copy(placedPdb, os.path.join( workdir, os.path.basename( placedPdb ) ) )
    
            #
            # SHELXE PROCESSING
            #
            # Now read the shelxe log to see how we did
            shelxeLog = os.path.join( resultDir, "build/shelxe/shelxe_run.log" )
            shelxePdb = os.path.join( resultDir, "build/shelxe", "shelxe_{0}_loc0_ALL_{1}_UNMOD.pdb".format( mrbumpResult.program, ensembleName ) )
            if not os.path.isfile( shelxePdb):
                print "NO SHELXE PDB"
                continue
            shutil.copy( shelxePdb , os.path.join(workdir, os.path.basename( shelxePdb ) ) )
 
            if os.path.isfile( shelxeLog ):
                shelxeP = analyse_run.ShelxeLogParser( shelxeLog )
                shelxeCC = shelxeP.CC
                shelxeAvgChainLength = shelxeP.avgChainLength
                print "SHELXE ",shelxeCC, shelxeAvgChainLength
 
            # Now calculate contacts
            ccalc = contacts.Contacts()
            ccalc.getContacts( nativePdb=nativePdb, placedPdb=placedPdb, resSeqMap=resSeqMap, nativeInfo=nativeInfo, shelxePdb=shelxePdb, workdir=workdir )
             
            if ccalc.best:
                good = ccalc.best.inregister + ccalc.best.ooregister
                bad = ccalc.best.numContacts - good
                print "GOT BEST ", good, bad, ccalc.best.backwards, ccalc.best.inregister, ccalc.best.ooregister
                 
                # Show origin stats
                oc = sorted(ccalc.originCompare.items(), key=lambda x: x[1], reverse=True )
                print "originCompare: ", oc
                duff=False
                if len(oc) > 1:
                    if oc[0][1] == oc[1][1]:
                        if len(oc) > 2:
                            if oc[2][1] >= oc[1][1]*.5:
                                duff=True
                    else:
                        if oc[1][1] >= oc[0][1]*.5:
                            duff=True
                    if duff:
                        print "OTHER ORIGINMATCHES ARE > 50%"
                 
                hfile = os.path.join( workdir, "{0}.helix".format( ensembleName ) )
                ccalc.writeHelixFile(filename=hfile, dsspP=dsspP )
            else:
                if shelxeCC >= 25 and shelxeAvgChainLength >= 10:
                    print "SUCCESS BUT NO BEST!"
            
            if True:
                placedInfo = pdbedit.get_info( placedPdb )
                rmsder = reforigin.ReforiginRmsd( nativePdb=nativePdb,
                                                  nativePdbInfo=nativeInfo,
                                                  placedPdb=placedPdb,
                                                  placedPdbInfo=placedInfo,
                                                  refModelPdb=refModelPdb,
                                                  cAlphaOnly=True )
                print "GOT REFORIGIN ",rmsder.rmsd
     
    # End loop over results
