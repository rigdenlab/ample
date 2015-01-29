'''
Created on 24 Oct 2014

@author: jmht
'''

# Python imports
import csv
import glob
import logging
import os
import shutil

# Our imports
import ample_util
import csymmatch
import maxcluster
import phenixer
import pdb_edit
import pdb_model
import reforigin
import residue_map
import rio
import rosetta_model

_logger=logging.getLogger()

def analyse(amoptd):
    
    if not os.path.isdir(amoptd['benchmark_dir']):
        raise RuntimeError,"Cannot find benchmark dir: {0}".format(amoptd['benchmark_dir'])
    os.chdir(amoptd['benchmark_dir'])

    analysePdb(amoptd)
    
    analyseModels(amoptd)
    
    _logger.info("Benchmark: generating naitive density map")
    # Generate map so that we can do origin searching
    amoptd['nativeDensityMap']=phenixer.generateMap(amoptd['mtz'],
                                                     amoptd['native_pdb'],
                                                     FP=amoptd['F'],
                                                     SIGFP=amoptd['SIGF'],
                                                     FREE=amoptd['FREE'],
                                                     directory=amoptd['benchmark_dir']
                                                     )
    
    data=[]
    # Only look at first cluster
    cluster=0
    
    # Get the ensembling data
    ensemble_results={} # Maps ensemble name to the data object
    if amoptd.has_key('ensemble_results'):
        ensemble_data = amoptd['ensemble_results'][cluster]
        if not len(ensemble_data):
            _logger.critical("Benchmark cannot find any ensemble data!")
            return

        # Get map of ensemble name -> ensemble result
        for e in ensemble_data:
            if ensemble_results.has_key( e.name ):
                raise RuntimeError, "Duplicate key: {0}".format( e.name )
            ensemble_results[ e.name ] = e
                    
    # Get mrbump_results for cluster
    mrbump_results = amoptd['mrbump_results'][cluster]
    if not len(mrbump_results):
        _logger.critical("Benchmark cannot find any mrbump results!")
        return

    for result in mrbump_results:
        
        d=mkDataDict(amoptd)
 
        # Get the ensemble data and add to the MRBUMP data
        edata=ensemble_results[result.ensembleName]
        d['ensembleNumModels'] =  edata.num_models
        d['ensembleNumResidues'] =  edata.num_residues
        d['ensembleSideChainTreatment'] = edata.side_chain_treatment
        d['ensembleRadiusThreshold'] = edata.radius_threshold
        d['ensembleTruncationThreshold'] =  edata.truncation_threshold
        d['ensemblePercentModel'] = int( ( float( edata.num_residues ) / float( amoptd['fasta_length'] ) ) * 100 )
        d['ensembleNumAtoms'] = edata.num_atoms
        d['ensembleCentroidModel'] = edata.centroid_model
        
        #ar.ensembleNativeRMSD = scoreP.rms( eP.centroidModelName )
        d['ensembleNativeTM'] = amoptd['maxComp'].tm(d['ensembleCentroidModel'])
        
        analyseSolution(amoptd,result,d)
        data.append(d)

    fileName=os.path.join(amoptd['benchmark_dir'],'results.csv' )
    writeCsv(fileName,data)
    amoptd['benchmark_results']=data
    return

def writeCsv(fileName,resultList):
    with open(fileName,'wb') as csvfile:
        csvwriter=csv.DictWriter(csvfile,
                                 fieldnames=sorted(resultList[0].keys()),
                                 delimiter=',',
                                 quotechar='"',
                                 quoting=csv.QUOTE_MINIMAL)

        csvwriter.writeheader()
        csvwriter.writerows(resultList)
    return
    
def mkDataDict(amoptd):

    d={}
    d['nativePdbCode']=amoptd['nativePdbCode']
    d['nativePdbTitle']=amoptd['nativePdbTitle']
    d['nativePdbResolution']=amoptd['nativePdbResolution']
    d['nativePdbSolventContent']=amoptd['nativePdbSolventContent']
    d['nativePdbMatthewsCoefficient']=amoptd['nativePdbMatthewsCoefficient']
    d['nativePdbSpaceGroup']=amoptd['nativePdbSpaceGroup']
    d['nativePdbNumAtoms']=amoptd['nativePdbNumAtoms']
    d['nativePdbNumResidues']=amoptd['nativePdbNumResidues']
    
    return d

def analyseSolution(amoptd,result,d):

    _logger.info("Benchmark: analysing result: {0}".format(result.ensembleName))

    if result.program=="phaser":
        mrPdb = result.phaserPdb
    elif result.program=="molrep":
        mrPdb = result.molrepPdb
    elif result.program=="unknown":
        return
    else:
        assert False,result

    if not mrPdb:
        if result.refmacPdb:
            mrPdb=result.refmacPdb
        if not mrPdb:
            return
    
    # Copy in attributes
    addMrbumpData(result, d)
    
    # debug - copy into work directory as reforigin struggles with long pathnames
    shutil.copy(mrPdb, os.path.join(amoptd['benchmark_dir'], os.path.basename(mrPdb)))
    
    mrPdbInfo=pdb_edit.get_info( mrPdb )
    
    d['numPlacedAtoms']=mrPdbInfo.numAtoms()
    d['numPlacedCA']=mrPdbInfo.numCalpha()

    # Get reforigin info
    rmsder = reforigin.ReforiginRmsd()
    rmsder.getRmsd(  nativePdbInfo=amoptd['nativePdbInfo'],
                     placedPdbInfo=mrPdbInfo,
                     refModelPdbInfo=amoptd['refModelPdbInfo'],
                     cAlphaOnly=True,
                     workdir=amoptd['benchmark_dir']
                     )
    d['reforiginRMSD']=rmsder.rmsd


    # 1. run reforigin to generate map for native with mtz (before this routine is called)
    # 2. run get_cc_mtz_pdb to calculate an origin
    mrOrigin=phenixer.ccmtzOrigin(nativeMap=amoptd['nativeDensityMap'], mrPdb=mrPdb)
    
    # offset.pdb is the mrModel shifted onto the new origin use csymmatch to wrap onto native
    csymmatch.Csymmatch().wrapModelToNative("offset.pdb",
                                            amoptd['native_pdb'],
                                            csymmatchPdb=os.path.join(amoptd['benchmark_dir'],
                                            "phaser_{0}_csymmatch.pdb".format(result.ensembleName))
                                            )

    # Score the origin with all-atom and rio
    rioData=rio.Rio().scoreOrigin(mrOrigin,
                                  mrPdbInfo=mrPdbInfo,
                                  nativePdbInfo=amoptd['nativePdbInfo'],
                                  resSeqMap=amoptd['resSeqMap'],
                                  workdir=amoptd['benchmark_dir']
                                  )

    # Set attributes
    d['aaNumContacts']  = rioData.aaNumContacts
    d['rioNumContacts'] = rioData.rioNumContacts
    d['rioInregister']  = rioData.rioInRegister
    d['rioOoRegister'] = rioData.rioOoRegister
    d['rioBackwards']   = rioData.rioBackwards
    d['rio']            = rioData.rioInRegister + rioData.rioOoRegister
    d['rioNoCat']       = rioData.rioNumContacts - ( rioData.rioInRegister + rioData.rioOoRegister )

#     # Now get the helix
#     helixSequence = contacts.Rio().helixFromContacts( contacts=rioData.contacts,
#                                                            dsspLog=dsspLog )
#     if helixSequence is not None:
#         ampleResult.rioHelixSequence = helixSequence
#         ampleResult.rioLenHelix      = len( helixSequence )
#         hfile = os.path.join( workdir, "{0}.helix".format( ampleResult.ensembleName ) )
#         with open( hfile, 'w' ) as f:
#             f.write( helixSequence+"\n" )

    #
    # This purely for checking and so we have pdbs to view
    # 
    # Wrap shelxe trace onto native using Csymmatch
    if not result.shelxePdb is None and os.path.isfile(result.shelxePdb):
        csymmatch.Csymmatch().wrapModelToNative( result.shelxePdb,
                                                 amoptd['native_pdb'],
                                                 origin=mrOrigin,
                                                 workdir=amoptd['benchmark_dir'])

    # Wrap parse_buccaneer model onto native
    if result.buccaneerPdb:
        # Need to rename Pdb as is just called buccSX_output.pdb
        csymmatchPdb = os.path.join(amoptd['benchmark_dir'], "buccaneer_{0}_csymmatch.pdb".format(result.ensembleName))

        csymmatch.Csymmatch().wrapModelToNative( result.buccaneerPdb,
                                                 amoptd['native_pdb'],
                                                 origin=mrOrigin,
                                                 csymmatchPdb=csymmatchPdb,
                                                 workdir=amoptd['benchmark_dir'])
        
    # Wrap parse_buccaneer model onto native
    if result.arpWarpPdb:
        # Need to rename Pdb as is just called buccSX_output.pdb
        csymmatchPdb = os.path.join(amoptd['benchmark_dir'], "arpwarp_{0}_csymmatch.pdb".format(result.ensembleName))

        csymmatch.Csymmatch().wrapModelToNative( result.arpWarpPdb,
                                                 amoptd['native_pdb'],
                                                 origin=mrOrigin,
                                                 csymmatchPdb=csymmatchPdb,
                                                 workdir=amoptd['benchmark_dir'])

    return

def addMrbumpData(result,d):
    """Copy data from mrbump result into dict"""
    
    attrs = [
             # 'jobDir',
             # 'mrDir',
             # 'name',
             # 'pdb',
             # 'ensembleName',
             'program',
             'solution',
             'rfact',
             'rfree',
             'buccRfact',
             'buccRfree',
             'arpWarpRfact',
             'arpWarpRfree',
             'shelxeCC',
             'shelxeACL',
             'buccFinalRfact',
             'buccFinalRfree',
             'arpWarpFinalRfact',
             'arpWarpFinalRfree',
             'phaserLLG',
             'phaserTFZ',
             'phaserTime',
             'phaserKilled',
             'molrepScore',
             'molrepTime'
             ]

    for a in attrs:
        if hasattr(result,a):
            d[a]=getattr(result, a)
        else:
            d[a]=None
    
    return


def analysePdb(amoptd):
    
    nativePdb=amoptd['native_pdb']
    nativePdbInfo = pdb_edit.get_info( nativePdb )
    
    # number atoms/residues
    natoms, nresidues = pdb_edit.num_atoms_and_residues(nativePdb)

    # Get information on the origins for this spaceGroup
    originInfo = pdb_model.OriginInfo( spaceGroupLabel=nativePdbInfo.crystalInfo.spaceGroup )

    # Do this here as a bug in pdbcur can knacker the CRYST1 data
    amoptd['nativePdbCode'] = nativePdbInfo.pdbCode
    amoptd['nativePdbTitle'] = nativePdbInfo.title
    amoptd['nativePdbResolution'] = nativePdbInfo.resolution
    amoptd['nativePdbSolventContent'] = nativePdbInfo.solventContent
    amoptd['nativePdbMatthewsCoefficient'] = nativePdbInfo.matthewsCoefficient
    amoptd['nativePdbSpaceGroup'] = originInfo.spaceGroup()
    amoptd['nativePdbNumAtoms'] = natoms
    amoptd['nativePdbNumResidues'] = nresidues
    
    # First check if the native has > 1 model and extract the first if so
    if len( nativePdbInfo.models ) > 1:
        _logger.info("nativePdb has > 1 model - using first")
        nativePdb1 = ample_util.filename_append( filename=nativePdb, astr="model1", directory=amoptd['benchmark_dir'] )
        pdb_edit.extract_model( nativePdb, nativePdb1, modelID=nativePdbInfo.models[0].serial )
        nativePdb = nativePdb1
        
    # Standardise the PDB to rename any non-standard AA, remove solvent etc
    nativePdbStd = ample_util.filename_append( filename=nativePdb, astr="std", directory=amoptd['benchmark_dir'] )
    pdb_edit.standardise( nativePdb, nativePdbStd )
    nativePdb = nativePdbStd
    
    # Get the new Info about the native
    nativePdbInfo = pdb_edit.get_info( nativePdb )
    
    # For maxcluster comparsion of shelxe model we need a single chain from the native so we get this here
    if len( nativePdbInfo.models[0].chains ) > 1:
        chainID = nativePdbInfo.models[0].chains[0]
        nativeChain1  = ample_util.filename_append( filename=nativePdbInfo.pdb,
                                                       astr="chain1".format( chainID ), 
                                                       directory=amoptd['benchmark_dir'])
        pdb_edit.to_single_chain( nativePdbInfo.pdb, nativeChain1 )
    else:
        nativeChain1 = nativePdbInfo.pdb
    
    # Additional data
    amoptd['nativePdbNumChains'] = len( nativePdbInfo.models[0].chains )
    amoptd['nativePdbInfo']=nativePdbInfo
    amoptd['nativePdbStd']=nativePdbStd
    amoptd['nativePdb1Chain']=nativeChain1
    amoptd['nativePdbOriginInfo']=originInfo
    
    return

def analyseModels(amoptd):
    
    # Get hold of a full model so we can do the mapping of residues
    refModelPdb = glob.glob(os.path.join(amoptd['models_dir'], "*.pdb"))[0]
    
    nativePdbInfo=amoptd['nativePdbInfo']
    
    resSeqMap = residue_map.residueSequenceMap()
    refModelPdbInfo = pdb_edit.get_info(refModelPdb)
    resSeqMap.fromInfo( refInfo=refModelPdbInfo,
                        refChainID=refModelPdbInfo.models[0].chains[0], # Only 1 chain in model
                        targetInfo=nativePdbInfo,
                        targetChainID=nativePdbInfo.models[0].chains[0]
                      )
    amoptd['resSeqMap']=resSeqMap
    amoptd['refModelPdbInfo']=refModelPdbInfo
    
    # Get the scores for the models - we use both the rosetta and maxcluster methods as maxcluster
    # requires a separate run to generate total RMSD
    #if False:
    _logger.info("Analysing RMSD scores for Rosetta models")
    try:
        amoptd['rosettaSP'] = rosetta_model.RosettaScoreParser(amoptd['models_dir'])
    except RuntimeError,e:
        print e
    amoptd['maxComp'] = maxcluster.Maxcluster(amoptd['maxcluster_exe'])
    _logger.info("Analysing Rosetta models with Maxcluster")
    amoptd['maxComp'].compareDirectory( nativePdbInfo=nativePdbInfo,
                              resSeqMap=resSeqMap,
                              modelsDirectory=amoptd['models_dir'],
                              workdir=amoptd['benchmark_dir'])
    
    return


def analyseSS(amoptd):
    # Secondary Structure assignments
    psipred_file = os.path.join( dataDir, "{0}.psipred_ss2".format(pdbCode)  )
    psipredP = PsipredParser( psipred_file )
    dsspLog = os.path.join( dataDir, "{0}.dssp".format( pdbCode ) )
    dsspP = dssp.DsspParser( dsspLog )
    return
