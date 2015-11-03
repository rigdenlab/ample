'''
Created on 24 Oct 2014

@author: jmht
'''

# Python imports
import copy
import cPickle
import glob
import logging
import os
import shutil
import unittest

# Our imports
import ample_exit
import ample_util
import csymmatch
import maxcluster
import mrbump_results
import pdb_edit
import pdb_model
import reforigin
import residue_map
import rio

_logger=logging.getLogger()

_oldroot=None
_newroot=None

def fixpath(path):
    # fix for analysing on a different machine
    if _oldroot and _newroot:
        return os.path.join(_newroot,path[len(_oldroot)+1:])
    else:
        return path

def analyse(amoptd, newroot=None):
    if newroot:
        #if newroot.endswith("/"): newroot=newroot[:-1]
        assert os.path.isdir(newroot)
        global _oldroot,_newroot
        _newroot=newroot
        _oldroot=amoptd['work_dir']
        # hack
        amoptd['maxcluster_exe']="/opt/maxcluster/maxcluster"
        amoptd['native_pdb']=os.path.join("/media/data/shared/testset/data",os.path.basename(amoptd['native_pdb']))
        amoptd['models_dir']=os.path.join("/media/data/shared/testset/models",amoptd['native_pdb_code'],"models")
    
    if not os.path.isdir(fixpath(amoptd['benchmark_dir'])):
        raise RuntimeError,"Cannot find benchmark dir: {0}".format(amoptd['benchmark_dir'])
    os.chdir(fixpath(amoptd['benchmark_dir']))

    analysePdb(amoptd)
    if not (amoptd['ideal_helices'] or amoptd['homologs']):
        analyseModels(amoptd)
    
#     _logger.info("Benchmark: generating naitive density map")
#     # Generate map so that we can do origin searching
#     amoptd['native_density_map']=phenixer.generateMap(amoptd['mtz'],
#                                                      amoptd['native_pdb'],
#                                                      FP=amoptd['F'],
#                                                      SIGFP=amoptd['SIGF'],
#                                                      FREE=amoptd['FREE'],
#                                                      directory=amoptd['benchmark_dir'])
    # Get the ensembling data
    if not len(amoptd['ensembles_data']):
        _logger.critical("Benchmark cannot find any ensemble data!")
        return

    # Get dict of ensemble name -> ensemble result
    ensemble_results = { e['name'] : e for e in amoptd['ensembles_data'] }
                    
    # Get mrbump_results for cluster
    mrbump_results = amoptd['mrbump_results']
    if not len(mrbump_results):
        _logger.critical("Benchmark cannot find any mrbump results!")
        return
    
    data=[]
    for result in mrbump_results:
        
        # use mrbump dict as basis for result object
        d = copy.copy(result)
        
        # Add in the data from the ensemble
        d.update(ensemble_results[d['ensemble_name']])
        assert d['ensemble_name'] == d['name'],d
        
        # General stuff
        d['ample_version'] = amoptd['ample_version']
        
        # Add in stuff we've cleaned from the pdb
        d['native_pdb_code'] = amoptd['native_pdb_code']
        d['native_pdb_title'] = amoptd['native_pdb_title']
        d['native_pdb_resolution'] = amoptd['native_pdb_resolution']
        d['native_pdb_solvent_content'] = amoptd['native_pdb_solvent_content']
        d['native_pdb_space_group'] = amoptd['native_pdb_space_group']
        d['native_pdb_num_chains'] = amoptd['native_pdb_num_chains']
        d['native_pdb_num_atoms'] = amoptd['native_pdb_num_atoms']
        d['native_pdb_num_residues'] = amoptd['native_pdb_num_residues']
 
        # Get the ensemble data and add to the MRBUMP data
        d['ensemble_percent_model'] = int((float(d['num_residues']) / float(amoptd['fasta_length'])) * 100)
        #ar.ensembleNativeRMSD = scoreP.rms( eP.centroidModelName )
        
        # Need to get the subcluster_centroid_model and then get the path to the original model
        if 'subcluster_centroid_model' in d:
            n=os.path.basename(d['subcluster_centroid_model'])
            cm=None
            for pdb in glob.glob(os.path.join(amoptd['models_dir'],"*.pdb")):
                if os.path.basename(pdb)==n:
                    cm=pdb
                    break
            if not cm:
                raise RuntimeError,"Cannot find model for subcluster_centroid_model {0}".format(d['subcluster_centroid_model'])
            #cm=d['cluster_centroid']
            d['ensemble_native_TM'] = amoptd['maxComp'].tm(cm)
            d['ensemble_native_RMSD'] = amoptd['maxComp'].rmsd(cm)
        analyseSolution(amoptd,d)
        data.append(d)

    fileName=os.path.join(fixpath(amoptd['benchmark_dir']),'results.csv' )
    writeCsv(fileName,data)
    amoptd['benchmark_results']=data
    return

def writeCsv(fileName,resultList):
    
    # List of all the keys we want to write out in order
    keylist= [
                'ample_version',
                
                # Native info
                'native_pdb_code',
                'native_pdb_title',
                'native_pdb_resolution',
                'native_pdb_solvent_content',
                'native_pdb_space_group',
                'native_pdb_num_atoms',
                'native_pdb_num_residues',
                'native_pdb_num_chains',
                
                # Get the ensemble data and add to the MRBUMP data
                'ensemble_name',
                'ensemble_percent_model',
                'ensemble_native_TM',
                'ensemble_native_RMSD',
                
                # cluster info
                'cluster_method',
                'num_clusters',
                'cluster_num',
                'cluster_centroid',
                'cluster_num_models',
                
                # truncation info
                'truncation_level',
                'percent_truncation',
                'truncation_method',
                'truncation_pruning',
                'truncation_variance',
                'num_residues',
                'pruned_residues',
                
                # subclustering info
                'subcluster_num_models',
                'subcluster_radius_threshold',
                'subcluster_centroid_model',
                
                # ensemble info
                #'name',
                'side_chain_treatment',
                'ensemble_num_atoms',
                
                # MR result info
                #'name',
                'MR_program',
                'Solution_Type',
                
                'PHASER_LLG',
                'PHASER_TFZ',
                'PHASER_RFZ',
                'PHASER_time',
                'PHASER_killed',
                'PHASER_version',
                'PHASER_errors',
                
                'MOLREP_score',
                'MOLREP_time',
                'MOLREP_version',
                
                'REFMAC_Rfact',
                'REFMAC_Rfree',
                'REFMAC_version',
                
                'BUCC_final_Rfact',
                'BUCC_final_Rfree',
                'BUCC_version',
                
                'ARP_final_Rfact',
                'ARP_final_Rfree',
                'ARP_version',
                
                'SHELXE_CC',
                'SHELXE_ACL',
                'SHELXE_MCL',
                'SHELXE_NC',
                'SHELXE_wMPE',
                'SHELXE_os',
                'SHELXE_time',
                'SHELXE_version',
                
                'SXRBUCC_version',
                'SXRBUCC_final_Rfact',
                'SXRBUCC_final_Rfree',
                
                'SXRARP_version',
                'SXRARP_final_Rfact',
                'SXRARP_final_Rfree',
                
                'num_placed_atoms',
                'reforigin_RMSD',
                
                'AA_num_contacts',
                'RIO_num_contacts',
                'RIO_in_register',
                'RIO_oo_register',
                'RIO_backwards',
                'RIO',
                'RIO_no_cat',
                'RIO_norm'
    
    ]
    
    #for d in resultList:
    #    for k in sorted(d.keys()):
    #        print "GOT ",k,d[k]
    
    with open(fileName,'wb') as csvfile:
        csvfile.write(",".join(keylist)+"\n")
        for d in resultList:
            #csvfile.write(",".join([d[k] for k in keylist]+"\n"))
            values=[]
            for k in keylist:
                if k in d and d[k] is not None:
                    values.append(str(d[k]).replace(",","^"))
                else:
                    values.append("N/A")
            csvfile.write(",".join(values)+"\n")
        csvfile.write("\n")

# Doesnt' seem to work       
#         csvwriter=csv.DictWriter(csvfile,
#                                  fieldnames=keylist,
#                                  extrasaction='ignore',
#                                  restval='N/A',
#                                  delimiter=',',
#                                  quotechar='"',
#                                  quoting=csv.QUOTE_MINIMAL)
#         csvwriter.writeheader()
#         csvwriter.writerows(resultList)
    return
    
def analyseSolution(amoptd,d):

    _logger.info("Benchmark: analysing result: {0}".format(d['ensemble_name']))

    mrPdb=None
    if d['MR_program']=="PHASER":
        mrPdb = d['PHASER_pdbout']
    elif d['MR_program']=="MOLREP":
        mrPdb = d['MOLREP_pdbout']
    elif d['MR_program']=="unknown":
        return

    if mrPdb is None or not os.path.isfile(mrPdb):
        #_logger.critical("Cannot find mrPdb {0} for solution {1}".format(mrPdb,d))
        return

    # debug - copy into work directory as reforigin struggles with long pathnames
    shutil.copy(mrPdb, os.path.join(fixpath(amoptd['benchmark_dir']), os.path.basename(mrPdb)))
    
    mrPdbInfo=pdb_edit.get_info( mrPdb )
    
    d['num_placed_atoms']=mrPdbInfo.numAtoms()
    d['num_placed_CA']=mrPdbInfo.numCalpha()

    # Get reforigin info
    rmsder = reforigin.ReforiginRmsd()
    try:
        rmsder.getRmsd(nativePdbInfo=amoptd['native_pdb_info'],
                       placedPdbInfo=mrPdbInfo,
                       refModelPdbInfo=amoptd['ref_model_pdb_info'],
                       cAlphaOnly=True,
                       workdir=fixpath(amoptd['benchmark_dir']))
        d['reforigin_RMSD']=rmsder.rmsd
    except Exception,e:
        _logger.critical("Error calculating RMSD: {0}".format(e))
        d['reforigin_RMSD']=999

    # Find the MR origin wrt to the native
    #mrOrigin=phenixer.ccmtzOrigin(nativeMap=amoptd['native_density_map'], mrPdb=mrPdb)
    #mrOrigin=shelxe.shelxe_origin(amoptd['shelxe_exe'],amoptd['native_pdb'],amoptd['mtz'],mrPdb=mrPdb)
    if not d['SHELXE_os']:
        _logger.critical("mrPdb {0} has no SHELXE_os origin shift.\n{1}".format(mrPdb,d))
        return
    mrOrigin=[c*-1 for c in d['SHELXE_os']]
    
    # Move pdb onto new origin
    originPdb=ample_util.filename_append(mrPdb, astr='offset',directory=fixpath(amoptd['benchmark_dir']))
    pdb_edit.translate(mrPdb, originPdb, mrOrigin)
    
    # offset.pdb is the mrModel shifted onto the new origin use csymmatch to wrap onto native
    csymmatch.Csymmatch().wrapModelToNative(originPdb,
                                            amoptd['native_pdb'],
                                            csymmatchPdb=os.path.join(fixpath(amoptd['benchmark_dir']),
                                            "phaser_{0}_csymmatch.pdb".format(d['ensemble_name']))
                                            )

    if not amoptd['ideal_helices']:
        # Score the origin with all-atom and rio
        rioData=rio.Rio().scoreOrigin(mrOrigin,
                                      mrPdbInfo=mrPdbInfo,
                                      nativePdbInfo=amoptd['native_pdb_info'],
                                      resSeqMap=amoptd['res_seq_map'],
                                      workdir=fixpath(amoptd['benchmark_dir'])
                                      )
    
        # Set attributes
        d['AA_num_contacts']  = rioData.aaNumContacts
        d['RIO_num_contacts'] = rioData.rioNumContacts
        d['RIO_in_register']  = rioData.rioInRegister
        d['RIO_oo_register']  = rioData.rioOoRegister
        d['RIO_backwards']    = rioData.rioBackwards
        d['RIO']              = rioData.rioInRegister + rioData.rioOoRegister
        d['RIO_no_cat']       = rioData.rioNumContacts - ( rioData.rioInRegister + rioData.rioOoRegister )
        d['RIO_norm']         = float(d['RIO']) / float(d['native_pdb_num_residues'])
    else:
        d['AA_num_contacts']  = None
        d['RIO_num_contacts'] = None
        d['RIO_in_register']  = None
        d['RIO_oo_register']  = None
        d['RIO_backwards']    = None
        d['RIO']              = None
        d['RIO_no_cat']       = None
        d['RIO_norm']         = None

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
    if not d['SHELXE_pdbout'] is None and os.path.isfile(fixpath(d['SHELXE_pdbout'])):
        csymmatch.Csymmatch().wrapModelToNative( fixpath(d['SHELXE_pdbout']),
                                                 amoptd['native_pdb'],
                                                 origin=mrOrigin,
                                                 workdir=fixpath(amoptd['benchmark_dir']))

    # Wrap parse_buccaneer model onto native
    if d['SXRBUCC_pdbout'] and os.path.isfile(fixpath(d['SXRBUCC_pdbout'])):
        # Need to rename Pdb as is just called buccSX_output.pdb
        csymmatchPdb = os.path.join(fixpath(amoptd['benchmark_dir']), "buccaneer_{0}_csymmatch.pdb".format(d['ensemble_name']))

        csymmatch.Csymmatch().wrapModelToNative( fixpath(d['SXRBUCC_pdbout']),
                                                 amoptd['native_pdb'],
                                                 origin=mrOrigin,
                                                 csymmatchPdb=csymmatchPdb,
                                                 workdir=fixpath(amoptd['benchmark_dir']))
        
    # Wrap parse_buccaneer model onto native
    if d['SXRARP_pdbout'] and os.path.isfile(fixpath(d['SXRARP_pdbout'])):
        # Need to rename Pdb as is just called buccSX_output.pdb
        csymmatchPdb = os.path.join(fixpath(amoptd['benchmark_dir']), "arpwarp_{0}_csymmatch.pdb".format(d['ensemble_name']))

        csymmatch.Csymmatch().wrapModelToNative( fixpath(d['SXRARP_pdbout']),
                                                 amoptd['native_pdb'],
                                                 origin=mrOrigin,
                                                 csymmatchPdb=csymmatchPdb,
                                                 workdir=fixpath(amoptd['benchmark_dir']))

    return

def analysePdb(amoptd):
    
    nativePdb=amoptd['native_pdb']
    nativePdbInfo = pdb_edit.get_info( nativePdb )
    
    # number atoms/residues
    natoms, nresidues = pdb_edit.num_atoms_and_residues(nativePdb)

    # Get information on the origins for this spaceGroup
    originInfo = pdb_model.OriginInfo( spaceGroupLabel=nativePdbInfo.crystalInfo.spaceGroup )

    # Do this here as a bug in pdbcur can knacker the CRYST1 data
    amoptd['native_pdb_code'] = nativePdbInfo.pdbCode
    amoptd['native_pdb_title'] = nativePdbInfo.title
    amoptd['native_pdb_resolution'] = nativePdbInfo.resolution
    amoptd['native_pdb_solvent_content'] = nativePdbInfo.solventContent
    amoptd['native_pdb_matthews_coefficient'] = nativePdbInfo.matthewsCoefficient
    amoptd['native_pdb_space_group'] = originInfo.spaceGroup()
    amoptd['native_pdb_num_atoms'] = natoms
    amoptd['native_pdb_num_residues'] = nresidues
    
    # First check if the native has > 1 model and extract the first if so
    if len( nativePdbInfo.models ) > 1:
        _logger.info("nativePdb has > 1 model - using first")
        nativePdb1 = ample_util.filename_append( filename=nativePdb, astr="model1", directory=fixpath(amoptd['benchmark_dir']))
        pdb_edit.extract_model( nativePdb, nativePdb1, modelID=nativePdbInfo.models[0].serial )
        nativePdb = nativePdb1
        
    # Standardise the PDB to rename any non-standard AA, remove solvent etc
    nativePdbStd = ample_util.filename_append( filename=nativePdb, astr="std", directory=fixpath(amoptd['benchmark_dir']))
    pdb_edit.standardise(nativePdb, nativePdbStd, del_hetatm=True)
    nativePdb = nativePdbStd
    
    # Get the new Info about the native
    nativePdbInfo = pdb_edit.get_info( nativePdb )
    
    # For maxcluster comparsion of shelxe model we need a single chain from the native so we get this here
    if len( nativePdbInfo.models[0].chains ) > 1:
        chainID = nativePdbInfo.models[0].chains[0]
        nativeChain1  = ample_util.filename_append( filename=nativePdbInfo.pdb,
                                                       astr="chain1".format( chainID ), 
                                                       directory=fixpath(amoptd['benchmark_dir']))
        pdb_edit.to_single_chain( nativePdbInfo.pdb, nativeChain1 )
    else:
        nativeChain1 = nativePdbInfo.pdb
    
    # Additional data
    amoptd['native_pdb_num_chains'] = len( nativePdbInfo.models[0].chains )
    amoptd['native_pdb_info'] = nativePdbInfo
    amoptd['native_pdb_std'] = nativePdbStd
    amoptd['native_pdb_1chain'] = nativeChain1
    amoptd['native_pdb_origin_info'] = originInfo
    
    return

def analyseModels(amoptd):
    
    # Get hold of a full model so we can do the mapping of residues
    refModelPdb = glob.glob(os.path.join(amoptd['models_dir'], "*.pdb"))[0]
    
    nativePdbInfo=amoptd['native_pdb_info']
    
    resSeqMap = residue_map.residueSequenceMap()
    refModelPdbInfo = pdb_edit.get_info(refModelPdb)
    resSeqMap.fromInfo( refInfo=refModelPdbInfo,
                        refChainID=refModelPdbInfo.models[0].chains[0], # Only 1 chain in model
                        targetInfo=nativePdbInfo,
                        targetChainID=nativePdbInfo.models[0].chains[0]
                      )
    amoptd['res_seq_map']=resSeqMap
    amoptd['ref_model_pdb_info']=refModelPdbInfo
    
    # Get the scores for the models - we use both the rosetta and maxcluster methods as maxcluster
    # requires a separate run to generate total RMSD
    #if False:
#     _logger.info("Analysing RMSD scores for Rosetta models")
#     try:
#         amoptd['rosettaSP'] = rosetta_model.RosettaScoreParser(amoptd['models_dir'])
#     except RuntimeError,e:
#         print e
    amoptd['maxComp'] = maxcluster.Maxcluster(amoptd['maxcluster_exe'])
    _logger.info("Analysing Rosetta models with Maxcluster")
    amoptd['maxComp'].compareDirectory( nativePdbInfo=nativePdbInfo,
                              resSeqMap=resSeqMap,
                              modelsDirectory=amoptd['models_dir'],
                              workdir=fixpath(amoptd['benchmark_dir']))
    
    return


def analyseSS(amoptd):
    # Secondary Structure assignments
    psipred_file = os.path.join( dataDir, "{0}.psipred_ss2".format(pdbCode)  )
    psipredP = PsipredParser( psipred_file )
    dsspLog = os.path.join( dataDir, "{0}.dssp".format( pdbCode ) )
    dsspP = dssp.DsspParser( dsspLog )
    return

def restartPkl(amoptd):
    """Currently a bit of a hack to restart from an existing pkl file.
    
    Assumption is that the pkl file we are given is for an AMPLE directory that
    hasn't been moved.
    We create a new AMPLE_X directory and work in there.
    """
    if not os.path.isfile(amoptd['restart_pkl']):
        msg = 'Cannot find pkl file: {0}'.format(amoptd['restart_pkl'])
        ample_exit.exit(msg)
    
    with open(amoptd['restart_pkl']) as f: amoptd_old = cPickle.load(f)
    
    # Update any variables that have changed - everything else uses the old paths
    amoptd_old['ample_log'] = amoptd['ample_log']
    if os.path.isfile(amoptd['native_pdb']): amoptd_old['native_pdb'] = amoptd['native_pdb']
    amoptd_old['run_dir'] = amoptd['run_dir']
    amoptd_old['results_path'] = os.path.join(amoptd['work_dir'],'resultsd.pkl')
    amoptd_old['work_dir'] = amoptd['work_dir']
    if amoptd['mrbump_dir']:
        if not os.path.isdir(amoptd['mrbump_dir']):
            msg = "Cannot find MRBUMP directory: {0}".format(amoptd['mrbump_dir'])
            ample_exit.exit(msg)
        amoptd_old['mrbump_dir'] = amoptd['mrbump_dir']
    
    # We can now replace the old dictionary with this new one
    amoptd = amoptd_old

    # Process the MRBUMP results and save to dictionary
    res_sum = mrbump_results.ResultsSummary()
    res_sum.extractResults(amoptd['mrbump_dir'])
    amoptd['mrbump_results'] = res_sum.results
    
    # Analse the full results
    amoptd['benchmark_dir'] = os.path.join(amoptd['work_dir'],"benchmark")
    os.mkdir(amoptd['benchmark_dir'])
    analyse(amoptd, newroot=None)
    ample_util.saveAmoptd(amoptd)
    
    ample_exit.exit("restart_pkl exiting")
    return


class Test(unittest.TestCase):
    def testBenchmark(self):
        pklfile="/home/jmht/ample-dev1/examples/toxd-example/ROSETTA_MR_0/resultsd.pkl"
        with open(pklfile) as f: d=cPickle.load(f)
        bd="/home/jmht/ample-dev1/python/foo"
        if not os.path.isdir(bd): os.mkdir(bd)
        d['benchmark_dir']=bd
        analyse(d)
        return
#
# Run unit tests
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(testSuite())

