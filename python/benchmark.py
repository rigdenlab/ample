'''
Created on 24 Oct 2014

@author: jmht
'''

# Python imports
import copy
import cPickle
import csv
import glob
import logging
import os
import shutil
import unittest

# Our imports
import ample_util
import csymmatch
import maxcluster
import pdb_edit
import pdb_model
import reforigin
import residue_map
import rio
import shelxe
#import rosetta_model

_logger=logging.getLogger()

def analyse(amoptd):
    
    if not os.path.isdir(amoptd['benchmark_dir']):
        raise RuntimeError,"Cannot find benchmark dir: {0}".format(amoptd['benchmark_dir'])
    os.chdir(amoptd['benchmark_dir'])

    analysePdb(amoptd)
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
        assert d['ensemble_name']==d['name'],d
        
        # Add in stuff we've cleaned from the pdb
        d['native_pdb_code']=amoptd['native_pdb_code']
        d['native_pdb_title']=amoptd['native_pdb_title']
        d['native_pdb_resolution']=amoptd['native_pdb_resolution']
        d['native_pdb_solvent_content']=amoptd['native_pdb_solvent_content']
        d['native_pdb_space_group']=amoptd['native_pdb_space_group']
        d['native_pdb_num_atoms']=amoptd['native_pdb_num_atoms']
        d['native_pdb_num_residues']=amoptd['native_pdb_num_residues']
 
        # Get the ensemble data and add to the MRBUMP data
        d['ensemble_percent_model'] = int( ( float( d['num_residues'] ) / float( amoptd['fasta_length'] ) ) * 100 )
        #ar.ensembleNativeRMSD = scoreP.rms( eP.centroidModelName )
        d['ensemble_native_TM'] = amoptd['maxComp'].tm(d['cluster_centroid'])
        d['ensemble_native_RMSD'] = amoptd['maxComp'].rmsd(d['cluster_centroid'])
        
        analyseSolution(amoptd,d)
        data.append(d)

    fileName=os.path.join(amoptd['benchmark_dir'],'results.csv' )
    writeCsv(fileName,data)
    amoptd['benchmark_results']=data
    return

def writeCsv(fileName,resultList):
    
    # Hack find all possible keys
    keys=set()
    for r in resultList:
        keys.update(r.keys())
    keys=list(sorted(keys))
    
    with open(fileName,'wb') as csvfile:
        csvwriter=csv.DictWriter(csvfile,
                                 fieldnames=keys,
                                 delimiter=',',
                                 quotechar='"',
                                 quoting=csv.QUOTE_MINIMAL)

        csvwriter.writeheader()
        csvwriter.writerows(resultList)
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
        #for k in sorted(d.keys()):
        #    print k,d[k]
        _logger.critical("Cannot find mrPdb {0} for solution {1}".format(mrPdb,d))
        return

    # debug - copy into work directory as reforigin struggles with long pathnames
    shutil.copy(mrPdb, os.path.join(amoptd['benchmark_dir'], os.path.basename(mrPdb)))
    
    mrPdbInfo=pdb_edit.get_info( mrPdb )
    
    d['num_placed_atoms']=mrPdbInfo.numAtoms()
    d['num_placed_CA']=mrPdbInfo.numCalpha()

    # Get reforigin info
    rmsder = reforigin.ReforiginRmsd()
    rmsder.getRmsd(nativePdbInfo=amoptd['native_pdb_info'],
                   placedPdbInfo=mrPdbInfo,
                   refModelPdbInfo=amoptd['ref_model_pdb_info'],
                   cAlphaOnly=True,
                   workdir=amoptd['benchmark_dir'])
    d['reforigin_RMSD']=rmsder.rmsd

    # Find the MR origin wrt to the native
    #mrOrigin=phenixer.ccmtzOrigin(nativeMap=amoptd['native_density_map'], mrPdb=mrPdb)
    mrOrigin=shelxe.shelxeOrigin(amoptd['shelxe_exe'],amoptd['native_pdb'],amoptd['mtz'],mrPdb=mrPdb)
    
    # Move pdb onto new origin
    originPdb=ample_util.filename_append(mrPdb, astr='offset')
    pdb_edit.translate(mrPdb, originPdb, mrOrigin)
    
    # offset.pdb is the mrModel shifted onto the new origin use csymmatch to wrap onto native
    csymmatch.Csymmatch().wrapModelToNative(originPdb,
                                            amoptd['native_pdb'],
                                            csymmatchPdb=os.path.join(amoptd['benchmark_dir'],
                                            "phaser_{0}_csymmatch.pdb".format(d['ensemble_name']))
                                            )

    # Score the origin with all-atom and rio
    rioData=rio.Rio().scoreOrigin(mrOrigin,
                                  mrPdbInfo=mrPdbInfo,
                                  nativePdbInfo=amoptd['native_pdb_info'],
                                  resSeqMap=amoptd['res_seq_map'],
                                  workdir=amoptd['benchmark_dir']
                                  )

    # Set attributes
    d['AA_num_contacts']  = rioData.aaNumContacts
    d['RIO_num_contacts'] = rioData.rioNumContacts
    d['RIO_in_register']  = rioData.rioInRegister
    d['RIO_oo_register'] = rioData.rioOoRegister
    d['RIO_backwards']   = rioData.rioBackwards
    d['RIO']            = rioData.rioInRegister + rioData.rioOoRegister
    d['RIO_no_cat']       = rioData.rioNumContacts - ( rioData.rioInRegister + rioData.rioOoRegister )

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
    if not d['SHELXE_pdbout'] is None and os.path.isfile(d['SHELXE_pdbout']):
        csymmatch.Csymmatch().wrapModelToNative( d['SHELXE_pdbout'],
                                                 amoptd['native_pdb'],
                                                 origin=mrOrigin,
                                                 workdir=amoptd['benchmark_dir'])

    # Wrap parse_buccaneer model onto native
    if d['SXRBUCC_pdbout']:
        # Need to rename Pdb as is just called buccSX_output.pdb
        csymmatchPdb = os.path.join(amoptd['benchmark_dir'], "buccaneer_{0}_csymmatch.pdb".format(d['ensemble_name']))

        csymmatch.Csymmatch().wrapModelToNative( d['SXRBUCC_pdbout'],
                                                 amoptd['native_pdb'],
                                                 origin=mrOrigin,
                                                 csymmatchPdb=csymmatchPdb,
                                                 workdir=amoptd['benchmark_dir'])
        
    # Wrap parse_buccaneer model onto native
    if d['SXRARP_pdbout']:
        # Need to rename Pdb as is just called buccSX_output.pdb
        csymmatchPdb = os.path.join(amoptd['benchmark_dir'], "arpwarp_{0}_csymmatch.pdb".format(d['ensemble_name']))

        csymmatch.Csymmatch().wrapModelToNative( d['SXRARP_pdbout'],
                                                 amoptd['native_pdb'],
                                                 origin=mrOrigin,
                                                 csymmatchPdb=csymmatchPdb,
                                                 workdir=amoptd['benchmark_dir'])

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
    amoptd['native_pdb_num_chains'] = len( nativePdbInfo.models[0].chains )
    amoptd['native_pdb_info']=nativePdbInfo
    amoptd['native_pdb_std']=nativePdbStd
    amoptd['native_pdb_1chain']=nativeChain1
    amoptd['native_pdb_origin_info']=originInfo
    
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
                              workdir=amoptd['benchmark_dir'])
    
    return


def analyseSS(amoptd):
    # Secondary Structure assignments
    psipred_file = os.path.join( dataDir, "{0}.psipred_ss2".format(pdbCode)  )
    psipredP = PsipredParser( psipred_file )
    dsspLog = os.path.join( dataDir, "{0}.dssp".format( pdbCode ) )
    dsspP = dssp.DsspParser( dsspLog )
    return


class Test(unittest.TestCase):

    def testBenchmark(self):
        pklfile="/opt/ample-dev1.testset/examples/toxd-example/ROSETTA_MR_4/resultsd.pkl"
        with open(pklfile) as f:
            d=cPickle.load(f)
        bd="/opt/ample-dev1.testset/python/foo"
        if not os.path.isdir(bd): os.mkdir(bd)
        d['benchmark_dir']=bd
        analyse(d)
        
        print d

        return

def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testBenchmark'))
    return suite

#
# Run unit tests
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(testSuite())

