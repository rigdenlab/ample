"""
17.02.2016

@author: jmht
"""

# System
import copy
import logging
import os
import shutil
import sys
import unittest

# Custom
import ample_util
import ensembler
import pdb_edit

# Inherit some variables defined in the ensembler module
SIDE_CHAIN_TREATMENTS = ensembler.SIDE_CHAIN_TREATMENTS

_logger = logging.getLogger(__name__)


def align_mustang(models, mustang_exe=None, work_dir=None):
    if not ample_util.is_exe(mustang_exe):
        raise RuntimeError, "Cannot find mustang executable: {0}".format(mustang_exe)
    
    if not work_dir: work_dir = os.getcwd()
    work_dir = os.path.abspath(work_dir)
    if not os.path.isdir(work_dir): os.mkdir(work_dir)
    os.chdir(work_dir)

    logfile = os.path.join(work_dir, 'mustang.log')
    basename = 'mustang'
    cmd = [mustang_exe, '-F', 'fasta', '-o', basename, '-i' ] + models
    rtn = ample_util.run_command(cmd, logfile=logfile, directory=work_dir)
    if not rtn == 0:
        raise RuntimeError, "Error running mustang. Check logfile: {0}".format(logfile)
    
    alignment_file = os.path.join(work_dir, basename + ".afasta")
    if not os.path.isfile(alignment_file): raise RuntimeError, "Could not find alignment file: {0} after running mustang!".format(alignment_file)
    return alignment_file

def align_gesamt(models, gesamt_exe=None, work_dir=None):
    if not ample_util.is_exe(gesamt_exe):
        raise RuntimeError, "Cannot find gesamt executable: {0}".format(gesamt_exe)
    
    if not work_dir: work_dir = os.getcwd()
    work_dir = os.path.abspath(work_dir)
    if not os.path.isdir(work_dir): os.mkdir(work_dir)
    os.chdir(work_dir)
    
    # Need to map chain name to pdb
    model2chain = {}
    for m in models:
        seqd = pdb_edit.sequence(m)
        if len(seqd) != 1: raise RuntimeError, "Model {0} does not contain a single chain, got: {1}".format(seqd.keys())
        model2chain[m] = seqd.keys()[0]
    
    basename = 'gesamt'
    logfile = os.path.join(work_dir, 'gesamt.log')
    alignment_file = os.path.join(work_dir, basename + ".afasta")
    
    # Build up command-line
    cmd = [gesamt_exe]
    # We iterate through the models to make sure the order stays the same
    for m in models: cmd += [ m, '-s', model2chain[m] ]
    cmd += ['-o', '{0}.pdb'.format(basename), '-a', alignment_file]
    
    rtn = ample_util.run_command(cmd, logfile=logfile, directory=work_dir)
    if not rtn == 0:
        raise RuntimeError, "Error running gesamt. Check logfile: {0}".format(logfile)
    
    if not os.path.isfile(alignment_file): raise RuntimeError, "Gesamt did not generate an alignment file.\nPlease check the logfile: {0}".format(logfile)
    return alignment_file
 

class Ensembler(ensembler.Ensembler):
    """Ensemble creator using on multiple distant homologous structures
    """
    
    def __init__(self):
        # Inherit all functions from Parent Ensembler
        ensembler.Ensembler.__init__(self)
        
        return
    
    def generate_ensembles(self,
                           models,
                           alignment_file=None,
                           ensembles_directory=None,
                           gesamt_exe=None,
                           homolog_aligner=None,
                           mustang_exe=None,
                           nproc=None,
                           percent_truncation=None,
                           side_chain_treatments=SIDE_CHAIN_TREATMENTS,
                           truncation_method=None,
                           work_dir=None):
        
        # Work dir set each time
        if not work_dir: raise RuntimeError, "Need to set work_dir!"
        self.work_dir = work_dir
        
        if not percent_truncation:
            percent_truncation = self.percent_truncation
        if not truncation_method:
            truncation_method = self.truncation_method
        if not ensembles_directory:
            self.ensembles_directory = os.path.join(work_dir, "ensembles")
        else:
            self.ensembles_directory = ensembles_directory
        
        if not len(models):
            raise RuntimeError, "Cannot find any models for ensembling!" 
        if not all([os.path.isfile(m) for m in models]):
            raise RuntimeError, "Problem reading models given to Ensembler: {0}".format(models) 
        
        _logger.info('Ensembling models in directory: {0}'.format(self.work_dir))
    
        # Create final ensembles directory
        if not os.path.isdir(self.ensembles_directory): os.mkdir(self.ensembles_directory)
        
        # standardise all the models
        std_models_dir = os.path.join(work_dir, "std_models")
        os.mkdir(std_models_dir)
        std_models = []
        for m in models:
            std_model = ample_util.filename_append(m, 'std', std_models_dir)
            pdb_edit.standardise(pdbin=m, pdbout=std_model, del_hetatm=True)
            std_models.append(std_model)
        
        # Get a structural alignment between the different models
        if not alignment_file:
            if homolog_aligner == 'mustang':
                _logger.info("Generating alignment file with mustang_exe: {0}".format(mustang_exe))
                alignment_file = align_mustang(std_models, mustang_exe=mustang_exe, work_dir=self.work_dir)
            elif homolog_aligner == 'gesamt':
                _logger.info("Generating alignment file with gesamt_exe: {0}".format(gesamt_exe))
                alignment_file = align_gesamt(std_models, gesamt_exe=gesamt_exe, work_dir=self.work_dir)
            else:
                raise RuntimeError, "Unknown homolog_aligner: {0}".format(homolog_aligner)
            _logger.info("Generated alignment file: {0}".format(alignment_file))
        else:
            _logger.info("Using alignment file: {0}".format(alignment_file))
            
        
        truncate_dir = os.path.join(self.work_dir,"homolog_truncate")
        if not os.path.isdir(truncate_dir): os.mkdir(truncate_dir)
            
        # Now truncate and create ensembles - as standard ample, but with no subclustering
        self.ensembles = []
        self.ensembles_data = []
        for truncated_models, truncated_models_data, truncated_model_dir in zip(*self.truncate_models(std_models,
                                                                                                      truncation_method=truncation_method,
                                                                                                      truncation_pruning=None,
                                                                                                      percent_truncation=percent_truncation,
                                                                                                      homologs=True,
                                                                                                      alignment_file=alignment_file,
                                                                                                      work_dir=truncate_dir)):
            tlevel = truncated_models_data['truncation_level']
            ensemble_dir = os.path.join(truncated_model_dir, "ensemble_{0}".format(tlevel))
            os.mkdir(ensemble_dir)
            os.chdir(ensemble_dir)
             
            # Need to create an alignment file for theseus
            basename = "e{0}".format(tlevel)
            pre_ensemble = self.superpose_models(truncated_models, basename=basename, work_dir=ensemble_dir, homologs=True)
            if not pre_ensemble:
                _logger.critical("Skipping ensemble {0} due to error with Theseus".format(basename))
                continue
            pre_ensemble_data = copy.copy(truncated_models_data)
             
            for ensemble, ensemble_data in zip(*self.edit_side_chains(pre_ensemble,
                                                                      pre_ensemble_data,
                                                                      side_chain_treatments,
                                                                      self.ensembles_directory,
                                                                      homologs=True)):
                self.ensembles.append(ensemble)
                self.ensembles_data.append(ensemble_data)
        
        return self.ensembles
    
class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        paths = cls.thisd.split(os.sep)
        cls.ample_dir = os.sep.join(paths[ :-1 ])
        cls.tests_dir = os.path.join(cls.ample_dir, "tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir, 'testfiles')
        
        cls.theseus_exe = ample_util.find_exe('theseus')
        cls.spicker_exe = ample_util.find_exe('spicker')
        cls.maxcluster_exe = ample_util.find_exe('maxcluster')


        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        formatter = logging.Formatter('%(message)s')
        ch.setFormatter(formatter)
        root.addHandler(ch)

        return
    
# hlfsimko: looks like it became redundant with the Integration Testing framework 
#
#    def testHomologs(self):
#        os.chdir(self.thisd)  # Need as otherwise tests that happen in other directories change os.cwd()
#        ensembler = Ensembler()
#        ensembler.theseus_exe = self.theseus_exe
#        
#        work_dir = os.path.join(self.tests_dir, "homologs_test")
#        if os.path.isdir(work_dir): shutil.rmtree(work_dir)
#        os.mkdir(work_dir)
#
#        pdb_list = [ '1ujbA.pdb', '2a6pA.pdb', '3c7tA.pdb']
#        models = [ os.path.join(self.ample_dir, 'examples', 'homologs', pdb) for pdb in pdb_list ]
#        alignment_file = os.path.join(self.ample_dir, 'examples', 'homologs', 'testthree.afasta')
#        ensembles = ensembler.generate_ensembles(models, alignment_file=alignment_file, work_dir=work_dir)
#        self.assertEqual(len(ensembles), 57)
#        shutil.rmtree(work_dir)
#        return
    
    def testGesamt(self):
        gesamt_exe = "/opt/ccp4-devtools/install/bin/gesamt"
        if not ample_util.is_exe(gesamt_exe): return
        
        pdb_list = [ '1ujb.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir, 'examples', 'homologs', pdb) for pdb in pdb_list ]
        work_dir = os.path.join(self.tests_dir, "gesamt_test")
        alignment_file = align_gesamt(models, gesamt_exe=gesamt_exe, work_dir=work_dir)
        self.assertTrue(os.path.isfile(alignment_file))
        shutil.rmtree(work_dir)
        return
    
    def testMustang(self):
        mustang_exe = "/opt/MUSTANG_v3.2.2/bin/mustang-3.2.1"
        if not ample_util.is_exe(mustang_exe): return
        
        pdb_list = [ '1ujb.pdb', '2a6pA.pdb', '3c7tA.pdb']
        models = [ os.path.join(self.ample_dir, 'examples', 'homologs', pdb) for pdb in pdb_list ]
        work_dir = os.path.join(self.tests_dir, "mustang_test")
        alignment_file = align_mustang(models, mustang_exe=mustang_exe, work_dir=work_dir)
        self.assertTrue(os.path.isfile(alignment_file))
        shutil.rmtree(work_dir)
        return

