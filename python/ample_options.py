'''
Class to hold the options for ample
'''
# python imports
import os

# our imports
import version

class AmpleOptions(object):

    def __init__(self):
        # The dictionary with all the options
        self.d = {}

        # dictionary with the default arguments - if any are paths add to paths in populate
        self.defaults = {
                            'alignment_file' : None,
                            'all_atom' : True,
                            'arpwarp_cycles' : 10,
                            'benchmark_mode' : False,
                            'blast_dir' : None,
                            'buccaneer_cycles' : 5,
                            'ccp4_jobid' : None,
                            'cluster_dir' : None,
                            'cluster_method' : 'spicker',
                            'constraints_file' : None,
                            'debug' : False,
                            'domain_all_chains_pdb' : None,
                            'domain_termini_distance' : 0,
                            'dry_run' : False,
                            'early_terminate' : True,
                            'ensembles_dir' : None,
                            'F' : None,
                            'fasta' : None,
                            'fast_protein_cluster_exe' : None,
                            'frags_3mers' : None,
                            'frags_9mers' : None,
                            'FREE' : None,
                            'gesamt_exe' : None,
                            'homologs' : False,
                            'homolog_aligner' : 'gesamt',
                            'ideal_helices' : False,
                            'import_cluster' : False,
                            'import_models' : False,
                            'import_ensembles' : False,
                            'improve_template' : None,
                            'LGA' : None,
                            'make_frags' : True,
                            'make_models' : True,
                            'maxcluster_exe' : None,
                            'max_ensemble_models' : 30,
                            'missing_domain' : False,
                            'models' : None,
                            'molrep_only' : False,
                            'mrbump_dir' : None,
                            'mr_keys' : None,
                            'mr_sequence' : None,
                            'mtz' : None,
                            'mustang_exe' : None,
                            'name' : 'ampl',
                            'native_pdb' : None,
                            'nmasu' : 0,
                            'nmodels' : 1000,
                            'nmr_model_in' : None,
                            'nmr_process' : None,
                            'nmr_remodel' : False,
                            'nmr_remodel_fasta' : None,
                            'no_gui' : False,
                            'nproc' : 1,
                            'nr' : None,
                            'num_clusters' : 1,
                            'output_pdb' : 'ample_output.pdb',
                            'percent' : 5,
                            'phaser_only' : True,
                            'phaser_kill' : 360, # This kills phaser after 6 hours
                            'phaser_rms' : 0.1,
                            'phenix_exe' : None,
                            'psipred_ss2' : None,
                            'purge' : False,
                            'restart_pkl' : None,
                            'rg_reweight' : None,
                            'ROSETTA_cluster' : None,
                            'rosetta_db' : None,
                            'rosetta_dir' : None,
                            'rosetta_fragments_exe' : None,
                            'rosetta_AbinitioRelax' : None,
                            'rosetta_version' : None,
                            'rcdir' : os.path.join( os.path.expanduser("~"), ".ample" ),
                            'run_dir' : os.getcwd(),
                            'scwrl_exe' : None,
                            'sf_cif' : None,
                            'shelx_cycles' : 15,
                            'shelxe_exe' : None,
                            'shelxe_rebuild' : False,
                            'shelxe_rebuild_arpwarp' : False,
                            'shelxe_rebuild_buccaneer' : False,
                            'SIGF' : None,
                            'spicker_exe' : None,
                            'submit_array' : True,
                            'submit_cluster' : False,
                            'submit_max_array' : None,
                            'submit_qtype' : None,
                            'submit_queue' : None,
                            'theseus_exe' : None,
                            'top_model_only' : False,
                            'transmembrane' : False,
                            'transmembrane2' : False,
                            'transmembrane_octopusfile' : None,
                            'transmembrane_lipofile' : None,
                            'transmembrane_spanfile' : None,
                            'truncation_method' : 'percent',
                            'truncation_pruning' : None,
                            'webserver_uri' : None,
                            'use_arpwarp' : True,
                            'use_buccaneer' : True,
                            'use_homs' : True,
                            'use_scwrl' : False,
                            'use_shelxe' : True,
                         }

        self.quick_mode = {
                           'max_ensemble_models' : 10,
                           'nmodels' : 200,
                           'percent' : 20,
                           'shelx_cycles' : 5,
                           'use_arpwarp' : False,
                           'use_buccaneer' : False,
                           'phaser_kill' : 15
                        }

        # Test use scrwl
        self.devel_mode = {
                           'early_terminate': False,
                           'benchmark_mode': True,
                           'shelxe_rebuild' : True,
                           'shelxe_rebuild_arpwarp' : True,
                           'shelxe_rebuild_buccaneer' : True,
                           'use_arpwarp' : False,
                           'use_buccaneer' : False,
                           #'mr_keys' : [ [ 'PKEY', 'KILL','TIME','360'  ] ],
                        }
        
        self.webserver_uri = {
                               'purge': True,
                               'shelxe_rebuild_buccaneer': True,
                               'submit_cluster' : True,
                               'submit_max_array' : 10,
                               'submit_qtype' : "SGE",
                               'submit_queue' : "all.q",
                               }

        # We have a debug mode as the logger isn't activated when we run
        self.debug = False

        return
    
    def populate(self, parser_args):
        """
        Fill ourselves with the options from the parser
        """

        tmpv = None
        for k, v in vars(parser_args).iteritems():
            #print "{} | {}".format(k, v)
            if isinstance(v,list):
                # All values are in a list
                tmpv  = v[0]
            else:
                tmpv = v

            # Bit of a hack for true/false and ccp4i "none" strings
            if isinstance(tmpv, str):
                if tmpv.lower() == "true":
                    tmpv = True
                elif tmpv.lower() == "false":
                    tmpv = False
                elif tmpv.lower() == "none":
                    tmpv = None

            self.d[k] = tmpv
        # end of loop

#        print "After populate"
#        for k, v in self.d.iteritems():
#            print "{} | {}".format( k, v )

        # Handle any defaults and any preset options
        self.process_options()

        return

    def process_options(self):
        """Check the options and process any preset defaults"""
        
        # Add the version
        self.d['ample_version']=version.__version__

        # First set anything that hasn't been set to its default option
        for k, v in self.defaults.iteritems():
            if k not in self.d:
                self.d[k] = v
            elif  self.d[k] == None:
                #if self.debug:
                #    print "Setting default value: {0} : {1}".format(k,v)
                self.d[k] = v
            else:
                if self.debug and self.d[k] != v:
                    print "Changed default value: {0} : {1}".format(k, self.d[k])

        # Any changes here
        if self.d['submit_qtype']:
            self.d['submit_qtype'] = self.d['submit_qtype'].upper()
        
        if self.d['shelxe_rebuild']:
            self.d['shelxe_rebuild_arpwap']=True
            self.d['shelxe_rebuild_buccaneer']=True

        # Convert all paths to absolute paths
        paths = [
                 'alignment_file',
                'blast_dir',
                'cluster_dir',
                'constraints_file',
                'domain_all_chains_pdb',
                'ensembles_dir',
                'fasta',
                'fast_protein_cluster_exe',
                'frags_3mers',
                'frags_9mers',
                'gesamt_exe',
                'import_cluster',
                'maxcluster_exe',
                'models',
                'mrbump_dir',
                'mr_sequence',
                'mtz',
                'mustang_exe',
                'native_pdb',
                'nmr_model_in',
                'nmr_remodel_fasta',
                'psipred_ss2',
                'restart_pkl',
                'rosetta_db',
                'rosetta_dir',
                'rosetta_fragments_exe',
                'rosetta_AbinitioRelax',
                'scwrl_exe',
                'sf_cif',
                'shelxe_exe',
                'spicker_exe',
                'theseus_exe',
                'transmembrane_octopusfile',
                'transmembrane_lipofile',
                'transmembrane_spanfile'
            ]
        for k, v in self.d.iteritems():
            if k in paths and isinstance(v, str):
                self.d[k] = os.path.abspath(v)

        # Check if using any preset options
        if self.d['devel_mode']: self.preset_options('devel_mode')
        if self.d['quick_mode']:self.preset_options('quick_mode')
        if self.d['webserver_uri']:self.preset_options('webserver_uri')
        return
    
    def preset_options(self,mode):
        assert hasattr(self,mode),"Unknown mode: {0}".format(mode)
        for k, v in getattr(self, mode).iteritems():
            # Set any that haven't been set
            if self.d[k] == None:
                self.d[k] = v
            else:
                # Already set - only overwrite if it's set to a default value, otherwise we
                # let the user go with what they've chosen but warn
                if self.d[k] != self.defaults[k] and self.d[k] != v  :
                    print "WARNING! Overriding {0} setting: {1} : {2} with user setting {3}".format(mode, k, v, self.d[k])
                else:
                    # We overwrite the default with our value
                    if self.debug:
                        print "Overriding default setting: {0} : {1} with {2} setting {3}".format(k, mode, self.defaults[k], v)
                    self.d[k] = v
        return
       
    def prettify_parameters(self):
        """
        Return the parameters nicely formated as a list of strings suitable for writing out to a file
        """
        pstr ='Parameters Used in this Run\n\n'
        for k, v in sorted(self.d.items()):
            pstr += "{0} : {1}\n".format(k, v)
        return pstr
