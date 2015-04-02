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
                            'ASU' : 0,
                            'benchmark_mode' : False,
                            'blast_dir' : None,
                            'buccaneer_cycles' : 5,
                            'ccp4_jobid' : None,
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
                            'import_models' : False,
                            'import_ensembles' : False,
                            'improve_template' : None,
                            'LGA' : None,
                            'make_frags' : False,
                            'make_models' : True,
                            'maxcluster_exe' : None,
                            'max_array_jobs' : None,
                            'max_ensemble_models' : 30,
                            'missing_domain' : False,
                            'models' : None,
                            'molrep_only' : False,
                            'mr_keys' : None,
                            'mr_sequence' : None,
                            'mtz' : None,
                            'name' : 'ampl',
                            'native_pdb' : None,
                            'nmodels' : 1000,
                            'NMR_model_in' : None,
                            'NMR_process' : None,
                            'NMR_remodel' : False,
                            'NMR_remodel_fasta' : None,
                            'nproc' : 1,
                            'nr' : None,
                            'num_clusters' : 1,
                            'output_pdb' : 'ample_output.pdb',
                            'percent' : 5,
                            'phaser_only' : False,
                            'phaser_kill' : 0,
                            'phenix_exe' : None,
                            'psipred_ss2' : None,
                            'purge' : False,
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
                            'split_mr' : False,
                            'submit_array' : True,
                            'submit_cluster' : False,
                            'submit_qtype' : None,
                            'submit_queue' : None,
                            'theseus_exe' : None,
                            'top_model_only' : False,
                            'transmembrane' : False,
                            'transmembrane_octopusfile' : None,
                            'transmembrane_lipofile' : None,
                            'transmembrane_spanfile' : None,
                            'truncation_method' : 'percent',
                            'truncation_pruning' : 'none',
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
                           'molrep_only' : False,
                           'phaser_only' : True,
                           'shelx_cycles' : 5,
                           'use_arpwarp' : False,
                           'use_buccaneer' : False,
                           'phaser_kill' : 15
                        }

        # Test use scrwl
        self.devel_mode = {
                           'early_terminate': False,
                           'benchmark_mode': True,
                           'phaser_only': True,
                           'use_shelxe' : True,
                           'shelxe_rebuild' : True,
                           'shelxe_rebuild_arpwarp' : True,
                           'shelxe_rebuild_buccaneer' : True,
                           'use_scwrl' : False,
                           'use_arpwarp' : False,
                           'use_buccaneer' : False,
                           # This kills phaser after 6 hours
                           'phaser_kill' : 360,
                           #'mr_keys' : [ [ 'PKEY', 'KILL','TIME','360'  ] ],
                        }
        
        self.webserver_uri = {
                               'max_array_jobs' : 10,
                               'purge': True,
                               'shelxe_rebuild_buccaneer': True,
                               'submit_cluster' : True,
                               'submit_qtype' : "SGE",
                               'submit_queue' : "all.q",
                               }

        # We have a debug mode as the logger isn't activated when we run
        self.debug = False

        return
    
    def populate( self, parser_args ):
        """
        Fill ourselves with the options from the parser
        """

        tmpv = None
        for k, v in vars(parser_args).iteritems():
            #print "{} | {}".format( k, v )
            if isinstance(v,list):
                # All values are in a list
                tmpv  = v[0]
            else:
                tmpv = v

            # Bit of a hack for true/false
            if isinstance( tmpv, str ):
                if tmpv.lower() == "true":
                    tmpv = True
                elif tmpv.lower() == "false":
                    tmpv = False

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
                'constraints_file',
                'domain_all_chains_pdb',
                'ensembles_dir',
                'fasta',
                'fast_protein_cluster_exe',
                'frags_3mers',
                'frags_9mers',
                'import_cluster',
                'maxcluster_exe',
                'models',
                'mr_sequence',
                'mtz',
                'native_pdb',
                'NMR_model_in',
                'NMR_remodel_fasta',
                'psipred_ss2',
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
            if k in paths and isinstance( v, str ):
                self.d[ k ] = os.path.abspath( v )

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
        pstr = ""
        pstr +='Params Used in this Run\n\n'

        keys1 = ['fasta','work_dir','mtz','name']
        pstr += '---input---\n'
        for k in keys1:
            pstr += "{0}: {1}\n".format(k, self.d[k])

        keys2 = ['make_frags','rosetta_fragments_exe','frags_3mers','frags_9mers']
        pstr+= '\n---fragments---\n'
        for k in keys2:
            pstr += "{0}: {1}\n".format(k, self.d[k])

        keys3 = ['make_models','rosetta_AbinitioRelax','rosetta_db']
        pstr+= '\n---modelling---\n'
        for k in keys3:
            pstr += "{0}: {1}\n".format(k, self.d[k])

        if self.d['use_scwrl']:
            pstr+= '\n---3rd party---\nSCWRL {0}\n'.format( self.d['scwrl_exe'] )

        keys4 = ['missing_domain','domain_all_chains_pdb']
        if keys4[0]:
            pstr+= '\n---Missing Domain---\n'
            for k in keys4:
                pstr += "{0}: {1}\n".format(k, self.d[k])

        # This only used for printing
        INSERT_DOMAIN = False
        if self.d['domain_termini_distance'] > 0:
            INSERT_DOMAIN = True
        pstr += '\nIs an Insert Domain {0} termini distance {1}\n'.format( INSERT_DOMAIN, self.d['domain_termini_distance'] )

        # Now print out everything else
        pstr += "\n---Other parameters---\n"

        done_keys = keys1 + keys2 + keys3 + keys4 + [ 'use_scwrl', 'domain_termini_distance'  ]
        for k, v in sorted(self.d.items()):
            if k not in done_keys:
                pstr += "{0} : {1}\n".format( k, v )

        return pstr
