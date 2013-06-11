'''
Class to hold the options for ample
'''
# python imports
import os

# Our imports
import printTable

class AmpleOptions(object):
    
    def __init__(self):
        

        
        # The dictionary with all the options
        self.d = {}
        
        # dictionary with the default arguments
        self.defaults = {
                            'alignment_file' : None,
                            'all_atom' : True,
                            'arpwarp_cycles' : 10,
                            'ASU' : 0,
                            'blast_dir' : None,
                            'buccaneer_cycles' : 5,
                            'CC' : None,
                            'debug' : None,
                            'domain_all_chains_pdb' : None,
                            'domain_termini_distance' : 0,
                            'early_terminate' : True,
                            'ensembler' : False,
                            'ensembles_dir' : None,
                            'F' : None,
                            'fasta' : None,
                            'frags_3mers' : None,
                            'frags_9mers' : None,
                            'FREE' : None,
                            'import_cluster' : False,
                            'improve_template' : None,
                            'LGA' : None,
                            'make_frags' : True,
                            'make_models' : True,
                            'maxcluster_exe' : None,
                            'max_ensemble_models' : 30,
                            'missing_domain' : False,
                            'models_dir' : None,
                            'molrep_only' : False,
                            'mr_keys' : None,
                            'mtz' : None,
                            'name' : None,
                            'nmodels' : 1000,
                            'NMR_model_in' : None,
                            'NMR_process' : None,
                            'NMR_remodel_fasta' : None,
                            'NMR_Truncate_only' : None,
                            'nproc' : 1,
                            'nr' : None,
                            'num_clusters' : 1,
                            'old_shelx' : False,
                            'percent' : 5,
                            'phaser_only' : False,
                            'phenix_exe' : None,
                            #'ROSETTA' : None,
                            'ROSETTA_cluster' : None,
                            'rosetta_db' : None,
                            'rosetta_dir' : None,
                            'rosetta_fragments_exe' : None,
                            'rosetta_path' : None,
                            'run_dir' : os.getcwd(),
                            'scwrl_exe' : None,
                            'sf_cif' : None,
                            'shelx_cycles' : 15,
                            'shelxe_exe' : None,
                            'SIGF' : None,
                            'spicker_exe' : None,
                            'split_mr' : False,
                            'submit_cluster' : False,
                            'submit_qtype' : None,
                            'theseus_exe' : None,
                            'top_model_only' : False,
                            'transmembrane' : False,
                            'transmembrane_lipofile' : None,
                            'transmembrane_spanfile' : None,
                            'use_arpwarp' : True,
                            'use_buccaneer' : True,
                            'use_homs' : True,
                            'use_scwrl' : False,
                            'use_shelxe' : False,

                         }
        
        self.quick_mode = {
                           'max_ensemble_models' : 10,
                           'nmodels' : 200,
                           'percent' : 20,
                           'phaser_only' : True,
                           'shelx_cycles' : 5,
                           'use_arpwarp' : True,
                           'use_buccaneer' : True,
                           # Needs to be a list of lists as there can be multiple mr_keys
                           # This kills phaser after 15 min - add when the CCP4 version of phaser supports it
                           # 'mr_keys' : [ [ 'PKEY', 'KILL','TIME','15'  ] ],
                        }
    
        # Test use scrwl
        self.devel_mode = {
                           'early_terminate': False,
                           'use_shelxe' : True,
                           'use_scwrl' : False,
                           'use_arpwarp' : False,
                           'use_buccaneer' : False,
                        }
    
        # We have a debug mode as the logger isn't activated when we run
        self.debug = False
        
        
    def final_summary(self, cluster=None):
        """Return a string summarising the results of the run.
        
        Args:
        cluster -- the number of the cluster to summarise (COUNTING FROM 0)
                   otherwise all clusters are summarised
        """
        
        if not cluster:
            # Get number of clusters from the length of the mrbump results lists
            clusters = [ i for i in range( len( self.d['mrbump_results'] ) ) ]
        else:
            if cluster >= len( self.d['mrbump_results'] ):
                raise RuntimeError, "Cluster number is not in results list"
            clusters = [ cluster ]
            
        # String to hold the results summary
        summary = ""
        for cluster in clusters:
            
            summary+="\n\nResults for cluster: {0}\n\n".format( cluster+1 )
            
            # Table for results with header
            results_table = []

            ensemble_results = None
            if self.d.has_key('ensemble_results'):
                ensemble_results = self.d['ensemble_results'][ cluster ]
            
                name2e = {}
                # Get map of name -> ensemble result
                for i, e in enumerate( ensemble_results ):
                    if name2e.has_key( e.name ):
                        raise RuntimeError, "Duplicate key: {0}".format( e.name )
                    name2e[ e.name ] = ensemble_results[ i ]

                results_table.append( ("Name", "MR_program", "Solution", "final_Rfact", "final_Rfree", "SHELXE_CC", "#Models", "#Residues") )
            else:
                results_table.append( ("Name", "MR_program", "Solution", "final_Rfact", "final_Rfree", "SHELXE_CC" ) )

            # Assume mrbump_results are already sorted
            mrbump_results = self.d['mrbump_results'][ cluster ]
            best=None
            for i, result in enumerate( mrbump_results ):
                
                # Remember best result
                if i == 0:
                    best = mrbump_results[i]
                
                result_summary = [ result.name,
                                   result.program,
                                   result.solution,
                                   result.rfact,
                                   result.rfree,
                                   result.shelxCC,
                                ]

                if ensemble_results:
                    # MRBUMP Results have loc0_ALL_ prepended and  _UNMOD appended
                    name = result.name[9:-6]
                    result_summary += [ name2e[ name ].num_models, name2e[ name ].num_residues ]
            
                results_table.append( result_summary )
                
            # Get nicely formatted string summarising the results
            table = printTable.Table()
            summary += table.pprint_table( results_table )
            
            # Show where it happened
            summary += '\nBest results so far are in :\n\n'
            summary +=  best.resultDir
        
        return summary
        
        
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
        
 
    def process_options(self):
        """Check the options and process any preset defaults"""
        
        # First set anything that hasn't been set to its default option
        for k, v in self.defaults.iteritems():
            if self.d[k] == None:
                #if self.debug:
                #    print "Setting default value: {0} : {1}".format(k,v)
                self.d[k] = v
            else:
                if self.debug and self.d[k] != v:
                    print "Changed default value: {0} : {1}".format(k, self.d[k])
        
        # Any changes here
        if self.d['submit_qtype']:
            self.d['submit_qtype'] = self.d['submit_qtype'].upper()
        
        
        # Check if using any preset options
        if self.d['devel_mode']:
            for k, v in self.devel_mode.iteritems():
                # Set any that haven't been set
                if self.d[k] == None:
                    self.d[k] = v
                else:
                    # Already set - only overwrite if it's set to a default value, otherwise we
                    # let the user go with what they've chosen but warn
                    if self.d[k] != self.defaults[k] and self.d[k] != v  :
                        print "WARNING! Overriding devel_mode setting: {0} : {1} with user setting {2}".format( k, v, self.d[k] )
                    else:
                        # We overwrite the default with our value
                        if self.debug:
                            print "Overriding default setting: {0} : {1} with devel_mode setting {2}".format( k, self.defaults[k], v )
                        self.d[k] = v

        if self.d['quick_mode']:
            for k, v in self.quick_mode.iteritems():
                # Set any that haven't been set
                if self.d[k] == None:
                    self.d[k] = v
                else:
                    # Already set - only overwrite if it's set to a default value, otherwise we
                    # let the user go with what they've chosen but warn
                    if self.d[k] != self.defaults[k] and self.d[k] != v  :
                        print "WARNING! Overriding quick_mode setting: {0} : {1} with user setting {2}".format( k, v, self.d[k] )
                    else:
                        # We overwrite the default with our value
                        if self.debug:
                            print "Overriding default setting: {0} : {1} with quick_mode setting {2}".format( k, self.defaults[k], v )
                        self.d[k] = v
        
        return
        
    def prettify_parameters(self):
        """
        Return the parameters nicely formated as a list of strings suitable for writing out to a file
        """
        pstr = ""
        pstr +='Params Used in this Run\n\n'
        
        keys1 = ['fasta','work_dir','mtz','pdb_code']
        pstr += '---input---\n'
        for k in keys1:
            pstr += "{0}: {1}\n".format(k, self.d[k])

        keys2 = ['make_frags','rosetta_fragments_exe','frags_3mers','frags_9mers']
        pstr+= '\n---fragments---\n'
        for k in keys2:
            pstr += "{0}: {1}\n".format(k, self.d[k])
            
        keys3 = ['make_models','rosetta_path','rosetta_db']
        pstr+= '\n---modelling---\n'
        for k in keys3:
            pstr += "{0}: {1}\n".format(k, self.d[k])
        
        if self.d['use_scwrl']:
            pstr+= '\n---3rd party---\nSCWRL {0}\n'.format( self.d['scwrl'] )

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
    

if __name__ == "__main__":
    
    import sys
    import cPickle
    
    if len(sys.argv) == 2:
        pfile = sys.argv[1]
    else:
        pfile = "/opt/ample-dev1/examples/toxd-example/ROSETTA_MR_15/ample_results.pkl"
    
    f = open(pfile)
    d = cPickle.load(f)
    
    AD = AmpleOptions()
    for k,v in d.iteritems():
        AD.d[k] = v
        
    print AD.final_summary()
