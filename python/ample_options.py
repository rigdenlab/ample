'''
Class to hold the options for ample
'''

class AmpleOptions(object):
    
    def __init__(self):
        
        # The dictionary with all the options
        self.d = {}
        pass
        
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
            
            #if v == False:
            #    self.d[k] = False
        
        # end of loop
        
        #for k, v in self.d.iteritems():
        #    print "{} | {}".format( k, v )
        
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

        keys2 = ['make_frags','rosetta_fragments_exe','frags3mers','frags9mers']
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
                pstr += "{} : {}\n".format( k, v )
                
        return pstr
