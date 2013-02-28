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
        pstr +='Params Used in this Run\n'
        
        tkeys = ['fasta','work_dir','mtz','pdb_code']
        pstr += '---input---\n'
        for k in tkeys:
            pstr += "{}: {}\n".format(k, self.d[k])

        tkeys = ['make_frags','rosetta_fragments_exe','frags3mers','frags9mers','pdb_code']
        pstr+= '---fragments---\n'
        for k in tkeys:
            pstr += "{}: {}\n".format(k, self.d[k])
            
        tkeys = ['make_models','rosetta_path','rosetta_db','pdb_code']
        pstr+= '---modelling---\n'
        for k in tkeys:
            pstr += "{}: {}\n".format(k, self.d[k])
        
        if self.d['use_scwrl']:
            pstr+= '\n---3rd party---\nSCWRL {}\n'.format( self.d['scwrl'] )

        tkeys = ['domain_all_chains_fasta','domain_all_chains_pdb']
        if tkeys[0] or tkeys[1]:
            pstr+= '---Missing Domain---\n'
            for k in tkeys:
                pstr += "{}: {}\n".format(k, self.d[k])
        
        # This only used for printing
        INSERT_DOMAIN = False
        if self.d['domain_termini_distance'] > 0:
            INSERT_DOMAIN = True
        pstr += 'Is an Insert Domain {} termini distance {}\n'.format( INSERT_DOMAIN, self.d['domain_termini_distance'] )
        
        return pstr
