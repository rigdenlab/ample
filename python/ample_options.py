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
        for k, v in vars(parser_args).iteritems():
            #print "{} | {}".format( k, v )   
            if isinstance(v,list):
                # All values are in a list
                self.d[k] = v[0]
            else:
                self.d[k] = v
            #if v == False:
            #    self.d[k] = False
        
        # end of loop
        
        for k, v in self.d.iteritems():
            print "{} | {}".format( k, v )   
        
    def write_parameter_logfile(self, filename=None):
        """
        Write out the parameters we contain
        """
        
        pass
        
#        if not filename:
#            filename="Params_used"
#         
#        Run_params = open( filename, "w")
#        
#        Run_params.write('input params\n')
#        if self.d['DEBUG']:
#        #if DEBUG == True:
#          #print var_args
#          for print_user in var_args:
#           if var_args[print_user] is not None:
#             print print_user +' : ' + str(var_args[print_user][0])
#        
#             Run_params.write(print_user +' : ' + str(var_args[print_user][0]) + '\n')
#        
#        Run_params.write('\nParams Used in this Run\n')
#        Run_params.write('\n---input---\nFasta '+FASTA+'\nRunDir '+RunDir+'\nMTZ '+MTZ+'\nname '+PDB_code+'\n')
#        Run_params.write('\n---fragments---\nMakeFrags '+str(MakeFrags)+'\n3mers '+frags_3_mers+'\n9mers '+frags_9_mers+'\n')
#        Run_params.write('\n---modelling---\nMakeModels '+str(MakeModels)+'\nROSETTA_PATH '+ROSETTA_PATH+'\n')
#        Run_params.write('ROSETTA_cluster '+ROSETTA_cluster+'\nROSETTA_DB '+ROSETTA_DB+'\nMake_fragments_exe '+Make_fragments_exe+'\n')
#        if USE_SCWRL == True:
#          Run_params.write('\n---3rd party---\nSCWRL '+SCWRL+'\n')
#        Run_params.write('\n---Missing Domain---\nall chains fasta '+domain_all_chain_fasta+'\nall chain pdb '+domain_all_chains_pdb+'\nMISSING DOMAINS='+str(MISSING_DOMAINS)+'\n')
#        Run_params.write('Is an Insert Domain '+str(INSERT_DOMAIN)+ ' termini distance '+ str(domain_termini_distance) +'\n')
#        
#        Run_params.close()
              
    #def __str__(self):
    #    return "Ample Options"