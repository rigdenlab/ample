'''
Created on 28 Feb 2013

@author: jmht
'''
def mrbump_cmd( adict, jobid=None, ensemble_pdb=None, fixed_iden=0.6 ):
    """
    Create MRBUMP input string based on given variables
    
    Args:
    adict -- amopt dictionary
    
    jmht - check fixed_iden - 0.6 if not specified
    """
    
    mrs = ""
    mrs+='mrbump HKLIN {0} SEQIN {1} HKLOUT OUT.mtz  XYZOUT OUT.pdb << eof\n'.format( adict['mtz'], adict['mr_sequence'] )
    mrs+='LABIN SIGF={0} F={1} FreeR_flag={2}\n'.format( adict['SIGF'], adict['F'], adict['FREE'] )
    mrs+='JOBID {0}_mrbump\n'.format( jobid )
    #mrs+='MRPROGRAM {0}\n'.format( adict['mrbump_programs'] )
    mrs+='MRPROGRAM {0}\n'.format( " ".join( adict['mrbump_programs'] ) )
    mrs+='LOCALFILE {0} CHAIN ALL RMS 0.1\n'.format( ensemble_pdb )
    mrs+='SCOPSEARCH False\n'
    mrs+='PQSSEARCH False\n'
    mrs+='SSMSEARCH False\n'
    mrs+='FAST False\n'
    mrs+='DOFASTA False\n'
    mrs+='MDLD False\n'
    mrs+='MDLC False\n'
    mrs+='MDLM False\n'
    mrs+='MDLP False\n'
    mrs+='MDLS False\n'
    mrs+='MDLU True\n'
    mrs+='UPDATE False\n'
    mrs+='BUCC  {0}\n'.format( adict['use_buccaneer'] )
    mrs+='BCYCLES  {0}\n'.format( adict['buccaneer_cycles'] )
    mrs+='ARPWARP  {0}\n'.format( adict['use_arpwarp'] )
    mrs+='ACYCLES  {0}\n'.format( adict['arpwarp_cycles'] )
    mrs+='SHELXE  {0}\n'.format( adict['use_shelxe'] )
    mrs+='SHLXEXE  {0}\n'.format( adict['shelxe_exe'] )
    mrs+='SCYCLES  {0}\n'.format( adict['shelx_cycles'] )
    mrs+='FIXSG True\n'
    mrs+='PJOBS 1\n'
    mrs+='CHECK False\n'
    mrs+='LITE True\n'
    mrs+='PICKLE False\n'
    mrs+='TRYALL True\n'
    mrs+='USEACORN False\n'
    mrs+='USEENSEM False\n'
    mrs+='CLEAN False\n'
    mrs+='DEBUG {0}\n'.format( adict['debug'] )
    #
    # Optional extras
    #
    if adict['shelxe_rebuild']:
        # Rebuild SHELXE trace with both Buccaneer and ArpWarp
        mrs+='SXREBUILD True\n'
        mrs+='SXRBUCC True\n'
        mrs+='SXRARPW True\n'
    if adict['ASU'] > 0:
        mrs+='NMASU  {0}\n'.format( adict['ASU'] )
    if adict['domain_all_chains_pdb']:
        mrs+='FIXED_XYZIN {0} IDEN {1}\n'.format( adict['domain_all_chains_pdb'], fixed_iden )
    
    # Extra keywords
    # This assumes everything in mr_keys is a list of [ KEYWORD, VALUE0, VALUE1, ...]
    for l in adict['mr_keys']:
        mrs += "  ".join( l ) + "\n"
    #for i in range( 0, len(adict['mr_keys']), 2 ):
    #    mrs+='{0}  {1}\n'.format( adict['mr_keys'][i], adict['mr_keys'][i+1]  )
        
    mrs+='END\n'
    mrs+='eof'
    
    return mrs
