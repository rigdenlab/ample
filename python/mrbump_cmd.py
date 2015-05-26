'''
Created on 28 Feb 2013

@author: jmht
'''
import os
import sys

def mrbump_cmd(adict,jobid,keyinFile):
    """Return the command to run mrbump
    Need to return the full path as on *!*windoze*!* the mrbump script isn't executable
    """
    if sys.platform.startswith("win"):
        ccp4python=os.path.join(os.environ["CCP4"],"bin","ccp4-python")
        mrbump=os.path.join(os.environ["CCP4"],"bin","mrbump")
        mrbump="{0} {1}".format(ccp4python,mrbump)
    else:
        mrbump='mrbump'
    return'{0} KEYIN {1} HKLIN {2} SEQIN {3} HKLOUT {4}.mtz  XYZOUT {4}.pdb'.format(mrbump,
                                                                                    keyinFile,
                                                                                    adict['mtz'],
                                                                                    adict['mr_sequence'],
                                                                                    jobid)

def keyword_dict(amoptd, extra_options={}):
    keywords = [
                'arpwarp_cycles',
                'buccaneer_cycles',
                'debug',
                'domain_all_chains_pdb',
                'F',
                'FREE',
                'mr_keys',
                'mrbump_programs',
                'native_pdb',
                'nmasu',
                'phaser_kill',
                'shelx_cycles',
                'shelxe_exe',
                'shelxe_rebuild_arpwarp',
                'shelxe_rebuild_buccaneer',
                'SIGF',
                'use_arpwarp',
                'use_buccaneer',
                'use_shelxe',
                ]
    
    # Pull out all mrbump options from the main ample dict
    key_dict = dict( (k,v) for k,v in amoptd.iteritems() if k in keywords )
    
    # Change any options for this ensemble
    for k,v in extra_options.iteritems(): key_dict[k] = v
    
    return key_dict

def mrbump_keywords(adict=None, jobid=None, ensemble_pdb=None, fixed_iden=0.6):
    """
    Create MRBUMP keywords
    
    Args:
    adict -- dictionary of options
    
    jmht - check fixed_iden - 0.6 if not specified
    """
    
    #mrs+='mrbump HKLIN {0} SEQIN {1} HKLOUT OUT.mtz  XYZOUT OUT.pdb << eof\n'.format( adict['mtz'], adict['mr_sequence'] )
    mrs='LABIN SIGF={0} F={1} FreeR_flag={2}\n'.format( adict['SIGF'], adict['F'], adict['FREE'] )
    mrs+='JOBID {0}_mrbump\n'.format( jobid )
    mrs+='MRPROGRAM {0}\n'.format( " ".join( adict['mrbump_programs'] ) )
    mrs+='LOCALFILE {0} CHAIN ALL RMS 0.1\n'.format( ensemble_pdb )
    #
    # Don't do any of the searches as we are providing a local file
    #
    mrs+='SCOPSEARCH False\n'
    mrs+='PQSSEARCH False\n'
    mrs+='SSMSEARCH False\n'
    mrs+='DOFASTA False\n'
    mrs+='DOPHMMER False\n'
    mrs+='DOHHPRED False\n'
    #
    mrs+='FAST False\n'
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
    if adict['shelxe_rebuild_arpwarp'] or adict['shelxe_rebuild_buccaneer']:
        # Rebuild SHELXE trace with both Buccaneer and ArpWarp
        mrs+='SXREBUILD True\n'
        if adict['shelxe_rebuild_buccaneer']: mrs+='SXRBUCC True\n'
        if adict['shelxe_rebuild_buccaneer']: mrs+='SXRARPW True\n'
    if adict['nmasu'] > 0:
        mrs+='NMASU  {0}\n'.format( adict['nmasu'] )
    if adict['domain_all_chains_pdb']:
        mrs+='FIXED_XYZIN {0} IDEN {1}\n'.format( adict['domain_all_chains_pdb'], fixed_iden )
    if adict['native_pdb']:
        mrs+='PDBNATIVE {0}\n'.format( adict['native_pdb'] )
    if adict['phaser_kill'] > 0:
        mrs+='PKEY KILL TIME {0}\n'.format(adict['phaser_kill'])
        
    
    # Extra keywords
    # This assumes everything in mr_keys is a list of [ KEYWORD, VALUE0, VALUE1, ...]
    for l in adict['mr_keys']:
        mrs += "  ".join( l ) + "\n"
    #for i in range( 0, len(adict['mr_keys']), 2 ):
    #    mrs+='{0}  {1}\n'.format( adict['mr_keys'][i], adict['mr_keys'][i+1]  )
        
    mrs+='END\n'
    #mrs+='eof'
    
    return mrs
