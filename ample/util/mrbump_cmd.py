'''
Created on 28 Feb 2013

@author: jmht
'''
import os
import sys

from ample.util import ample_util

def mrbump_cmd(name, mtz, mr_sequence, keyword_file):
    """Return the command to run mrbump"""
    if sys.platform.startswith("win"):
        mrbump = os.path.join(os.environ["CCP4"], "bin", "mrbump" + ample_util.SCRIPT_EXT)
    else:
        mrbump = os.path.join(os.environ["CCP4"], "bin", "mrbump")
    cmd = [
        mrbump,
        "KEYIN", "{0}".format(keyword_file),
        "HKLIN", "{0}".format(mtz),
        "SEQIN", "{0}".format(mr_sequence),
        "HKLOUT", "{0}.mtz".format(name),
        "XYZOUT", "{0}.pdb".format(name),
    ]
    return " ".join(cmd)

def keyword_dict(ensemble_pdb, name, amoptd, extra_options={}):
    """Extract the mrbump keywords from the main ample dictionary and add/change any from
    the extra_options dict"""
    keywords = [
                'arpwarp_cycles',
                'buccaneer_cycles',
                'debug',
                'domain_all_chains_pdb',
                'F',
                'FREE',
                'mr_keys',
                'mr_sg_all',
                'mrbump_programs',
                'native_pdb',
                'nmasu',
                'phaser_kill',
                'phaser_rms',
                'shelx_cycles',
                'shelxe_exe',
                'shelxe_rebuild_arpwarp',
                'shelxe_rebuild_buccaneer',
                'SIGF',
                'refine_rebuild_arpwarp',
                'refine_rebuild_buccaneer',
                'use_shelxe',
                ]
    
    # Pull out all mrbump options from the main ample dict
    key_dict = dict((k, v) for k, v in amoptd.iteritems() if k in keywords)
    
    # Change any/add options for this ensemble
    for k, v in extra_options.iteritems(): key_dict[k] = v
    
    # Add ensemble_pdb and name
    key_dict['name'] = name
    key_dict['ensemble_pdb'] = ensemble_pdb
    return key_dict

def mrbump_keyword_file(odict, fixed_iden=0.6):
    """
    Create MRBUMP keywords
    
    Args:
    odict -- dictionary of options
    
    jmht - check fixed_iden - 0.6 if not specified
    """
    mrs = 'LABIN SIGF={0} F={1} FreeR_flag={2}\n'.format(odict['SIGF'], odict['F'], odict['FREE'])
    mrs += 'JOBID {0}_mrbump\n'.format(odict['name'])
    mrs += 'MRPROGRAM {0}\n'.format(" ".join(odict['mrbump_programs']))
    mrs += 'LOCALFILE {0} CHAIN ALL RMS {1}'.format((odict['ensemble_pdb']), odict['phaser_rms'])
    if 'ncopies' in odict and odict['ncopies'] > 0: mrs += ' COPIES {0}'.format(odict['ncopies'])
    mrs += '\n'
    #
    # Don't do any of the searches as we are providing a local file
    #
    mrs += 'SCOPSEARCH False\n'
    mrs += 'PQSSEARCH False\n'
    mrs += 'SSMSEARCH False\n'
    mrs += 'DOFASTA False\n'
    mrs += 'DOPHMMER False\n'
    mrs += 'DOHHPRED False\n'
    #
    mrs += 'FAST False\n'
    mrs += 'MDLD False\n'
    mrs += 'MDLC False\n'
    mrs += 'MDLM False\n'
    mrs += 'MDLP False\n'
    mrs += 'MDLS False\n'
    mrs += 'MDLU True\n'
    mrs += 'UPDATE False\n'
    mrs += 'BUCC  {0}\n'.format(odict['use_buccaneer'])
    mrs += 'BCYCLES  {0}\n'.format(odict['buccaneer_cycles'])
    mrs += 'ARPWARP  {0}\n'.format(odict['refine_rebuild_arpwarp'])
    mrs += 'ACYCLES  {0}\n'.format(odict['arpwarp_cycles'])
    mrs += 'SHELXE  {0}\n'.format(odict['use_shelxe'])
    mrs += 'SHLXEXE  {0}\n'.format(odict['shelxe_exe'])
    mrs += 'SCYCLES  {0}\n'.format(odict['shelx_cycles'])
    mrs += 'FIXSG True\n'
    mrs += 'PJOBS 1\n'
    mrs += 'CHECK False\n'
    mrs += 'LITE True\n'
    mrs += 'PICKLE False\n'
    mrs += 'TRYALL True\n'
    mrs += 'USEACORN False\n'
    mrs += 'USEENSEM False\n'
    mrs += 'CLEAN False\n'
    mrs += 'DEBUG {0}\n'.format(odict['debug'])
    
    #
    # Optional extras
    #
    if odict['shelxe_rebuild_arpwarp'] or odict['shelxe_rebuild_buccaneer']:
        # Rebuild SHELXE trace with both Buccaneer and ArpWarp
        mrs += 'SXREBUILD True\n'
        if odict['shelxe_rebuild_buccaneer']: mrs += 'SXRBUCC True\n'
        if odict['shelxe_rebuild_buccaneer']: mrs += 'SXRARPW True\n'
    if odict['nmasu'] > 0:
        mrs += 'NMASU  {0}\n'.format(odict['nmasu'])
    if odict['domain_all_chains_pdb']:
        mrs += 'FIXED_XYZIN {0} IDEN {1}\n'.format(odict['domain_all_chains_pdb'], fixed_iden)
    if odict['native_pdb']:
        mrs += 'PDBNATIVE {0}\n'.format(odict['native_pdb'])
    if odict['phaser_kill'] > 0:
        mrs += 'PKEY KILL TIME {0}\n'.format(odict['phaser_kill'])
    if odict['mr_sg_all']:
        mrs += 'PKEY SGALTERNATIVE SELECT ALL\n'
    
    # Extra keywords
    # This assumes everything in mr_keys is a list of [ KEYWORD, VALUE0, VALUE1, ...]
    if odict['mr_keys']:
        for l in odict['mr_keys']: mrs += "  ".join(l) + "\n"
    mrs += 'END\n'
    return mrs
