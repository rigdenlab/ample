'''
Created on 28 Feb 2013

@author: jmht
'''
import unittest

def mrbump_cmd( adict, jobid=None, ensemble_pdb=None, fixed_iden=0.6 ):
    """
    Return a MRBUMP input string based on the variables in the adict dictionary and given keyword parameters
    
    jmht - check fixed_iden - 0.6 if not specified
    """
    
    if amopt.d['domain_all_chains_fasta']:
        fasta = amopt.d['domain_all_chains_fasta']
    else:
        fasta = amopt.d['fasta']
    
    mrs = ""
    mrs+='mrbump HKLIN {} SEQIN {} HKLOUT OUT.mtz  XYZOUT OUT.pdb << eof\n'.format( adict['mtz'], fasta )
    mrs+='LABIN SIGF={} F={} FreeR_flag={}\n'.format( adict['SIGF'], adict['F'], adict['FREE'] )
    mrs+='JOBID {}_mrbump\n'.format( jobid )
    mrs+='MRPROGRAM {}\n'.format( adict['mrbump_programs'] )
    mrs+='LOCALFILE {} CHAIN ALL RMS 1.2\n'.format( ensemble_pdb )
    mrs+='SCOPSEARCH False\n'
    mrs+='PQSSEARCH False\n'
    mrs+='SSMSEARCH False\n'
    mrs+='{}\n'.format( adict['ASU'] )
    mrs+='FAST False\n'
    mrs+='DOFASTA False\n'
    mrs+='MDLD False\n'
    mrs+='MDLC False\n'
    mrs+='MDLM False\n'
    mrs+='MDLP False\n'
    mrs+='MDLS False\n'
    mrs+='MDLU True\n'
    mrs+='UPDATE False\n'
    mrs+='BUCC  {}\n'.format( adict['use_buccaneer'] )
    mrs+='BCYCLES  {}\n'.format( adict['buccaneer_cycles'] )
    mrs+='ARPWARP  {}\n'.format( adict['use_arpwarp'] )
    mrs+='ACYCLES  {}\n'.format( adict['arpwarp_cycles'] )
    mrs+='SHELXE  {}\n'.format( adict['use_shelxe'] )
    mrs+='SCYCLES  {}\n'.format( adict['shelx_cycles'] )
    mrs+='FIXSG True\n'
    mrs+='PJOBS 1\n'
    mrs+='CHECK False\n'
    mrs+='LITE True\n'
    mrs+='PICKLE False\n'
    mrs+='TRYALL True\n'
    mrs+='USEACORN False\n'
    mrs+='USEENSEM False\n'
    mrs+='CLEAN False\n'
    mrs+='DEBUG True\n'
    
    for k in adict['mr_keys']:
        mrs+='{}\n'.format(k)
        
    if adict['domain_all_chains_pdb']:
        mrs+='FIXED_XYZIN {} IDEN {}\n'.format( adict['domain_all_chains_pdb'], fixed_iden )
        
    mrs+='END\n'
    mrs+='eof'
    
    return mrs

class Test(unittest.TestCase):


    def testName(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()