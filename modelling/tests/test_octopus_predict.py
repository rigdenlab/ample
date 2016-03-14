
import os
import unittest
from ample.modelling import octopus_predict


class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = cls.thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -2 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"testing")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        return

    def test_get_predict(self):
        os.chdir(self.thisd) # Need as otherwise tests that happen in other directories change os.cwd()        
         
        fastafile = os.path.join(self.testfiles_dir ,"2uui.fasta")
        octo = octopus_predict.OctopusPredict()
        fasta = octo.getFasta(fastafile)
        octo.getPredict("2uui",fasta)
     
        self.assertIsNotNone(octo.topo, "Error getting topo file")
    
        ref="""##############################################################################
OCTOPUS result file
Generated from http://octopus.cbr.su.se/ at 2014-11-17 18:15:38
Total request time: 2.67 seconds.
##############################################################################


Sequence name: 2uui_A; molId:1; molType:protein; unp:Q16873; molName:LEUKOTRIENE C4 SY...
Sequence length: 156 aa.
Sequence:
MHHHHHHKDEVALLAAVTLLGVLLQAYFSLQVISARRAFRVSPPLTTGPPEFERVYRAQV
NCSEYFPLFLATLWVAGIFFHEGAAALCGLVYLFARLRYFQGYARSAQLRLAPLYASARA
LWLLVALAALGLLAHFLPAALRAALLGRLRTLLPWA

OCTOPUS predicted topology:
oooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiiiiiiii
iiiiiiiMMMMMMMMMMMMMMMoMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiMMMMM
MMMMMMMMMMMMMMMMoooooooooooooooooooo
"""
        
        with open(octo.topo) as f:
            lines=[l.strip() for l in f]
        
        self.assertEqual(lines[7:],ref.split("\n")[7:])
        
        os.unlink("2uui.topo")
        os.unlink("2uui.nnprf")

if __name__ == "__main__":
    unittest.main()
