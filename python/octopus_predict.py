'''
Created on 18 Feb 2013

@author: jmht

Query the octopus server http://octopus.cbr.su.se to get transmembrane predictions
'''

import sys
import os

from HTMLParser import HTMLParser
import logging
import urllib
import urllib2
import unittest


class ParseFileUrl(HTMLParser):
    """
    Parse the page returned by an octopus search to get the links to the files
    """
    
    
    def __init__(self):
        self.topo = None
        self.nnprf = None
        HTMLParser.__init__(self)
    
    def handle_starttag(self, tag, attrs):
        """Set recording whenever we encounter a tag so handle_data can process it"""
        if tag == 'a' and attrs:
            for name, value in attrs:
                if name == "href" and str(value).endswith(".topo"): 
                    self.topo = str(value)
                elif name == "href" and str(value).endswith(".nnprf"):
                    self.nnprf = str(value)
# End ParseFileUrl

class OctopusPredict(object):
    """
    Query the octopus server http://octopus.cbr.su.se to get transmembrane predictions
    """
    
    def __init__(self):
        
        
        self.logger = logging.getLogger()
        
        self.octopus_url = "http://octopus.cbr.su.se/"
        
        # The fasta sequence to query with
        self.fasta = None
        # path to the topo file
        self.topo = None
        # path to the nnprf file
        self.nnprf = None
        
        # url of topo & nnprf files on server
        self.topo_url = None
        self.nnprf_url = None
        
    
    def getPredict(self, name, fasta, directory=None):
        """
        Get the octopus prediction for the given fasta sequence string
        
        Args:
        fasta -- fasta sequence string
        name -- name for the files
        directory -- directory to write files to
        
        Sets the topo and nnprf attributes
        """
        
        if not directory:
            directory = os.getcwd()
        
        self.octopusFileUrls(fasta)
        self.writeFiles(name, directory )
        
    def octopusFileUrls(self, fasta ):
        """
        Query the server for a prediction for the given fasta sequence.
        
        Args:
        fasta -- a single fasta sequence as a string
        
        Sets the urls of the topo and nnprf files
        """
    
        data = { 'do' : 'Submit OCTOPUS',  'sequence' : fasta }
        edata = urllib.urlencode( data )
        
        #try:
        req = urllib2.urlopen( self.octopus_url, edata )
        #except Exception, e:
        #    # Need to encode to deal with possible dodgy characters 
        #    print "Error accessing url: %s\n%s" % (url.encode('ascii','replace'),e)
        #    sys.exit(1)
        
        
        m = ParseFileUrl()
        # Calls handle_starttag, _data etc.
        html = req.read()
        m.feed( html )
        
        if not m.topo:
            f = open("OCTOPUS_ERROR.html", "w")
            f.write( html )
            f.close()
            msg = "Error getting prediction for fasta:{}\nCheck file: OCTOPUS_ERROR.html".format( fasta.splitlines()[0] )
            self.logger.critical(msg)
            raise RuntimeError, msg
        

        self.topo_url = self.octopus_url + m.topo
        self.nnprf_url = self.octopus_url + m.nnprf
        
        return
    ##End octopusFileUrls
    
    
    def writeFiles(self, name, directory ):
        """
        Write the files on the server to disk
        
        Args:
        directory: where to write files
        name: name for files (with suffix .topo and .nnprf)
        """
        # Get full file path
        try:
            topo_req = urllib2.urlopen( self.topo_url )
        except urllib2.HTTPError,e:
            msg = "Error accessing topo file: {}\n{}".format(self.topo_url,e)
            self.logger.critical(msg)
            raise RuntimeError, msg
        
        try:
            nnprf_req = urllib2.urlopen( self.nnprf_url )
        except urllib2.HTTPError,e:
            msg ="Error accessing nnprf file: {}\n{}\nTransmembrane prediction may have failed!".format(self.nnprf_url,e)
            self.logger.warn(msg)
        
        fname = os.path.join(directory, name + ".topo")
        f = open( fname, "w" )
        f.writelines( topo_req.readlines() )
        self.logger.debug("Wrote topo file: {}".format( fname ) )
        f.close()
        self.topo=fname
        
        fname = os.path.join(directory, name + ".nnprf")
        f = open( fname, "w" )
        f.writelines( nnprf_req.readlines() )
        self.logger.debug("Wrote nnprf file: {}".format( fname ) )
        f.close()
        self.nnprf = fname
        
        
    def getFasta(self, fastafile ):
        """
        Given a fastafile, extract the first sequence
        and return it as \n separated string
        """
        fasta = []
        header=False
        with open( fastafile, "r" ) as f:
            for line in f:
                line = line.strip()
                # skip empty
                if not len(line):
                    continue
                # Only read one sequence
                if line.startswith(">"):
                    if header:
                        break
                    header=True
                fasta.append(line)
            
        if not len(fasta):
            return None
        
        return "\n".join(fasta)


  
class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        cls.ample_dir = os.sep.join( paths[ : -1 ] )
        cls.tests_dir=os.path.join(cls.ample_dir,"tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir,'testfiles')
        return

    def testGetPredict(self):
        """See we can get the prediction"""
        
        logging.basicConfig()
        logging.getLogger().setLevel(logging.DEBUG)
         
        fastafile = os.path.join(self.ample_dir,"tests","testfiles","2uui.fasta")
        octo = OctopusPredict()
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
        
        #class FOO():
        #    pass
        #octo=FOO()
        #octo.topo="2uui.topo"
        with open(octo.topo) as f:
            lines=[l.strip() for l in f]
        
        self.assertEqual(lines[7:],ref.split("\n")[7:])
        
        os.unlink("2uui.topo")
        os.unlink("2uui.nnprf")
        
        return 

def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testGetPredict'))
    return suite
    
#
# Run unit tests
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(testSuite())
