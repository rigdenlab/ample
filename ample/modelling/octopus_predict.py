"""
Created on 18 Feb 2013

@author: jmht

Query the octopus server http://octopus.cbr.su.se to get transmembrane predictions
"""

import logging
import os
import sys
import urllib

if sys.version_info.major < 3:
    from HTMLParser import HTMLParser
    from urllib2 import urlopen, HTTPError
else:
    from html.parser import HTMLParser
    from urllib.error import HTTPError
    from urllib.request import urlopen

TIMEOUT = 10.0  # in seconds


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
        self.writeFiles(name, directory)

    def octopusFileUrls(self, fasta):
        """
        Query the server for a prediction for the given fasta sequence.
        
        Args:
        fasta -- a single fasta sequence as a string
        
        Sets the urls of the topo and nnprf files
        """

        data = dict(do='Submit OCTOPUS', sequence=fasta)
        edata = urllib.urlencode(data)

        req = urlopen(self.octopus_url, edata, timeout=TIMEOUT)

        m = ParseFileUrl()
        # Calls handle_starttag, _data etc.
        html = req.read()
        m.feed(html)

        if not m.topo:
            with open("OCTOPUS_ERROR.html", "w") as f:
                f.write(html)
            raise RuntimeError(
                "Error getting prediction for fasta:{}\nCheck file: OCTOPUS_ERROR.html".format(fasta.splitlines()[0])
            )

        self.topo_url = self.octopus_url + m.topo
        self.nnprf_url = self.octopus_url + m.nnprf

    def writeFiles(self, name, directory):
        """
        Write the files on the server to disk
        
        Args:
        directory: where to write files
        name: name for files (with suffix .topo and .nnprf)
        """
        try:
            topo_req = urlopen(self.topo_url, timeout=TIMEOUT)
        except HTTPError as e:
            raise RuntimeError("Error accessing topo file: {}\n{}".format(self.topo_url, e))

        try:
            nnprf_req = urlopen(self.nnprf_url, timeout=TIMEOUT)
        except uHTTPError as e:
            msg = "Error accessing nnprf file: {}\n{}\nTransmembrane prediction may have failed!".format(
                self.nnprf_url, e
            )
            self.logger.warn(msg)

        fname = os.path.join(directory, name + ".topo")
        with open(fname, "w") as f:
            f.writelines(topo_req.readlines())
        self.logger.debug("Wrote topo file: {}".format(fname))
        self.topo = fname

        fname = os.path.join(directory, name + ".nnprf")
        with open(fname, "w") as f:
            f.writelines(nnprf_req.readlines())
        self.logger.debug("Wrote nnprf file: {}".format(fname))
        self.nnprf = fname

    def getFasta(self, fastafile):
        """Given a fastafile, extract the first sequence and return it as \n separated string"""
        fasta = []
        header = False
        with open(fastafile, "r") as f:
            for line in f:
                line = line.strip()
                # skip empty
                if not len(line):
                    continue
                # Only read one sequence
                if line.startswith(">"):
                    if header:
                        break
                    header = True
                fasta.append(line)

        if not len(fasta):
            return None

        return "\n".join(fasta)
