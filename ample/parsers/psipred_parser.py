"""A parser object for Psipred Ss2 files"""

from __future__ import print_function

__author__ = "Felix Simkovic & Jens Thomas"
__date__ = "13 Jan 2016"
__version__ = "0.1"

import collections
import logging
import warnings

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


class PsipredSs2Parser(object):
    """Parser for psipred ss2 file"""

    def __init__(self, ss2file=None):
        self.residues = None
        if ss2file:
            self.parse(ss2file)

    @property
    def secondary_structure(self):
        """The secondary structure

        Returns
        -------
        str
           The secondary structure one-letter codes

        """
        if not self.residues:
            return None
        return "".join([i.ss for i in self.residues])

    def parse(self, ss2file):
        """Parse a secondary structure file

        Parameters
        ----------
        ss2file : str
           The path to the Psipred ss2 file

        """
        PSIPredResidueInfo = collections.namedtuple(
            "PSIPredResidueInfo", ["rank", "residue", "ss", "coil", "helix", "strand"]
        )
        residues = []

        with open(ss2file, 'r') as fhin:
            for line in iter(fhin.readline, ''):
                if line[0] == '#' or not line.strip():
                    continue

                line = line.split()

                rank = int(line[0])
                residue = line[1]
                ss = line[2]
                coil, helix, strand = map(float, line[3:6])

                residues.append(
                    PSIPredResidueInfo(rank=rank, residue=residue, ss=ss, coil=coil, helix=helix, strand=strand)
                )

        self.residues = tuple(residues)
        return

    def check_content(self):
        """Check the secondary structure composition"""

        H = len([i for i in self.residues if i.ss == "H"])
        E = len([i for i in self.residues if i.ss == "E"])

        if H > 0 and E > 0:
            logging.info('Your protein is predicted to be mixed alpha beta, your chances of success are intermediate')
        if H == 0 and E > 0:
            logging.info('Your protein is predicted to be all beta, your chances of success are low')
        if H > 0 and E == 0:
            logging.info('Your protein is predicted to be all alpha, your chances of success are high')
        if H == 0 and E == 0:
            logging.info('Your protein is has no predicted secondary structure, your chances of success are low')
        return

    def checkContent(self):
        """Check the secondary structure composition"""
        warnings.warn(
            DeprecationWarning, "This function will be removed in a future release - use check_content() instead"
        )
        return self.check_content()

    def getSecondaryStructure(self):
        """Get the secondary structure content"""
        warnings.warn(
            DeprecationWarning,
            "This function will be removed in a future release - use attribute secondary_structure instead",
        )
        return self.secondary_structure
