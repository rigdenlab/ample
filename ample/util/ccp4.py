"""CCP4 wrapper"""

__author__ = "Felix Simkovic"
__date__ = "29 Jan 2018"
__version__ = "1.0"

from distutils.version import StrictVersion

import os


class CCP4(object):
    """Wrapper class for CCP4 installation"""
    def __init__(self):
        self.email = "ccp4@stfc.ac.uk"
        self.root = CCP4RootDirectory()
        self.version = CCP4Version()


class CCP4RootDirectory(object):
    """The CCP4 root directory"""
    def __init__(self):
        if "CCP4" not in os.environ:
            raise KeyError("Cannot find CCP4 installation - please make sure CCP4 "
                           + "is installed and the setup scripts have been run!")
        elif "CCP4_SCR" not in os.environ:
            raise KeyError("$CCP4_SCR environment variable not set - please make sure "
                           + "CCP4 is installed and the setup scripts have been run!")
        elif not os.path.isdir(os.environ['CCP4_SCR']):
            raise ValueError("Cannot find the $CCP4_SCR directory: {0}".format(os.environ["CCP4_SCR"]))
        else:
            self._root = os.environ['CCP4']

    def __str__(self):
        return self._root

    def __repr__(self):
        return "{}: {}".format(self.__class__.__name__, self._root)


class CCP4Version(StrictVersion):
    """The CCP4 version class"""
    def __init__(self):
        StrictVersion.__init__(self)
        ccp4_major_minor = os.path.join(os.environ["CCP4"], "lib", "ccp4", "MAJOR_MINOR")
        if os.path.isfile(ccp4_major_minor):
            with open(ccp4_major_minor, "r") as f_in:
                tversion = f_in.read().strip()
        else:
            stdout = cexec(['pdbcur' + EXE_EXT], permit_nonzero=True)
            tversion = None
            for line in stdout.split(os.linesep):
                if line.startswith(' ### CCP4'):
                    tversion = line.split()[2].rstrip(':')
            if tversion is None:
                raise RuntimeError("Cannot determine CCP4 version")
        self.parse(tversion)
