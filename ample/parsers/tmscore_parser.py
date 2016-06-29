

__author__ = "Felix Simkovic"
__data__ = "29.06.2016"
__version__ = "1.0"


class TMscoreLogParser(object):
    """
    Class to mine information from a tmscore log

    Attributes
    ----------
    tm : float
       TemplateModelling score
    maxsub : float
       MaxSub score
    gdtts : float
       GDT-TS score
    gdtha : float
       GDT-HA score
    rmsd : float
       RMSD score
    nr_residues_common : int
       Number of residues in common between structures

    """

    def __init__(self):
        self.tm = 0.0
        self.maxsub = 0.0
        self.gdtts = 0.0
        self.gdtha = 0.0
        self.rmsd = 0.0
        self.nr_residues_common = 0
        return

    def parse(self, logfile):
        """
        Read a TMscore logfile to extract the scores

        Parameters
        ----------
        logfile : str
           Path to the TMscore logfile
        """
        with open(logfile, 'r') as f:
            for line in iter(f.readline, ''):
                if line.startswith("Number of residues in common"):
                    self.nrResiduesCommon=int(line.split()[5])
                if line.startswith("RMSD"):
                    self.rmsd = float(line.split()[5])
                if line.startswith("TM-score"):
                    self.tm = float(line.split()[2])
                if line.startswith("MaxSub-score"):
                    self.maxsub = float(line.split()[1])
                if line.startswith("GDT-TS-score"):
                    self.gdtts = float(line.split()[1])
                if line.startswith("GDT-HA-score"):
                    self.gdtha = float(line.split()[1])

        return

    def reset(self):
        """
        Reset the TMscore parser
        """
        self.tm = 0.0
        self.maxsub = 0.0
        self.gdtts = 0.0
        self.gdtha = 0.0
        self.rmsd = 0.0
        self.nrResiduesCommon = 0
