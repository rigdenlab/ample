

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
                    nr_residues_common = int(line.split()[5])
                if line.startswith("RMSD"):
                    rmsd = float(line.split()[5])
                if line.startswith("TM-score"):
                    tm = float(line.split()[2])
                if line.startswith("MaxSub-score"):
                    maxsub = float(line.split()[1])
                if line.startswith("GDT-TS-score"):
                    gdtts = float(line.split()[1])
                if line.startswith("GDT-HA-score"):
                    gdtha = float(line.split()[1])

        self.set(tm, maxsub, gdtts, gdtha, rmsd, nr_residues_common)

        return

    def set(self, tm, maxsub, gdtts, gdtha, rmsd, nr_residues_common):
        """
        Manually set the TMscore parser parameters
        """
        self.tm = tm
        self.maxsub = maxsub
        self.gdtts = gdtts
        self.gdtha = gdtha
        self.rmsd = rmsd
        self.nr_residues_common = nr_residues_common
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
        self.nr_residues_common = 0
        return
