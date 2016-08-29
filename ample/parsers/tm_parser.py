

__author__ = "Felix Simkovic"
__data__ = "29.06.2016"
__version__ = "1.0"


class TMalignLogParser(object):
    """
    Class to mine information from a TMalign log

    Attributes
    ----------
    tm : float
       TemplateModelling score
    rmsd : float
       RMSD score
    nr_residues_common : int
       Number of residues in common between structures
    seq_id : float
       Sequence identity [n_identical/n_aligned]

    Warnings
    --------
    TMalign parser cannot distinguish which TM-score to use. Please run TMalign with the ``-a`` flag.
    """

    def __init__(self):
        self.tm = 0.0
        self.rmsd = 0.0
        self.nr_residues_common = 0
        self.seq_id = 0.0
        return

    def parse(self, logfile):
        """
        Read a TMalign logfile to extract the scores

        Parameters
        ----------
        logfile : str
           Path to the TMalign logfile
        """
        with open(logfile, 'r') as f:
            for line in iter(f.readline, ''):
                if line.startswith("Aligned length"):
                    _line = line.split()
                    nr_residues_common = int(_line[2][:-1])
                    rmsd = float(_line[4][:-1])
                    seq_id = float(_line[-1])
                if line.startswith("TM-score"):
                    tm = float(line.split()[1])

        self.set(tm, seq_id, rmsd, nr_residues_common)
        return

    def set(self, tm, seq_id, rmsd, nr_residues_common):
        """
        Manually set the TMalign parser parameters
        """
        self.tm = tm
        self.seq_id = seq_id
        self.rmsd = rmsd
        self.nr_residues_common = nr_residues_common
        return

    def reset(self):
        """
        Reset the TMalign parser
        """
        self.tm = 0.0
        self.seq_id = 0.0
        self.rmsd = 0.0
        self.nr_residues_common = 0
        return


class TMscoreLogParser(object):
    """
    Class to mine information from a TMscore log

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
