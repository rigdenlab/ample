#!/usr/bin/env ccp4-python

class LsqmanLogParser(object):
    """ Class to mine information from a lsqman log """

    def __init__(self):
        self.mc_rho = 0.
        self.rmsd = 0.
        self.rmsd_norm = 0.
        self.rmsd_rel = 0.
        self.version = None

    def parse(self, logfile):
        with open(logfile, 'r') as f:
            for line in iter(f.readline, ''):
                line = line.strip() # Odd spaces at beginning of lines
                if line.startswith("Maiorov-Crippen"): self.mc_rho=float(line.split()[-1])
                if line.startswith("Relative RMSD"): self.rmsd_rel=float(line.split()[-1])
                if line.startswith("Normalised RMSD"): self.rmsd=float(line.split()[-2])
                if line.startswith("RMSD / Nalign"): self.rmsd_norm=float(line.split()[-2])
                if line.startswith("Version"): self.version=line.split()[-1]
        return
##End LsqmanLogParser
