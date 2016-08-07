#!/usr/bin/env ccp4-python

import os

from ample.util import ample_util
from ample.util import mtz_util

class Cphasematch( object ):
    
    def __init__(self):
        self.resolution = None
        self.F = None
        self.SIGF= None
        self.name = None
        self.cad_logfile = None
        self.cphasematch_logfile = None
        self.before_origin = None
        self.after_origin = None
        
    def run(self, native_MTZ, placed_MTZ):
        self.run_cad(native_MTZ, placed_MTZ)
        self.run_cphasematch()
        self.return_results()
        return
    
    def run_cad(self, native_MTZ, placed_MTZ, cleanup=True):
        
        self.F, self.SIGF, _ = mtz_util.get_labels(native_MTZ)
        self.resolution = mtz_util.get_resolution(native_MTZ)
        self.name = os.path.basename(placed_MTZ)[16:33]
        self.cad_logfile = self.name + "_cad.log"
        
        
        cmd= [ 'cad',
              "hklin1",
              native_MTZ,
              "hklin2",
              placed_MTZ,
              "hklout",
              "all_phases.mtz"
              ]
        
        key = """LABIN FILE_NUMBER 1 E1={0} E2={1} E3=FC E4=PHIC
LABIN FILE_NUMBER 2 E1=FC E2=PHIC
LABOUT FILE_NUMBER 2 E1=FCalc E2=PHICalc
RESOLUTION OVERALL 1000.0 {2}
END""".format(self.F, self.SIGF, self.resolution)

        retcode = ample_util.run_command(cmd=cmd, stdin=key, logfile=self.cad_logfile)
         
        if retcode != 0: raise RuntimeError, "Error running command: {0}".format( " ".join(cmd) )
        
        if cleanup: os.unlink(self.cad_logfile)
        return
    
    def run_cphasematch(self, cleanup=False):
        
        self.cphasematch_logfile = self.name + "_cphasematch.log"
        
        cmd= [ 'cphasematch',
              "-stdin",
              ]
        
        key = """mtzin all_phases.mtz
mtzout cphasematch.mtz
colin-fo /*/*/[{0},{1}]
colin-fc-1 /*/*/[FC,PHIC]
colin-fc-2 /*/*/[FCalc,PHICalc]
resolution-bins 12""".format(self.F, self.SIGF)

        retcode = ample_util.run_command(cmd=cmd, stdin=key, logfile=self.cphasematch_logfile)
         
        if retcode != 0: raise RuntimeError, "Error running command: {0}".format( " ".join(cmd) )
        
        if cleanup: os.unlink(self.cphasematch_logfile)
        return
    
    def return_results(self):
        with open(self.cphasematch_logfile, 'r') as f:
            for line in f:
                if line.startswith(' Mean phase error before origin fixing:'):
                    self.before_origin = line.split()[-1]
                elif line.startswith(' Mean phase error after  origin fixing:'):
                    self.after_origin = line.split()[-1]
        return
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Run cphasematch on two mtz files', prefix_chars="-")

    group = parser.add_argument_group()
    group.add_argument('-native_mtz',
                       help="Input native MTZ file")
    group.add_argument('-placed_mtz',
                       help="Input placed MTZ file")

    args = parser.parse_args()

    if args.native_mtz and args.placed_mtz:
        cp = Cphasematch()
        cp.run_cad(args.native_mtz, args.placed_mtz)
        cp.run_cphasematch()
        cp.return_results()
