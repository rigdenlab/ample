#!/usr/bin/env ccp4-python

'''
03.07.2015

@author: hlfsimko
'''

import argparse
import os
import subprocess
import sys

if not "CCP4" in os.environ.keys(): raise RuntimeError('CCP4 not found')
sys.path.insert(0, os.path.join(os.environ['CCP4'], "share", "ample", "python"))
sys.path.insert(0, os.path.join(os.environ['CCP4'], "share", "ample", "parsers"))

# Custom
import ample_util
import parse_alignment
import parse_cdhit


def cluster(cd_hit, alnfile, wdir):
    '''Run cd-hit for your alignment file'''

    in_fname    = alnfile
    out_fname   = alnfile.rsplit('.',1)[0] + '.cd'

    # -n 4 for thresholds 0.6 ~ 0.7
    # Taken from http://weizhong-lab.ucsd.edu/cd-hit/wiki/doku.php?id=cd-hit_user_guide
    # Command with options taken from https://github.com/weizhongli/cdhit/tree/master/psi-cd-hit
    command_62 = [cd_hit,
                    '-i', in_fname,
                    '-d', str(0),
                    '-o', out_fname,
                    '-c', str(0.62),
                    '-n', str(4),
                    '-G', str(1),
                    '-g', str(1),
                    '-b', str(20),
                    '-aL', str(0.0),
                    '-aS', str(0.0),
                    '-s', str(0.0)
                    ]
    
    log = os.path.join(wdir, "cdhit_62.log")
    p = ample_util.run_command(command_62, logfile=log, directory=wdir)

    return out_fname + '.clstr'
##End cluster()

def main(args):
    cd_hit  = os.path.abspath(args['cd_hit'])
    alnfile = os.path.abspath(args['alnfile'])
    wdir    = os.path.abspath(args['wdir'])
   
    # Removing gaps in MSA - CDHIT does not like them
    alnfile_formatted = alnfile.rsplit(".", 1)[0] + ".formatted.fasta"
    ap = parse_alignment.AlignmentParser()
    ap.resetA3M(alnfile, alnfile_formatted)
    
    # Cluster the sequences in the formatted alignment
    cluster_fname = cluster(cd_hit, alnfile_formatted, wdir)
    
    # Analyse the cluster output file
    cp = parse_cdhit.CDhitLogParser()
    cp.parse(cluster_fname)
    return cp.nrEffSeqs

if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument('cd_hit', type=str, help='CD-HIT executable')
    p.add_argument('alnfile', type=str, help='File containing Multiple Sequence Alignment')
    p.add_argument('-wdir', type=str, default=os.getcwd(), help='Working directory')
    args=vars(p.parse_args())
    neff = main(args)

    print neff
