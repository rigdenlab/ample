#!/usr/bin/env ccp4-python

"""
11.11.2015

@author: hlfsimko
"""

import argparse
import logging
import os
import sys


if not "CCP4" in os.environ.keys(): raise RuntimeError('CCP4 not found')
sys.path.insert(0, os.path.join(os.environ['CCP4'], "share", "ample", "parsers"))
#sys.path.append(os.path.join(os.environ["HOME"], "opt", "ample-dev1", "parsers"))


import parse_bbcontacts
import parse_casprr
import parse_ccmpred
import parse_epcmap
import parse_evfold
import parse_gremlin
import parse_pconsc


def main(args):
    contactfile = os.path.abspath(args['contactfile'])

    outfile = os.path.abspath(args['outfile']) \
                    if args['outfile'] \
                    else os.path.abspath(os.path.join(os.getcwd(), 
                                         os.path.basename(contactfile).rsplit('.', 1)[0] + ".CASPRR"))

    if args['bbcontacts']:
        cp = parse_bbcontacts.BBcontactsContactParser()
        cp.read(contactfile)
    elif args['ccmpred']:
        cp = parse_ccmpred.CCMpredContactParser()
        cp.read(contactfile)
    elif args['epcmap']:
        cp = parse_epcmap.EPCMapContactParser()
        cp.read(contactfile)
    elif args['evfold']:
        cp = parse_evfold.EVfoldContactParser()
        cp.read(contactfile)
    elif args['gremlin']:
        cp = parse_gremlin.GremlinContactParser()
        cp.read(contactfile)
    elif args['pconsc']:
        cp = parse_pconsc.PconscContactParser()
        cp.read(contactfile)
        
    contacts = cp.getContacts()

    op = parse_casprr.CaspContactParser()
    op.setContacts(contacts)
    op.sortContacts('res1_index', descending=False)
    op.write(outfile)
##End main()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-o', metavar="FILE", dest="outfile", type=str, 
                        help="output filename")
    parser.add_argument('contactfile', type=str, help="contact filename")

    contacts = parser.add_mutually_exclusive_group(required=True)
    contacts.add_argument("-bb", action="store_true", dest="bbcontacts",
                          help="bbcontacts contact file")
    contacts.add_argument("-ccm", action="store_true", dest="ccmpred",
                          help="CCMpred contact matrix")
    contacts.add_argument("-epc", action="store_true", dest="epcmap",
                          help="EPC-Map contact file")
    contacts.add_argument("-evf", action="store_true", dest="evfold",
                          help="EVFold contact file")
    contacts.add_argument("-gre", action="store_true", dest="gremlin",
                          help="GREMLIN contact file")
    contacts.add_argument("-pco", action="store_true", dest="pconsc",
                          help="PconsC1/2/3 contact file")

    args = vars(parser.parse_args())
    main(args)
