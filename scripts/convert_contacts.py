#!/usr/bin/env ccp4-python

"""
11.11.2015

@author: hlfsimko
"""

import argparse
import logging
import os
import sys

from ample.parsers import bbcontacts_parser
from ample.parsers import casprr_parser
from ample.parsers import ccmpred_parser
from ample.parsers import epcmap_parser
from ample.parsers import evfold_parser
from ample.parsers import gremlin_parser
from ample.parsers import pconsc_parser

def main(args):
    contactfile = os.path.abspath(args['contactfile'])

    outfile = os.path.abspath(args['outfile']) \
                    if args['outfile'] \
                    else os.path.abspath(os.path.join(os.getcwd(), 
                                         os.path.basename(contactfile).rsplit('.', 1)[0] + ".CASPRR"))

    if args['bbcontacts']:
        cp = bbcontacts_parser.BBcontactsContactParser()
        cp.read(contactfile)
    elif args['ccmpred']:
        cp = ccmpred_parser.CCMpredContactParser()
        cp.read(contactfile)
    elif args['epcmap']:
        cp = epcmap_parser.EPCMapContactParser()
        cp.read(contactfile)
    elif args['evfold']:
        cp = evfold_parser.EVfoldContactParser()
        cp.read(contactfile)
    elif args['gremlin']:
        cp = gremlin_parser.GremlinContactParser()
        cp.read(contactfile)
    elif args['pconsc']:
        cp = pconsc_parser.PconscContactParser()
        cp.read(contactfile)
        
    contacts = cp.getContacts()

    op = casprr_parser.CaspContactParser()
    op.setContacts(contacts)
    op.sortContacts('res1_index', descending=False)
    op.write(outfile)

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
