#!/usr/bin/env ccp4-python

"""
11.11.2015

@author: hlfsimko
"""

import argparse
import os
import sys

from ample.parsers import bbcontacts_parser
from ample.parsers import bclcontact_parser
from ample.parsers import casprr_parser
from ample.parsers import ccmpred_parser
from ample.parsers import epcmap_parser
from ample.parsers import evfold_parser
from ample.parsers import gremlin_parser
from ample.parsers import pconsc_parser

OPTIONS = {'bbcontacts' : bbcontacts_parser.BBcontactsContactParser,
           'bclcontact' : bclcontact_parser.BCLContactsContactParser,
           'ccmpred' : ccmpred_parser.CCMpredContactParser,
           'epcmap' : epcmap_parser.EPCMapContactParser,
           'evfold' : evfold_parser.EVfoldContactParser,
           'gremlin' : gremlin_parser.GremlinContactParser,
           'pconsc' : pconsc_parser.PconscContactParser}

def main(args):
    contactfile = os.path.abspath(args['contactfile'])

    if args['outfile']:
        outfile = os.path.abspath(args['outfile'])
    else:
        fname = os.path.basename(contactfile).rsplit('.', 1)[0] + ".CASPRR"
        outfile = os.path.join(os.getcwd(), fname)
                    
    contact_parser = OPTIONS.get(args['method'])
    if not contact_parser: 
        print "No parser defined for this method: {0}".format(args['method'])
        sys.exit(1)
    
    cp = contact_parser()
    cp.read(contactfile)
     
    contacts = cp.get_contacts()
 
    op = casprr_parser.CaspContactParser()
    op.set_contacts(contacts)
    op.sort_contacts('res1_index', descending=False)
    op.write(outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-m",  metavar="OPTION", dest="method", type=str,
                        help="options are [ {0} ]".format(" | ".join(OPTIONS.keys())))
    parser.add_argument('-o', metavar="FILE", dest="outfile", type=str, 
                        help="output filename")
    parser.add_argument('contactfile', type=str, help="contact filename")
    args = vars(parser.parse_args())
    main(args)
