#!/usr/local/bin/python2.7
# encoding: utf-8
'''

Can pass in both ample dictionary or command-line args

Parse main cli with ensemble-cli with:

ensemble_group = parser.add_argument_group("Ensembler Options")

# Have list of arguments each with
ensemble_opts = OrderedDict(
help_opts = {
'short_opt' : '-h',
'long_opt' : '--help_opts', # Must be same as key with two underscores appended
'dest' : None,
'type' : 'str',
'nargs' : 1,
'help' : 'Additional bbcontacts file. Requires normal contactfile'
},
)

for k,v in amopt.d:
    if k in ensemble_opts:
        optd = ensemble_opts[k]
        # Need to get first two options and then **kw the rest
        parser.add_argument(
'''

from ample.ensembler.constants import SIDE_CHAIN_TREATMENTS
from ample.util.argparse_util import add_general_options

def add_ensembler_options(parser, standalone=False):
    """Function to add the ensemble-specific options"""
    
    if standalone:
        # Add options for running as a standalone module
        add_general_options(parser)
    else:
        parser = parser.add_argument_group('Ensemble Options')

    # Ensenble-specific options here
    parser.add_argument('-cluster_dir', type=str, nargs=1,
                       help='Path to directory of pre-clustered models to import')
    
    parser.add_argument('-cluster_method', type=str, nargs=1,
                       help='How to cluster the models for ensembling (spicker|fast_protein_cluster')

    parser.add_argument('-ensembler_timeout', type=int, nargs=1,
                       help='Time in seconds before timing out ensembling')

    parser.add_argument('-gesamt_exe', metavar='gesamt_exe', type=str, nargs=1,
                       help='Path to the gesamt executable')

    parser.add_argument('-homologs', metavar='True/False', type=str, nargs=1,
                       help='Generate ensembles from homologs models (requires -alignment_file)')
    
    parser.add_argument('-homolog_aligner', metavar='homolog_aligner', type=str, nargs=1,
                       help='Program to use for structural alignment of homologs (gesamt|mustang)')

    parser.add_argument('-max_ensemble_models', type=str, nargs=1,
                       help='Maximum number of models permitted in an ensemble')

    parser.add_argument('-maxcluster_exe', type=str, nargs=1,
                       help='Path to Maxcluster executable')
    
    parser.add_argument('-mustang_exe', metavar='mustang_exe', type=str, nargs=1,
                       help='Path to the mustang executable')
    
    parser.add_argument('-num_clusters', type=int, nargs=1,
                       help='The number of Spicker clusters of the original decoys that will be sampled [1]')

    parser.add_argument('-percent', metavar='percent_truncation', type=str, nargs=1,
                       help='percent interval for truncation')
            
    parser.add_argument('-score_matrix', type=str, nargs=1,
                       help='Path to score matrix for spicker')
    
    parser.add_argument('-score_matrix_file_list', type=str, nargs=1,
                       help='File with list of ordered model names for the score_matrix')
     
    parser.add_argument('-side_chain_treatments', type=str, nargs='+', action='append',
                       help='The side chain treatments to use. Default: {0}'.format(SIDE_CHAIN_TREATMENTS))
    
    parser.add_argument('-subcluster_program', type=str, nargs=1,
                       help='Program for subclustering models [maxcluster]')
        
    parser.add_argument('-theseus_exe', metavar='Theseus exe (required)', type=str, nargs=1,
                       help='Path to theseus executable')

    parser.add_argument('-truncation_method', type=str, nargs=1,
                       help='How to truncate the models for ensembling percent|thresh|focussed|scores')
    
    parser.add_argument('-truncation_pruning', type=str, nargs=1,
                       help='Whether to remove isolated residues (single)')
    
    parser.add_argument('-truncation_scorefile', type=str, nargs=1,
                        help="CSV file containing per residue scores - COLUMN ONE MUST BE RESIDUE INDEX STARTING FROM 1")
    
    parser.add_argument('-truncation_scorefile_header', type=str, nargs='+', action='append',
                        help="column headers to be used to create ensembles")
              
    return
