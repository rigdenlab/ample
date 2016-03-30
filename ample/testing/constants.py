
__all__ = ["CLUSTER_ARGS", "EXAMPLE_DIRS", "EXTRA_ARGS"]

################################################################################
# Arguments used on the cluster
CLUSTER_ARGS = [ 
    [ '-submit_cluster', 'True' ],
    [ '-submit_qtype', 'SGE' ],
    [ '-submit_array', 'True' ],
    [ '-no_gui', 'True' ],
    #[ '-submit_max_array', None ],
    #[ '-submit_queue', None ],
]

################################################################################
# List of which test directories to process
# TODO: 'transmembrane.3LBW'
#    'missing-domain.1k04',
EXAMPLE_DIRS = [
    'contact-example',
    'homologs',
    'ideal-helices',
    'import-data',
    'nmr.remodel',
    'nmr.truncate',
    'single-model',
    'toxd-example',
]

################################################################################
# Any args that are to be added/updated
EXTRA_ARGS = [ 
    ['-no_gui','True' ],
    #[ '-do_mr','False'],
]
