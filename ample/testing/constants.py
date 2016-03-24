
import os

__all__ = ["AMPLE_DIR", "CLUSTER_ARGS", "EXAMPLE_DIRS", "EXTRA_ARGS"]

################################################################################
# AMPLE's root directory
AMPLE_DIR = os.sep.join(os.path.abspath(os.path.dirname(__file__)).split(os.sep)[ :-1 ])
SHARE_DIR = os.path.join(os.environ["CCP4"], "share", "ample")
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
