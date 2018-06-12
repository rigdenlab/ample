
__all__ = ["CLUSTER_ARGS", "EXAMPLE_DIRS", "EXTRA_ARGS"]

################################################################################
# Arguments used on the cluster
# Any set to None are ignored and only used for getting the key names in _update_cluster_args
CLUSTER_ARGS = {
                'submit_cluster' : None,
                'submit_qtype' : 'SGE',
                'submit_array' : None,
                'submit_pe_lsf' : None,
                'submit_pe_sge' : 'smp',
                'submit_queue' : None,
                #'-submit_max_array', None ],
                #'-submit_queue', None ],
}

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
    [ '-nproc',1], # Each test case needs to be run on a single processor
    [ '-do_mr','False'],
    [ '-purge','2'],
]
