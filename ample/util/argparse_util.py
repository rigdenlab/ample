"""
Below is code for creating a boolean argument parser with a default, and a range of choices

class BoolAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values is None: values = self.default
        if values in [ '0', 'f', 'F', 'false', 'False', False ]:
            values=False
        elif values in [ '1', 't', 'T', 'true', 'True', True ]:
            values=True
        else:
            msg = 'Unrecognised True/False value: {0}'.format(values)
            raise argparse.ArgumentError(self, msg)
        setattr(namespace, self.dest, values)

parser = argparse.ArgumentParser()

parser.add_argument('-d', '--debug', action=BoolAction, default=True, metavar='True/False', nargs='?', help="debug")
# print(parser.parse_args('-d 2'.split()))
print(parser.parse_args('-d false'.split()))

"""

__author__ = "Jens Thomas & Felix Simkovic"
__date__ = "03 Apr 2016"
__version__ = "1.0"

import argparse


def add_core_options(parser=None):
    """Function to add any arguments required by all runtypes"""
    if parser is None:
        parser = argparse.ArgumentParser()
    parser.add_argument('-config_file', help="user configuration file")
    parser.add_argument('-debug', metavar='True/False', help=argparse.SUPPRESS)
    parser.add_argument('-nproc', type=int, help="Number of processors [1]. For local, serial runs the jobs will be split across nproc processors. For cluster submission, this should be the number of processors on a node.")
    parser.add_argument('-work_dir', help='Path to the directory where the job will run (will be created if it doesn\'t exist)')
    return parser


def add_cluster_submit_options(parser=None):
    """Add the options for submission to a cluster queuing system"""
    if parser is None:
        parser = argparse.ArgumentParser()
    submit_group = parser.add_argument_group('Cluster queue submission options')
    submit_group.add_argument('-submit_array', metavar='True/False', help='Submit SGE jobs as array jobs')
    submit_group.add_argument('-submit_cluster', metavar='True/False', help='Submit jobs to a cluster - need to set -submit_qtype flag to specify the batch queue system.')
    submit_group.add_argument('-submit_max_array', type=int, help='The maximum number of jobs to run concurrently with SGE array job submission')
    submit_group.add_argument('-submit_num_array_jobs', type=int, help='The number of jobs to run concurrently with SGE array job submission')
    submit_group.add_argument('-submit_pe_lsf', help='Cluster submission: string to set number of processors for LSF queueing system')
    submit_group.add_argument('-submit_pe_sge', help='Cluster submission: string to set number of processors for SGE queueing system')
    submit_group.add_argument('-submit_queue', help='The queue to submit to on the cluster.')
    submit_group.add_argument('-submit_qtype', help='Cluster submission queue type - currently support SGE and LSF')
    return parser


def add_general_options(parser=None):
    from ample.util import version
    if parser is None:
        parser = argparse.ArgumentParser()
    add_core_options(parser)
    parser.add_argument('-alignment_file', help='Alignment file in fasta format. For homologues the first line of each sequence must be the pdb file name')
    parser.add_argument('-allow_his_tag', metavar='True/False', help='Allow HIS tags in the input sequence')
    parser.add_argument('-blast_dir', help='Directory where ncbi blast is installed (binaries in expected in bin subdirectory)')
    parser.add_argument('-classic_mode', metavar='True/False', help='Preset options to run the original AMPLE clustering/truncation options (1 cluster, 3 subclustering radii, 3 sidechains)')
    parser.add_argument('-ccp4i2_xml', help='Path to CCP4I2 XML file - if not None indicates we are running under CCP4I2')
    parser.add_argument('-devel_mode', metavar='devel_mode', help='Preset options to run in development mode - takes longer')
    parser.add_argument('-dry_run', metavar='True/False', help='Check if input files and supplied options are valid.')
    parser.add_argument('-early_terminate', metavar='True/False', help='Stop the run as soon as a success has been found.')
    parser.add_argument('-ensembles', help='Path to directory containing existing ensembles')
    parser.add_argument('-fasta', help='protein fasta file. (required)')
    parser.add_argument('-fast_protein_cluster_exe', help='path to fast_protein_cluster executable')
    parser.add_argument('-F', metavar='flag for F', help='Flag for F column in the MTZ file')
    parser.add_argument('-FREE', metavar='flag for FREE', help='Flag for FREE column in the MTZ file')
    parser.add_argument('-ideal_helices', metavar='True/False', help='Use ideal polyalanine helices to solve structure (8 helices: from 5-40 residues)')
    parser.add_argument('-improve_template', metavar='improve_template', help='Path to a template to improve - NMR, homolog')
    parser.add_argument('-LGA', metavar='path_to_LGA dir', help=argparse.SUPPRESS)
    parser.add_argument('-make_models', metavar='True/False', help= 'run rosetta modeling, set to False to import pre-made models (required if making models locally default True)')
    parser.add_argument('-max_array_jobs', help='Maximum number of array jobs to run')
    parser.add_argument('-missing_domain', metavar='True/False', help='Modelling a missing domain - requires domain_all_chains_pdb argument')
    parser.add_argument('-models', metavar='models', help='Path to a folder of PDB decoys, or a tarred and gzipped/bziped, or zipped collection of decoys')
    parser.add_argument('-mr_sequence', help="sequence file for crystal content (if different from what's given by -fasta)")
    parser.add_argument('-mtz', metavar='MTZ in', help='The MTZ file with the reflection data.')
    parser.add_argument('-name', metavar='job_name', help='4-letter identifier for job [ampl]')
    parser.add_argument('-native_pdb', metavar='native_pdb', help='Path to the crystal structure PDB for benchmarking.')
    parser.add_argument('-native_mtz', metavar='native_pdb', help='Path to the native MTZ containing FC and PHIC calculated phases for benchmarking.')
    parser.add_argument('-nmr_model_in', metavar='nmr_model_in', help='PDB with NMR models')
    parser.add_argument('-nmr_process', type=int, help='number of times to process the NMR models')
    parser.add_argument('-nmr_remodel', metavar='True/False', help='Remodel the NMR structures')
    parser.add_argument('-nmr_remodel_fasta', help='The FASTA sequence to be used for remodelling the NMR ensemble if different from the default FASTA sequence')
    parser.add_argument('-purge', metavar='purge_level', type=int, choices=[0, 1, 2], help='Delete intermediate files and failed MRBUMP results: 0 - None, 1 - Some, 2 - All possible')
    parser.add_argument('-psipred_ss2', metavar='PSIPRED_FILE', help='Psipred secondary structure prediction file')
    parser.add_argument('-quick_mode', metavar='True/False', help='Preset options to run quickly, but less thoroughly')
    parser.add_argument('-restart_pkl', help='Rerun a job using the pickled ample dictionary')
    parser.add_argument('-run_dir', metavar='run_directory', help='Directory where the AMPLE work directory will be created [current dir]')
    parser.add_argument('-rvapi_document', help='Path to an existing rvapi document (for running under jscofe)')
    parser.add_argument('-scwrl_exe', metavar='path to scwrl', help='Path to Scwrl4 executable')
    parser.add_argument('-show_gui', metavar='True/False', help='Pop up and display a stand-alone GUI')
    parser.add_argument('-single_model', help='Single structure model to be used to create ensembles')
    parser.add_argument('-sf_cif', help='Path to a structure factor CIF file (instead of MTZ file)')
    parser.add_argument('-SIGF', help='Flag for SIGF column in the MTZ file')
    parser.add_argument('-top_model_only', metavar='True/False', help='Only process the top model in each ensemble')
    parser.add_argument('--version', action='version', version='%(prog)s {0}'.format(version.__version__))
    parser.add_argument('-webserver_uri', help='URI of the webserver directory - also indicates we are running as a webserver')
    return parser


def add_contact_options(parser=None):
    """Contact prediction related options"""
    if parser is None:
        parser = argparse.ArgumentParser()
    contact_group = parser.add_argument_group("Contact Restraints Options")
    contact_group.add_argument('-bbcontacts_file', help='Additional bbcontacts file. Requires normal contactfile')
    contact_group.add_argument('-bbcontacts_format',
        help='Residue contact file format. For available formats refer to the AMPLE documentation')
    contact_group.add_argument('-contact_file', help='Residue contact file')
    contact_group.add_argument('-contact_format', help='Residue contact file format. For available formats refer to the AMPLE documentation')
    contact_group.add_argument('-disulfide_constraints_file', help='Disulfide residue constraints for ab initio modelling')
    contact_group.add_argument('-distance_to_neighbour', type=int, help="Min. distance between residue pairs for contact (default=5)")
    contact_group.add_argument('-energy_function', help='Rosetta energy function for contact restraint conversion (default=FADE)')
    contact_group.add_argument('-native_cutoff', type=float, help='Distance cutoff for reference contacts in native structure (default=8A)')
    contact_group.add_argument('--no-contact-prediction', default=False, action="store_true", help="Do not predict contacts")
    contact_group.add_argument('-restraints_factor', type=float, help='Factor (* Sequence length) determining number of contact restraints to use (default=1.0)')
    contact_group.add_argument('-restraints_file', help='Residue restraints for ab initio modelling')
    contact_group.add_argument('-restraints_weight', type=float, help="Additional energy weighting of restraints in Rosetta")
    contact_group.add_argument('-subselect_mode', help="Long-range decoy satisfaction subselection mode - one of [{0}]".format(" | ".join(["linear", "scaled", "cutoff"])))
    return parser


def add_mr_options(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser()
    mr_group = parser.add_argument_group('MRBUMP/Molecular Replacement Options')
    mr_group.add_argument('-arpwarp_cycles', type=int, help='The number of ArpWarp cycles to run')
    mr_group.add_argument('-buccaneer_cycles', type=int, help='The number of Bucanner rebuilding cycles to run')
    mr_group.add_argument('-do_mr', metavar='True/False', help='Run or skip the Molecular Replacement step')
    mr_group.add_argument('-domain_all_chains_pdb', help='Fixed input to mr bump')
    mr_group.add_argument('-domain_termini_distance', help='distance between termini for insert domains')
    mr_group.add_argument('-early_terminate_SHELXE_CC', type=float, help='SHELXE_CC criteria for when a job has succeeeded')
    mr_group.add_argument('-early_terminate_SHELXE_ACL', type=int, help='SHELXE_ACL criteria for when a job has succeeeded')
    mr_group.add_argument('-molrep_only', metavar='True/False', help='Only use Molrep for Molecular Replacement step in MRBUMP')
    mr_group.add_argument('-mrbump_dir', help='Path to a directory of MRBUMP jobs (see restart_pkl)')
    mr_group.add_argument('-mr_keys', nargs='+', action='append', help='Additional keywords for MRBUMP - are passed through without editing')
    mr_group.add_argument('-mr_sg_all', metavar='True/False', help='Try all possible space groups in PHASER Molecular Replacement step in MRBUMP')
    mr_group.add_argument('-nmasu', type=int, help='Manually specify the number of molecules in the asymmetric unit - sets the NMASu MRBUMP flag')
    mr_group.add_argument('-phaser_kill', metavar='phaser_kill', type=int, help='Time in minutes after which phaser will be killed (0 to leave running)')
    mr_group.add_argument('-phaser_only', metavar='True/False', help='Only use Phaser for Molecular Replacement step in MRBUMP')
    mr_group.add_argument('-phaser_rms', metavar='phaser_rms', help='RMS value for phaser')
    mr_group.add_argument('-refine_rebuild_arpwarp', metavar='True/False', help='True to use ARPWARP to rebuild the REFMAC-refined MR result.')
    mr_group.add_argument('-refine_rebuild_buccaneer', metavar='True/False', help='True to use Buccaneer to rebuild the REFMAC-refined MR result.')
    mr_group.add_argument('-shelx_cycles', help='The number of shelx cycles to run when rebuilding.')
    mr_group.add_argument('-shelxe_exe', metavar='path to shelxe executable', help='Path to the shelxe executable')
    mr_group.add_argument('-shelxe_rebuild', metavar='True/False', help='Rebuild shelxe traced pdb with buccaneer and arpwarp')
    mr_group.add_argument('-shelxe_rebuild_arpwarp', metavar='True/False', help='Rebuild shelxe traced pdb with arpwarp')
    mr_group.add_argument('-shelxe_rebuild_buccaneer', metavar='True/False', help='Rebuild shelxe traced pdb with buccaneer')
    mr_group.add_argument('-use_scwrl', metavar='True/False', help='Remodel sidechains of the decoy models using Scwrl4')
    mr_group.add_argument('-use_shelxe', metavar='True/False', help='True to use shelxe')
    return parser


def add_rosetta_options(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser()
    rosetta_group = parser.add_argument_group('ROSETTA Modelling Options')
    rosetta_group.add_argument('-all_atom', metavar='True/False', help="Do all-atom Rosetta modelling (adds \"-return_full_atom true\" to rosetta arguments")
    rosetta_group.add_argument('-frags_3mers', help='Path to file with pre-existing Rosetta 3mer fragments')
    rosetta_group.add_argument('-frags_9mers', help='Path to file with pre-existing Rosetta 3mer fragments')
    rosetta_group.add_argument('-make_frags', metavar='True/False', help='set True to generate Rosetta 3mers and 9mers locally, False to import fragments')
    rosetta_group.add_argument('-nmodels', default=1000, metavar='number of models', type=int, help='number of models to make (default: 1000)')
    rosetta_group.add_argument('-nr', metavar='nr', help='Path to the NR non-redundant sequence database')
    rosetta_group.add_argument('-rg_reweight', metavar='radius of gyration reweight', type=float, help='Set the Rosetta -rg_reweight flag to specify the radius of gyration reweight.')
    rosetta_group.add_argument('-rosetta_AbinitioRelax', help='Path to Rosetta AbinitioRelax executable')
    rosetta_group.add_argument('-ROSETTA_cluster', help='location of rosetta cluster')
    rosetta_group.add_argument('-rosetta_db', help='Path to the Rosetta database directory')
    rosetta_group.add_argument('-rosetta_dir', help='The Rosetta install directory')
    rosetta_group.add_argument('-rosetta_fragments_exe', help='Location of the Rosetta make_fragments.pl script')
    rosetta_group.add_argument('-rosetta_version', type=float, help='The version number of Rosetta')
    rosetta_group.add_argument('-transmembrane', metavar='True/False', help='Do Rosetta modelling for transmembrane proteins (Ovchinnikov protocol)')
    rosetta_group.add_argument('-transmembrane_old', metavar='True/False', help='Do Rosetta modelling for transmembrane proteins (Yarov-Yarovoy protocol)')
    rosetta_group.add_argument('-transmembrane_octopusfile', help='Octopus transmembrane topology predicition file')
    rosetta_group.add_argument('-transmembrane_spanfile', help='Span file for modelling transmembrane proteins')
    rosetta_group.add_argument('-transmembrane_lipofile', help='Lips4 file for modelling transmembrane proteins')
    rosetta_group.add_argument('-use_homs', metavar='True/False', help="Select ROSETTA fragments from homologous models")
    return parser


def add_ensembler_options(parser=None):
    # --------------------------------------------------------------------------------------------- #
    # sphinx-argparse ignores Mock imports and thus cannot find iotbx.pdb when generating the docs. #
    try:
        from ample.ensembler.constants import ALLOWED_SIDE_CHAIN_TREATMENTS
        from ample.ensembler.truncation_util import TRUNCATION_METHODS
    except ImportError:
        allowed_side_chain_treatments = ['polyala', 'reliable', 'allatom', 'unmod']
        truncation_methods = ['percent']
    else:
        allowed_side_chain_treatments = ALLOWED_SIDE_CHAIN_TREATMENTS[:]
        truncation_methods = [t.value for t in TRUNCATION_METHODS]
    # --------------------------------------------------------------------------------------------- #
    if parser is None:
        import argparse
        parser = argparse.ArgumentParser()
    ensembler_group = parser.add_argument_group('Ensemble Options')
    ensembler_group.add_argument('-cluster_dir', help='Path to directory of pre-clustered models to import')
    ensembler_group.add_argument('-cluster_method', help='How to cluster the models for ensembling (spicker|fast_protein_cluster')
    ensembler_group.add_argument('-ensembler_timeout', type=int, help='Time in seconds before timing out ensembling')
    ensembler_group.add_argument('-gesamt_exe', metavar='gesamt_exe', help='Path to the gesamt executable')
    ensembler_group.add_argument('-homologs', metavar='True/False', help='Generate ensembles from homologs models (requires -alignment_file)')
    ensembler_group.add_argument('-homolog_aligner', metavar='homolog_aligner', help='Program to use for structural alignment of homologs (gesamt|mustang)')
    ensembler_group.add_argument('-ensemble_max_models', help='Maximum number of models permitted in an ensemble')
    ensembler_group.add_argument('-maxcluster_exe', help='Path to Maxcluster executable')
    ensembler_group.add_argument('-mustang_exe', metavar='mustang_exe', help='Path to the mustang executable')
    ensembler_group.add_argument('-num_clusters', type=int, help='The number of Spicker clusters of the original decoys that will be sampled [1]')
    ensembler_group.add_argument('-percent', metavar='percent_truncation', help='percent interval for truncation')
    ensembler_group.add_argument('-percent_fixed_intervals', nargs='+', type=int, help='list of integer percentage intervals for truncation')
    ensembler_group.add_argument('-score_matrix', help='Path to score matrix for spicker')
    ensembler_group.add_argument('-score_matrix_file_list', help='File with list of ordered model names for the score_matrix')
    ensembler_group.add_argument('-side_chain_treatments', type=str, nargs='+', help='The side chain treatments to use. Options: ' + '|'.join(allowed_side_chain_treatments))
    ensembler_group.add_argument('-spicker_exe', help='Path to spicker executable')
    ensembler_group.add_argument('-subcluster_radius_thresholds', type=float, nargs='+', help='The radii to use for subclustering the truncated ensembles')
    ensembler_group.add_argument('-subcluster_program', help='Program for subclustering models [maxcluster]')
    ensembler_group.add_argument('-theseus_exe', metavar='Theseus exe (required)', help='Path to theseus executable')
    ensembler_group.add_argument('-thin_clusters', metavar='True/False', help='Create ensembles from 10 clusters with 1 + 3A subclustering and polyAlanine sidechains')
    ensembler_group.add_argument('-truncation_method', help='How to truncate the models for ensembling: ' + '|'.join(truncation_methods))
    ensembler_group.add_argument('-truncation_pruning', help='Whether to remove isolated residues (single)')
    ensembler_group.add_argument('-truncation_scorefile', help="CSV file containing per residue scores - COLUMN ONE MUST BE RESIDUE INDEX STARTING FROM 1")
    ensembler_group.add_argument('-truncation_scorefile_header', nargs='+', help="column headers to be used to create ensembles")
    return parser

def process_command_line(args=None, contacts=True, modelling=True, mol_rep=True):
    """Process the command-line for the main AMPLE program.
    :args: optional argument that can hold the command-line arguments if we
    have been called from within python for testing
    """
    parser = argparse.ArgumentParser(
        description="AMPLE: Ab initio Modelling of Proteins for moLEcular replacement", prefix_chars="-")
    add_general_options(parser)
    add_cluster_submit_options(parser)
    add_ensembler_options(parser)

    if contacts:
        add_contact_options(parser)
    if mol_rep:
        add_mr_options(parser)
    if modelling:
        add_rosetta_options(parser)

    return parser.parse_args(args)
