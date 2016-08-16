"""
@author: jmht, hlfsimko


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
#print parser.parse_args('-d 2'.split())
print parser.parse_args('-d false'.split())


"""

from ample.util import version

def add_core_options(parser):
    """Function to add any arguments required by all runtypes"""
    
    parser.add_argument('-config_file', type=str, help="user configuration file")
        
    parser.add_argument('-debug', metavar='True/False', type=str, nargs=1,
                       help='Run in debug mode (CURRENTLY UNUSED)')
    
    parser.add_argument('-nproc', type=int, default=1,
                       help="Number of processors [1]. For local, serial runs the jobs will be split across nproc processors." + \
                        "For cluster submission, this should be the number of processors on a node.")

    parser.add_argument('-work_dir', type=str, nargs=1,
                       help='Path to the directory where the job will run (will be created if it doesn\'t exist)')
    
    return

def add_general_options(parser):
    
    # Always add the core option
    add_core_options(parser)
    
    parser.add_argument('-alignment_file', type=str, nargs=1,
                       help='Alignment file in fasta format. For homologues the first line of each sequence must be the pdb file name')
    
    parser.add_argument('-allow_his_tag', metavar='True/False', type=str, nargs=1,
                       help='Allow HIS tags in the input sequence')
    
    parser.add_argument('-blast_dir', type=str, nargs=1,
                       help='Directory where ncbi blast is installed (binaries in expected in bin subdirectory)')
    
    parser.add_argument('-ccp4_jobid', type=int, nargs=1,
                       help='Set the CCP4 job id - only needed when running from the CCP4 GUI')
    
    parser.add_argument('-devel_mode', metavar='devel_mode', type=str, nargs=1,
                       help='Preset options to run in development mode - takes longer')
    
    parser.add_argument('-dry_run', metavar='True/False', type=str, nargs=1,
                         help='Check if input files and supplied options are valid.')
    
    parser.add_argument('-early_terminate', metavar='True/False', type=str, nargs=1,
                         help='Stop the run as soon as a success has been found.')
    
    parser.add_argument('-ensembles', type=str, nargs=1,
                       help='Path to directory containing existing ensembles')
    
    parser.add_argument('-fasta', type=str, nargs=1,
                       help='protein fasta file. (required)')
    
    parser.add_argument('-fast_protein_cluster_exe', type=str, nargs=1,
                       help='path to fast_protein_cluster executable')
    
    parser.add_argument('-F', metavar='flag for F', type=str, nargs=1,
                       help='Flag for F column in the MTZ file')
    
    parser.add_argument('-FREE', metavar='flag for FREE', type=str, nargs=1,
                       help='Flag for FREE column in the MTZ file')
            
    parser.add_argument('-ideal_helices', metavar='True/False', type=str, nargs=1,
                       help='Use ideal polyalanine helices to solve structure (8 helices: from 5-40 residues)')

    parser.add_argument('-improve_template', metavar='improve_template', type=str, nargs=1,
                       help='Path to a template to improve - NMR, homolog')
    
    parser.add_argument('-LGA', metavar='path_to_LGA dir', type=str, nargs=1,
                       help='pathway to LGA folder (not the exe) will use the \'lga\' executable. UNUSED')

    parser.add_argument('-make_models', metavar='True/False', type=str, nargs=1,
                       help='run rosetta modeling, set to False to import pre-made models (required if making models locally default True)')
    
    parser.add_argument('-max_array_jobs', type=str, nargs=1,
                       help='Maximum number of array jobs to run')
    
    parser.add_argument('-missing_domain', metavar='True/False', type=str, nargs=1,
                       help='Modelling a missing domain - requires domain_all_chains_pdb argument')
    
    parser.add_argument('-models', metavar='models', type=str, nargs=1,
                       help='Path to a folder of PDB decoys, or a tarred and gzipped/bziped, or zipped collection of decoys')

    parser.add_argument('-mr_sequence', type=str, nargs=1,
                       help="sequence file for crystal content (if different from what's given by -fasta)")

    parser.add_argument('-mtz', metavar='MTZ in', type=str, nargs=1,
                       help='The MTZ file with the reflection data.')

    parser.add_argument('-name', metavar='job_name', type=str, nargs=1,
                       help='4-letter identifier for job [ampl]')
    
    parser.add_argument('-native_pdb', metavar='native_pdb', type=str, nargs=1,
                       help='Path to the crystal structure PDB for benchmarking.')
    
    parser.add_argument('-native_mtz', metavar='native_pdb', type=str, nargs=1,
                       help='Path to the native MTZ containing FC and PHIC calculated phases for benchmarking.')
    
    parser.add_argument('-nmodels', metavar='number of models', type=int, nargs=1,
                       help='number of models to make (default: 1000)')
    
    parser.add_argument('-nr', metavar='nr', type=str, nargs=1,
                       help='Path to the NR non-redundant sequence database')
    
    parser.add_argument('-nmr_model_in', metavar='nmr_model_in', type=str, nargs=1,
                       help='PDB with NMR models')
    
    parser.add_argument('-nmr_process', type=int, nargs=1,
                       help='number of times to process the NMR models')
    
    parser.add_argument('-nmr_remodel', metavar='True/False', type=str, nargs=1,
                       help='Remodel the NMR structures')
    
    parser.add_argument('-nmr_remodel_fasta', type=str, nargs=1,
                       help='The FASTA sequence to be used for remodelling the NMR ensemble if different from the default FASTA sequence')
    
    parser.add_argument('-no_gui', metavar='True/False', type=str, nargs=1,
                       help='Do not display the AMPLE gui.')
    
    parser.add_argument('-output_pdb', type=str, nargs=1,
                       help='Name of the final result pdb to output [ample_output.pdb]')
    
    parser.add_argument('-purge', metavar='True/False', type=str, nargs=1,
                       help='Delete all intermediate files and failed MRBUMP results')

    parser.add_argument('-psipred_ss2', metavar='PSIPRED_FILE', type=str, nargs=1,
                       help='Psipred secondary structure prediction file')
    
    parser.add_argument('-quick_mode', metavar='True/False', type=str, nargs=1,
                       help='Preset options to run quickly, but less thoroughly')
    
    parser.add_argument('-restart_pkl', type=str, help='Rerun a job using the pickled ample dictionary')
    
    parser.add_argument('-run_dir', metavar='run_directory', type=str, nargs=1,
                       help='Directory where the AMPLE work directory will be created [current dir]')
    
    parser.add_argument('-scwrl_exe', metavar='path to scwrl', type=str, nargs=1,
                       help='Path to Scwrl4 executable')
    
    parser.add_argument('-single_model', type=str, nargs=1,
                       help='Single structure model to be used to create ensembles')
    
    parser.add_argument('-sf_cif', type=str, nargs=1,
                       help='Path to a structure factor CIF file (instead of MTZ file)')
    
    parser.add_argument('-SIGF', type=str, nargs=1,
                       help='Flag for SIGF column in the MTZ file')
    
    parser.add_argument('-spicker_exe', type=str, nargs=1,
                       help='Path to spicker executable')
    
    parser.add_argument('-submit_array', metavar='True/False', type=str, nargs=1,
                       help='Submit SGE jobs as array jobs')
    
    parser.add_argument('-submit_cluster', metavar='True/False', type=str, nargs=1,
                       help='Submit jobs to a cluster - need to set -submit_qtype flag to specify the batch queue system.')
    
    parser.add_argument('-submit_qtype', type=str, nargs=1,
                       help='Cluster submission queue type - currently support SGE and LSF')

    parser.add_argument('-submit_pe_lsf', type=str, nargs=1,
                       help='Cluster submission: string to set number of processors for LSF queueing system')

    parser.add_argument('-submit_pe_sge', type=str, nargs=1,
                       help='Cluster submission: string to set number of processors for SGE queueing system')
    
    parser.add_argument('-submit_queue', type=str, nargs=1,
                       help='The queue to submit to on the cluster.')
    
    parser.add_argument('-top_model_only', metavar='True/False', type=str, nargs=1,
                       help='Only process the top model in each ensemble')
    
    parser.add_argument('--version', action='version', version='%(prog)s {0}'.format(version.__version__))
    
    parser.add_argument('-webserver_uri', type=str, nargs=1,
                       help='URI of the webserver directory - also indicates we are running as a webserver')
    return

def add_contact_options(parser):
    
    contact_group = parser.add_argument_group("Contact Restraints Options")
    
    contact_group.add_argument('-bbcontacts_file', type=str, nargs=1,
                       help='Additional bbcontacts file. Requires normal contactfile')
    
    contact_group.add_argument('-contact_file', type=str, nargs=1,
                       help='Residue contact file in CASP RR format')
    
    contact_group.add_argument('-disulfide_constraints_file', type=str, nargs=1,
                       help='Disulfide residue constraints for ab initio modelling')
    
    contact_group.add_argument('-distance_to_neighbour', type=int, nargs=1,
                       help="Min. distance between residue pairs for contact (default=5)")
    
    contact_group.add_argument('-energy_function', type=str, nargs=1,
                       help='Rosetta energy function for contact restraint conversion (default=FADE)')
    
    contact_group.add_argument('-native_cutoff', type=float, nargs=1,
                       help='Distance cutoff for reference contacts in native structure (default=8A)')
    
    contact_group.add_argument('-restraints_factor', type=float, nargs=1,
                       help='Factor (* Sequence length) determining number of contact restraints to use (default=1.0)')
    
    contact_group.add_argument('-restraints_file', type=str, nargs=1,
                       help='Residue restraints for ab initio modelling')
    
    contact_group.add_argument('-restraints_weight', type=float, nargs=1,
                       help="Additional energy weighting of restraints in Rosetta")
    return

def add_mr_options(parser): 
    
    mr_group = parser.add_argument_group('MRBUMP/Molecular Replacement Options')
    
    mr_group.add_argument('-arpwarp_cycles', type=int, nargs=1,
                       help='The number of ArpWarp cycles to run') 
    
    mr_group.add_argument('-buccaneer_cycles', type=int, nargs=1,
                       help='The number of Bucanner rebuilding cycles to run')

    mr_group.add_argument('-do_mr', type=str, metavar='True/False', nargs=1,
                       help='Run or skip the Molecular Replacement step')

    mr_group.add_argument('-domain_all_chains_pdb', type=str, nargs=1,
                       help='Fixed input to mr bump')
    
    mr_group.add_argument('-domain_termini_distance', type=str, nargs=1,
                       help='distance between termini for insert domains')

    mr_group.add_argument('-molrep_only', metavar='True/False', type=str, nargs=1,
                       help='Only use Molrep for Molecular Replacement step in MRBUMP')
    
    mr_group.add_argument('-mrbump_dir', type=str, nargs=1,
                       help='Path to a directory of MRBUMP jobs (see restart_pkl)')
    
    mr_group.add_argument('-mr_keys', type=str, nargs='+', action='append',
                       help='Additional keywords for MRBUMP - are passed through without editing')
    
    mr_group.add_argument('-mr_sg_all', metavar='True/False', type=str, nargs=1,
                       help='Try all possible space groups in PHASER Molecular Replacement step in MRBUMP')

    mr_group.add_argument('-nmasu', type=int, nargs=1,
                       help='Manually specify the number of molecules in the asymmetric unit - sets the NMASu MRBUMP flag')
    
    mr_group.add_argument('-phaser_kill', metavar='phaser_kill', type=int, nargs=1,
                       help='Time in minutes after which phaser will be killed (0 to leave running)')

    mr_group.add_argument('-phaser_only', metavar='True/False', type=str, nargs=1,
                       help='Only use Phaser for Molecular Replacement step in MRBUMP')

    mr_group.add_argument('-phaser_rms', metavar='phaser_rms', type=str, nargs=1,
                       help='RMS value for phaser')

    mr_group.add_argument('-shelx_cycles', type=str, nargs=1,
                         help='The number of shelx cycles to run when rebuilding.')
    
    mr_group.add_argument('-shelxe_exe', metavar='path to shelxe executable', type=str, nargs=1,
                       help='Path to the shelxe executable')
    
    mr_group.add_argument('-shelxe_rebuild', metavar='True/False', type=str, nargs=1,
                       help='Rebuild shelxe traced pdb with buccaneer and arpwarp')
    
    mr_group.add_argument('-shelxe_rebuild_arpwarp', metavar='True/False', type=str, nargs=1,
                       help='Rebuild shelxe traced pdb with arpwarp')
    
    mr_group.add_argument('-shelxe_rebuild_buccaneer', metavar='True/False', type=str, nargs=1,
                       help='Rebuild shelxe traced pdb with buccaneer')

    mr_group.add_argument('-use_arpwarp', metavar='True/False', type=str, nargs=1,
                       help='True to use arpwarp to rebuild.')
    
    mr_group.add_argument('-use_buccaneer', metavar='True/False', type=str, nargs=1,
                       help='True to use Buccaneer')
    
    mr_group.add_argument('-use_scwrl', metavar='True/False', type=str, nargs=1,
                       help='Remodel sidechains of the decoy models using Scwrl4')
    
    mr_group.add_argument('-use_shelxe', metavar='True/False', type=str, nargs=1,
                       help='True to use shelxe')
    return

def add_rosetta_options(parser):
    
    rosetta_group = parser.add_argument_group('ROSETTA Modelling Options')

    rosetta_group.add_argument('-all_atom', metavar='True/False', type=str, nargs=1,
                       help="Do all-atom Rosetta modelling (adds \"-return_full_atom true\" to rosetta arguments")
    
    rosetta_group.add_argument('-frags_3mers', type=str, nargs=1,
                       help='Path to file with pre-existing Rosetta 3mer fragments')
    
    rosetta_group.add_argument('-frags_9mers', type=str, nargs=1,
                       help='Path to file with pre-existing Rosetta 3mer fragments')
    
    rosetta_group.add_argument('-make_frags', metavar='True/False', type=str, nargs=1,
                       help='set True to generate Rosetta 3mers and 9mers locally, False to import fragments')        
    
    rosetta_group.add_argument('-rg_reweight', metavar='radius of gyration reweight', type=float, nargs=1,
                       help='Set the Rosetta -rg_reweight flag to specify the radius of gyration reweight.')
    
    rosetta_group.add_argument('-rosetta_AbinitioRelax', type=str, nargs=1,
                       help='Path to Rosetta AbinitioRelax executable')
    
    rosetta_group.add_argument('-ROSETTA_cluster', type=str, nargs=1,
                       help='location of rosetta cluster')
    
    rosetta_group.add_argument('-rosetta_db', type=str, nargs=1,
                       help='Path to the Rosetta database directory')
    
    rosetta_group.add_argument('-rosetta_dir', type=str, nargs=1,
                       help='The Rosetta install directory')

    rosetta_group.add_argument('-rosetta_fragments_exe', type=str, nargs=1,
                       help='Location of the Rosetta make_fragments.pl script')
    
    rosetta_group.add_argument('-rosetta_version', type=float, nargs=1,
                       help='The version number of Rosetta')

    rosetta_group.add_argument('-transmembrane', type=str, nargs=1,
                       help='Do Rosetta modelling for transmembrane proteins')
    
    rosetta_group.add_argument('-transmembrane2', type=str, nargs=1,
                       help='Do Rosetta modelling for transmembrane proteins (NEW PROTOCOL)')
    
    rosetta_group.add_argument('-transmembrane_octopusfile', type=str, nargs=1,
                       help='Octopus transmembrane topology predicition file')
    
    rosetta_group.add_argument('-transmembrane_spanfile', type=str, nargs=1,
                       help='Span file for modelling transmembrane proteins')
    
    rosetta_group.add_argument('-transmembrane_lipofile', type=str, nargs=1,
                       help='Lips4 file for modelling transmembrane proteins')

    rosetta_group.add_argument('-use_homs', metavar='True/False', type=str, nargs=1,
                       help='Select ROSETTA fragments from homologous models')
    return


    
