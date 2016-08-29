.. _ample_options:

*************
AMPLE Options
*************
The following list contains all the options that AMPLE takes as arguments. It is divided into several sections depending on the part of the AMPLE pipeline.

.. _general_options:

General Options
---------------
=============================   ===========================
Command-line flag               Description
=============================   ===========================
``-config_file``                user configuration file
``-debug``                      run in debug mode (CURRENTLY UNUSED)
``-nproc``                      number of processors [1]. For local, serial runs the jobs will be split across nproc processors.For cluster submission, this should be the number of processors on a node.
``-work_dir``                   path to the directory where the job will run (will be created if it doesn't exist)
``-alignment_file``             alignment file in fasta format. For homologues the first line of each sequence must be the pdb file name
``-allow_his_tag``              allow HIS tags in the input sequence
``-blast_dir``                  directory where ncbi blast is installed (binaries in expected in bin subdirectory)
``-ccp4_jobid``                 set the CCP4 job id - only needed when running from the CCP4 GUI
``-devel_mode``                 preset options to run in development mode - takes longer
``-dry_run``                    check if input files and supplied options are valid.
``-early_terminate``            stop the run as soon as a success has been found.
``-ensembles``                  path to directory containing existing ensembles
``-fasta``                      protein fasta file. (required)
``-fast_protein_cluster_exe``   path to fast_protein_cluster executable
``-F``                          flag for F column in the MTZ file
``-FREE``                       flag for FREE column in the MTZ file
``-ideal_helices``              use ideal polyalanine helices to solve structure (8 helices: from 5-40 residues)
``-improve_template``           path to a template to improve - NMR, homolog
``-LGA path_to_LGA``            pathway to LGA folder (not the exe) will use the 'lga' executable. UNUSED
``-make_models``                run rosetta modeling, set to False to import pre-made models (required if making models locally default True)
``-max_array_jobs``             maximum number of array jobs to run
``-missing_domain``             modelling a missing domain - requires domain_all_chains_pdb argument
``-models models``              path to a folder of PDB decoys, or a tarred and gzipped/bziped, or zipped collection of decoys
``-mr_sequence``                sequence file for crystal content (if different from what's given by -fasta)
``-mtz``                        the MTZ file with the reflection data.
``-name``                       4-letter identifier for job [ampl]
``-native_pdb``                 path to the crystal structure PDB for benchmarking.
``-nmodels``                    number of models to make (default: 1000)
``-nr``                         path to the NR non-redundant sequence database
``-nmr_model_in``               PDB with NMR models
``-nmr_process``                number of times to process the NMR models
``-nmr_remodel``                remodel the NMR structures
``-nmr_remodel_fasta``          FASTA sequence to be used for remodelling the NMR ensemble if different from the default FASTA sequence
``-no_gui``                     do not display the AMPLE gui.
``-output_pdb``                 name of the final result pdb to output [ample_output.pdb]
``-purge``                      delete all intermediate files and failed MRBUMP results
``-psipred_ss2``                psipred secondary structure prediction file
``-quick_mode``                 preset options to run quickly, but less thoroughly
``-restart_pkl``                rerun a job using the pickled ample dictionary
``-run_dir``                    directory where the AMPLE work directory will be created [current dir]
``-scwrl_exe``                  path to Scwrl4 executable
``-single_model``               single structure model to be used to create ensembles
``-sf_cif``                     path to a structure factor CIF file (instead of MTZ file)
``-SIGF``                       flag for SIGF column in the MTZ file
``-spicker_exe``                path to spicker executable
``-submit_array``               submit SGE jobs as array jobs
``-submit_cluster``             submit jobs to a cluster - need to set -submit_qtype flag to specify the batch queue system.
``-submit_qtype``               cluster submission queue type - currently support SGE and LSF
``-submit_pe_lsf``              cluster submission: string to set number of processors for LSF queueing system
``-submit_pe_sge``              cluster submission: string to set number of processors for SGE queueing system
``-submit_queue``               the queue to submit to on the cluster.
``-top_model_only``             only process the top model in each ensemble
``--version``                   show program's version number and exit
``-webserver_uri``              URI of the webserver directory - also indicates we are running as a webserver
=============================   ===========================

.. _ensemble_options:

Ensemble Options
----------------
================================  ===========================
Command-line flag                 Description
================================  ===========================
``-cluster_dir``                  path to directory of pre-clustered models to import
``-cluster_method``               how to cluster the models for ensembling (spicker|fast_protein_cluster)
``-gesamt_exe``                   path to the gesamt executable
``-homologs``                     generate ensembles from homologs models (requires ``-alignment_file``)
``-homolog_aligner``              program to use for structural alignment of homologs (gesamt|mustang)
``-max_ensemble_models``          maximum number of models permitted in an ensemble
``-maxcluster_exe``               path to Maxcluster executable
``-mustang_exe``                  path to the mustang executable
``-num_clusters``                 the number of Spicker clusters of the original decoys that will be sampled [1]
``-percent``                      percent interval for truncation
``-score_matrix``                 path to score matrix for spicker
``-score_matrix_file_list``       file with list of ordered model names for the score_matrix
``-side_chain_treatments``        side chain treatments to use. Default: ['polyAla', 'reliable', 'allatom']
``-subcluster_program``           program for subclustering models [maxcluster]
``-theseus_exe``                  path to theseus executable
``-truncation_method``            how to truncate the models for ensembling [percent|thresh|focussed|scores]
``-truncation_pruning``           whether to remove isolated residues (single)
``-truncation_scorefile``         CSV file containing per residue scores - COLUMN ONE MUST BE RESIDUE INDEX STARTING FROM 1
``-truncation_scorefile_header``  column headers to be used to create ensembles
================================  ===========================

.. _contact_options:

Contact Restraints Options
--------------------------
===============================   ===========================
Command-line flag                 Description
===============================   ===========================
``-bbcontacts_file``              additional bbcontacts file. Requires normal contactfile
``-contact_file``                 residue contact file in `CASP RR`_ format
``-disulfide_constraints_file``   disulfide residue constraints for ab initio modelling
``-distance_to_neighbour``        minimum distance between residue pairs for contact (default: 5)
``-energy_function``              Rosetta energy function for contact restraint conversion (default: *FADE*)
``-native_cutoff``                distance cutoff for reference contacts in native structure (default: 8A)
``-restraints_factor``            factor (* Sequence length) determining number of contact restraints to use (default: 1.0)
``-restraints_file``              residue restraints for ab initio modelling
``-restraints_weight``            additional energy weighting of restraints in Rosetta
===============================   ===========================

.. _mrbump_options:

MRBUMP/Molecular Replacement Options
------------------------------------
=============================   ===========================
Command-line flag               Description
=============================   ===========================
``-arpwarp_cycles``             the number of ArpWarp cycles to run
``-buccaneer_cycles``           the number of Bucanner rebuilding cycles to run
``-do_mr``                      run or skip the Molecular Replacement step
``-domain_all_chains_pdb``      fixed input to mr bump
``-domain_termini_distance``    distance between termini for insert domains
``-molrep_only``                only use Molrep for Molecular Replacement step in MRBUMP
``-mrbump_dir``                 path to a directory of MRBUMP jobs (see restart_pkl)
``-mr_keys``                    additional keywords for MRBUMP - are passed through without editing
``-mr_sg_all``                  try all possible space groups in PHASER Molecular Replacement step in MRBUMP
``-nmasu``                      manually specify the number of molecules in the asymmetric unit - sets the ``NMASu`` MRBUMP flag
``-phaser_kill``                time in minutes after which phaser will be killed (0 to leave running)
``-phaser_only``                only use Phaser for Molecular Replacement step in MRBUMP
``-phaser_rms``                 rms value for phaser (default-0.1)
``-shelx_cycles``               number of shelx cycles to run when rebuilding.
``-shelxe_exe``                 path to the shelxe executable
``-shelxe_rebuild``             rebuild shelxe traced pdb with buccaneer and arpwarp
``-shelxe_rebuild_arpwarp``     rebuild shelxe traced pdb with arpwarp
``-shelxe_rebuild_buccaneer``   rebuild shelxe traced pdb with buccaneer
``-use_arpwarp``                ``True`` to use arpwarp to rebuild.
``-use_buccaneer``              ``True`` to use Buccaneer
``-use_scwrl``                  remodel sidechains of the decoy models using Scwrl4
``-use_shelxe``                 ``True`` to use shelxe
=============================   ===========================

.. _rosetta_options:

ROSETTA Modelling Options
-------------------------
==============================   ===========================
Command-line flag                Description
==============================   ===========================
``-all_atom``                    all-atom Rosetta modelling (adds ``-return_full_atom true`` to arguments)
``-frags_3mers``                 path to file with pre-existing Rosetta 3mer fragments
``-frags_9mers``                 path to file with pre-existing Rosetta 3mer fragments
``-make_frags``                  set ``True`` to generate Rosetta 3mers and 9mers locally, False to import fragments
``-rg_reweight``                 set the Rosetta ``-rg_reweight`` flag to specify the radius of gyration reweight.
``-rosetta_AbinitioRelax``       path to Rosetta AbinitioRelax executable
``-ROSETTA_cluster``             location of rosetta cluster
``-rosetta_db``                  path to the Rosetta database directory
``-rosetta_dir``                 Rosetta install directory
``-rosetta_fragments_exe``       location of the Rosetta make_fragments.pl script
``-rosetta_version``             the version number of Rosetta
``-transmembrane``               Rosetta modelling for transmembrane proteins
``-transmembrane2``              Rosetta modelling for transmembrane proteins (NEW PROTOCOL)
``-transmembrane_octopusfile``   Octopus transmembrane topology predicition file
``-transmembrane_spanfile``      span file for modelling transmembrane proteins
``-transmembrane_lipofile``      Lips4 file for modelling transmembrane proteins
``-use_homs``                    select ROSETTA fragments from homologous models
==============================   ===========================


.. _CASPRR: http://predictioncenter.org/casproll/index.cgi?page=format
