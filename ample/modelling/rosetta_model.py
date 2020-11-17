"""
Created on 21 Feb 2013

@author: jmht
"""

# Python modules
import copy
import glob
import logging
import os
import random
import shutil

# Our modules
from ample.modelling import energy_functions, octopus_predict, multimer_definitions
from ample.parsers import psipred_parser
from ample.util import ample_util, pdb_edit, sequence_util
from pyjob.factory import TaskFactory
from pyjob.script import ScriptCollector, Script

logger = logging.getLogger(__name__)


def align_mafft(query_seq, template_seq, logger, mafft_exe=None):
    if not mafft_exe:
        mafft_exe = os.path.join(os.environ['CCP4'], 'libexec', 'mafft')
        if not ample_util.is_exe(mafft_exe):
            raise RuntimeError("Cannot find CCP4 mafft binary: {0}".format(mafft_exe))

    logger.info(
        "Running mafft binary {0} to generate alignment between target fasta and the NMR model sequence.".format(
            mafft_exe
        )
    )

    name = "{0}__{1}".format(query_seq.name, template_seq.name)
    mafft_input = "{0}_concat.fasta".format(name)
    # Join two sequences
    mafft_seq = copy.copy(query_seq)
    mafft_seq += template_seq
    mafft_seq.write_fasta(mafft_input)
    cmd = [mafft_exe, '--maxiterate', '1000', '--localpair', '--quiet', mafft_input]
    logfile = os.path.abspath('mafft.out')

    # Due to a bug in the CCP4 mafft installation, we need to set MAFFT_BINARIES
    env = copy.copy(os.environ)
    env['MAFFT_BINARIES'] = os.path.join(env['CCP4'], 'libexec')
    ret = ample_util.run_command(cmd, logfile=logfile, env=env)
    if ret != 0:
        raise RuntimeError("Error running mafft for alignnment - check logfile: {0}".format(logfile))

    seq_align = sequence_util.Sequence()
    seq_align.from_fasta(logfile, canonicalise=False)

    logger.info("Got Alignment:\n{0}\n{1}".format(seq_align.sequences[0], seq_align.sequences[1]))
    logger.info("If you want to use a different alignment, import using -alignment_file")

    alignment_file = "{0}_align.grishin".format(name)
    with open(alignment_file, 'w') as w:
        # First line is query and template name - must match the name of the template pdb
        w.write(
            """## {0}  {1}
# hhsearch\n
scores_from_program: 0 1.00\n'+
0 {2}
0 {3}
""".format(
                query_seq.name, template_seq.name, seq_align.sequences[0], seq_align.sequences[1]
            )
        )
    return os.path.abspath(alignment_file)


def align_clustalw(query_seq, template_seq, logger, clustalw_exe=None):
    if not clustalw_exe:
        mafft_exe = os.path.join(os.environ['CCP4'], 'libexec', 'clustalw2')
        if not ample_util.is_exe(mafft_exe):
            raise RuntimeError("Cannot find CCP4 clustalw2 binary: {}".format(mafft_exe))

    name = "{0}__{1}".format(query_seq.name, template_seq.name)
    clustalw_input = "{0}_concat.fasta".format(name)
    query_seq.concat(template_seq, clustalw_input)
    align_out = "{0}_align.fasta".format(name)
    logfile = os.path.abspath("clustalw2.log")
    cmd = [
        clustalw_exe,
        '-align',
        '-outorder=input',
        '-output=fasta',
        '-infile={0}'.format(clustalw_input),
        '-outfile={0}'.format(align_out),
    ]

    ret = ample_util.run_command(cmd, logfile=logfile)
    if ret != 0:
        raise RuntimeError("Error running clustalw2 for alignnment - check logfile: {0}".format(logfile))
    seq_align = sequence_util.Sequence()
    seq_align.from_fasta(align_out, canonicalise=False)

    logger.info("Got Alignment:\n{0}\n{1}".format(seq_align.sequences[0], seq_align.sequences[1]))
    logger.info("If you want to use a different alignment, import using -alignment_file")

    alignment_file = "{0}_align.grishin".format(name)
    with open(alignment_file, 'w') as w:
        # First line is query and template name - must match the name of the template pdb
        w.write(
            """## {0}  {1}
# hhsearch\n
scores_from_program: 0 1.00\n'+
0 {2}
0 {3}
""".format(
                query_seq.name, template_seq.name, seq_align.sequences[0], seq_align.sequences[1]
            )
        )
    return os.path.abspath(alignment_file)


class RosettaModel(object):
    """
    Class to run Rosetta modelling
    """

    def __init__(self, optd=None, rosetta_dir=None):

        self.debug = None
        self.nproc = None
        self.nmodels = None
        self.work_dir = None  # Where the modelling happens - can be deleted on exit
        self.ample_dir = None
        self.models_dir = None
        self.rosetta_dir = None
        self.rosetta_bin = None
        self.rosetta_relax_exe = None
        self.rosetta_minirosetta_exe = None
        self.rosetta_mr_protocols = None
        self.rosetta_idealize_jd2 = None
        self.rosetta_db = None
        self.rosetta_version = None

        self.fasta = None
        self.sequence_length = None
        self.all_atom = None

        # Fragment variables
        self.name = None
        self.frags_3mers = None
        self.frags_9mers = None
        self.use_homs = None
        self.fragments_directory = None
        self.fragments_exe = None

        # Transmembrane variables
        self.transmembrane = None
        self.transmembrane_old = None
        self.tm_patch_file = None
        self.octopus2span = None
        self.run_lips = None
        self.align_blast = None
        self.nr = None
        self.blastpgp = None
        self.octopusTopology = None
        self.spanfile = None
        self.lipofile = None

        # Extra options
        self.psipred_ss2 = None
        self.domain_termini_distance = None
        self.rad_gyr_reweight = None
        self.improve_template = None
        self.multimer_modelling = None
        self.num_chains = None
        self.nativePdbStd = None
        self.rosetta_flagsfile = None
        self.restraints_file = None
        self.restraints_weight = None
        self.disulfide_constraints_file = None

        self.nmr_remodel = None
        self.nmr_process_ntimes = None
        self.nmr_alignment_file = None
        self.mnr_remodel_fasta = None

        if optd:
            self.set_paths(optd=optd, rosetta_dir=rosetta_dir)
            self.set_from_dict(optd)
        if self.work_dir and not os.path.isdir(self.work_dir):
            os.mkdir(self.work_dir)
        return

    def ab_initio_cmd(self):
        """Return the command to run rosetta as a list suitable for subprocess"""
        cmd = [
            '-database',
            self.rosetta_db,
            '-in::file::fasta',
            self.fasta,
            '-in:file:frag3',
            self.frags_3mers,
            '-in:file:frag9',
            self.frags_9mers,
            '-out:pdb',
            '-out:file:silent',
            'silent.out',
            '-run:constant_seed',
            '-abinitio:relax',
            '-relax::fast',
        ]

        if self.rosetta_version >= 3.4:
            # Recommended default paramenters - see also Radius of gyration reweight
            # we omit the -use_filters true option as we want to continue with refinement
            # even for 'failed' models as the failed models will be used anyway
            cmd += [
                "-abinitio::increase_cycles",
                "10",
                "-abinitio::rsd_wt_helix",
                "0.5",
                "-abinitio::rsd_wt_loop",
                "0.5",
            ]

            if self.psipred_ss2:  # not sure if this works < 3.4
                cmd += ["-psipred_ss2", self.psipred_ss2]

        if self.all_atom:
            cmd += ['-return_full_atom', 'true']
        else:
            cmd += ['-return_full_atom', 'false']

        if self.transmembrane_old:
            cmd += [
                '-in:file:spanfile',
                self.spanfile,
                '-in:file:lipofile',
                self.lipofile,
                '-abinitio:membrane',
                '-membrane:no_interpolate_Mpair',
                '-membrane:Menv_penalties',
                '-score:find_neighbors_3dgrid',
                '-membrane:normal_cycles',
                '40',
                '-membrane:normal_mag',
                '15',
                '-membrane:center_mag',
                '2',
                '-mute core.io.database',
                '-mute core.scoring.MembranePotential',
            ]
        elif self.transmembrane:
            cmd += ['-score:patch', self.tm_patch_file]
            if self.restraints_file and os.path.isfile(self.restraints_file):
                self.restraints_weight = 3

        # Radius of gyration reweight
        if self.rad_gyr_reweight is not None:
            cmd += ['-rg_reweight', str(self.rad_gyr_reweight)]
        else:
            cmd += ['-rg_reweight', "0.5"]

        # Add any restraints
        cmd = self.cmd_add_restraints(cmd)

        # Improve Template
        if self.improve_template:
            cmd += [
                '-in:file:native',
                self.improve_template,
                '-abinitio:steal_3mers',
                'true',
                '-abinitio:steal9mers',
                'true',
                '-abinitio:start_native',
                'true',
                '-templates:force_native_topology',
                'true',
            ]
        return cmd

    def ab_initio_model(self, processed_models=None):
        """Run the ab initio modelling and return a list of models.

        Parameters
        ----------
        processed_models : list, optional
           A list of pdb models that should have been checked for suitability.
           NB - currently only required by NMR remodelling
 
        Returns
        -------
        list
           The list of created ab initio pdb models
        
        """

        if self.nmr_remodel:
            return self.do_nmr_remodel(processed_models)
        elif self.multimer_modelling:
            return self.do_multimer_modelling()

        rosetta_executable = self.rosetta_relax_exe
        flagsfile = self.rosetta_flagsfile
        if not flagsfile:
            # Run any setup requried
            if self.transmembrane_old:
                rosetta_executable = self.transmembrane_exe
                self.tm_make_files()
            elif self.transmembrane:
                self.tm2_make_patch(self.work_dir)

            if self.domain_termini_distance > 0 and not self.multimer_modelling:
                if self.restraints_file:
                    raise RuntimeError(
                        "Cannot set up domain restraints with existing restraints file: {}!".format(
                            self.restraints_file
                        )
                    )
                self.restraints_file = self.setup_domain_restraints()

            # create the flags file with the rosetta directives
            flagsfile = os.path.join(self.work_dir, 'rosetta.flags')
            flags = self.process_cmd_list(self.ab_initio_cmd())
            with open(flagsfile, 'w') as w:
                w.write(flags)

        return self.model_from_flagsfile(flagsfile, rosetta_executable=rosetta_executable)

    def cmd_add_restraints(self, cmd):
        """Add any restraints and files to the ROSETTA command-line options"""
        if self.restraints_file:
            cmd += ['-constraints:cst_file', self.restraints_file, '-constraints:cst_fa_file', self.restraints_file]
            if self.restraints_weight is not None:
                cmd += [
                    '-constraints:cst_weight',
                    str(self.restraints_weight),
                    '-constraints:cst_fa_weight',
                    str(self.restraints_weight),
                ]

        # Add compatibility for extra disulfide restraints
        if self.disulfide_constraints_file and os.path.isfile(self.disulfide_constraints_file):
            cmd += ['-in::fix_disulf', str(self.disulfide_constraints_file)]
        return cmd

    def do_multimer_modelling(self):
        symmetry_file = self.create_multimer_symmetry_file()
        constraints_file = self.create_multimer_constraints_file()
        broker_file = self.create_broker_definition_file()
        flagsfile = self.create_multimer_flagsfile(
            broker_file=broker_file, symmetry_file=symmetry_file, constraints_file=constraints_file
        )
        models = self.model_from_flagsfile(
            flagsfile=flagsfile, rosetta_executable=self.rosetta_minirosetta_exe, consolidate_pdbs=False
        )
        return self.process_multimer_models(models)

    def process_multimer_models(self, modelsin):
        """Merge the multi-chain pdbs into a single chain"""
        # Work out how many chains we're using
        chaind = sequence_util.sequence_data(modelsin[0])
        chains = None
        if self.num_chains and self.num_chains < len(chaind):
            # for now we just assume we want continguous ones
            chains = list(chaind.keys())[: self.num_chains]
            logger.info("Selecting chains %s from multimer models", chains)
        models = []
        for i, pdbin in enumerate(modelsin):
            pdbout = os.path.join(self.models_dir, "multimermodel_{}.pdb".format(i))
            if chains and len(chains) == 1:
                pdb_edit.extract_chain(pdbin, pdbout, chainID=chains[0])
            else:
                pdb_edit.merge_chains(pdbin, pdbout, chains=chains)
            models.append(pdbout)
        return models

    def create_broker_definition_file(self):
        broker_file = os.path.join(self.work_dir, 'setup_init.tpb')
        with open(broker_file, 'w') as w:
            w.write(multimer_definitions.BROKER_SETUP_STR)
        return broker_file

    def create_multimer_symmetry_file(self):
        if self.multimer_modelling == multimer_definitions.DIMER:
            symfile_name = 'symdef_dimer.dat'
            symdef = multimer_definitions.SYMFILE_DIMER
        elif self.multimer_modelling == multimer_definitions.TRIMER:
            symfile_name = 'symdef_trimer.dat'
            symdef = multimer_definitions.SYMFILE_TRIMER
        elif self.multimer_modelling == multimer_definitions.TETRAMER:
            symfile_name = 'symdef_tetramer.dat'
            symdef = multimer_definitions.SYMFILE_TETRAMER
        else:
            raise RuntimeError("Unrecognised multimer_modelling mode: {}".format(self.multimer_modelling))

        symfile_path = os.path.join(self.work_dir, symfile_name)
        with open(symfile_path, 'w') as w:
            w.write(symdef)

        return symfile_path

    def create_multimer_constraints_file(self):
        restraints_file = os.path.join(self.work_dir, 'multimer_constraints.txt')
        if self.multimer_modelling == multimer_definitions.DIMER:
            nchains = 2
        elif self.multimer_modelling == multimer_definitions.TRIMER:
            nchains = 3
        elif self.multimer_modelling == multimer_definitions.TETRAMER:
            nchains = 4

        restraint_str = energy_functions.RosettaFunctionConstructs().FLAT_HARMONIC
        domain_termini_distance = self.sequence_length * 1.5
        import string
        import itertools

        chains = list(string.ascii_uppercase)[:nchains]
        fstr = ""
        for chain1, chain2 in itertools.combinations(chains, 2):
            # First add the terminal restraints
            optd = {
                'atom1': 'CA',
                'res1_seq': 1,
                'atom2': 'CA',
                'res2_seq': self.sequence_length,
                'x0': domain_termini_distance,
                'stddev': 3.0,
                'tol': 5.0,
            }
            fstr += restraint_str.format(**optd) + os.linesep
            # Now add the restraints between chains
            for resseq in range(1, self.sequence_length + 1):
                resseq = str(resseq)
                optd = {
                    'atom1': 'CA',
                    'res1_seq': resseq + chain1,
                    'atom2': 'CA',
                    'res2_seq': resseq + chain2,
                    'x0': 10,
                    'stddev': 3.0,
                    'tol': 5.0,
                }
                fstr += restraint_str.format(**optd) + os.linesep
        with open(restraints_file, 'w') as w:
            w.write(fstr)
        return restraints_file

    def create_multimer_flagsfile(self, broker_file=None, symmetry_file=None, constraints_file=None):
        flagsfile = os.path.join(self.work_dir, 'multimer_flagsfile.txt')
        optd = {
            'broker_file': broker_file,
            'symdef_file': symmetry_file,
            'constraint_file': constraints_file,
            'fasta_file': self.fasta,
            'frags3': self.frags_3mers,
            'frags9': self.frags_9mers,
        }
        with open(flagsfile, 'w') as w:
            w.write(multimer_definitions.FLAGSFILE_STR.format(**optd))
        return flagsfile

    def find_binary(self, name):
        """
        Find a rosetta binary on different platforms
        separate from object as it's currently used by the NMR stuff - which is in dire need of refactoring.

        """
        if not self.rosetta_bin and os.path.isdir(self.rosetta_bin):
            raise RuntimeError("Cannot find rosetta_bin directory: {}".format(self.rosetta_bin))
        binaries = glob.glob(os.path.join(self.rosetta_bin, "{0}.*".format(name)))
        if not len(binaries):
            return False
        # Could check for shortest - for now just return the first
        binary = os.path.abspath(binaries[0])
        if os.path.isfile(binary):
            return binary
        return False

    def generate_seeds(self, nseeds, start=1000000, end=4000000):
        """Generate a list of nseed seeds

        Parameters
        ----------
        nseeds : int
           The number of seeds required
        start : int
           Beginning of random range
        nseeds : int
           End of random range
            """
        assert 0 < nseeds < end - start, "Invalid seed count: {0}".format(nseeds)
        seeds = set()
        while len(seeds) < nseeds:
            seeds.add(random.randint(start, end))
        return list(seeds)

    def fragment_cmd(self):
        """
        Return the command to make the fragments as a list

        """
        # It seems that the script can't tolerate "-" in the directory name leading to the fasta file,
        # so we need to copy the fasta file into the fragments directory and just use the name here
        fasta = os.path.split(self.fasta)[1]

        cmd = [self.fragments_exe, '-rundir', self.fragments_directory, '-id', self.name, fasta]

        if self.transmembrane_old:
            # cmd += [ '-noporter', '-nopsipred','-sam']
            pass
        else:
            # version dependent flags
            if self.rosetta_version == 3.3:
                # jmht the last 3 don't seem to work with 3.4
                cmd += ['-noporter', '-nojufo', '-nosam', '-noprof']
            elif self.rosetta_version >= 3.4:
                cmd += ['-noporter']

        # Whether to exclude homologs
        if not self.use_homs:
            cmd.append('-nohoms')

        # Use user-supplied psipred_ss2 file (thanks for Felix for spotting this omission)
        if self.psipred_ss2:
            cmd += ['-nopsipred', '-psipredfile', self.psipred_ss2]

        # Be 'chatty'
        cmd.append('-verbose')

        return cmd

    def generate_fragments(self, amoptd):
        """
        Run the script to generate the fragments
        """
        logger.info('----- making fragments--------')
        if not os.path.exists(self.fragments_directory):
            os.mkdir(self.fragments_directory)

        # It seems that the script can't tolerate "-" in the directory name leading to the fasta file,
        # so we need to copy the fasta file into the fragments directory
        fasta = os.path.split(self.fasta)[1]
        shutil.copy2(self.fasta, os.path.join(self.fragments_directory, fasta))

        collector = ScriptCollector(None)
        script = Script(directory=self.fragments_directory, stem="mkfrags")
        script.append(" ".join(self.fragment_cmd()))
        collector.add(script)

        success = self.run_scripts(collector, run_dir=self.fragments_directory, job_time=21600, monitor=None)
        logfile = "{0}.log".format(os.path.splitext(script)[0])
        if not success:
            logfile = "{0}.log".format(os.path.abspath(os.path.splitext(os.path.basename(script))[0]))
            msg = "Error generating fragments!\nPlease check the logfile {0}".format(logfile)
            raise RuntimeError(msg)

        if self.rosetta_version >= 3.4:
            # new name format: $options{runid}.$options{n_frags}.$size" . "mers
            self.frags_3mers = os.path.join(self.fragments_directory, self.name + '.200.3mers')
            self.frags_9mers = os.path.join(self.fragments_directory, self.name + '.200.9mers')
        else:
            # old_name_format: aa$options{runid}$fragsize\_05.$options{n_frags}\_v1_3"
            self.frags_3mers = os.path.join(self.fragments_directory, 'aa' + self.name + '03_05.200_v1_3')
            self.frags_9mers = os.path.join(self.fragments_directory, 'aa' + self.name + '09_05.200_v1_3')

        if not os.path.exists(self.frags_3mers) or not os.path.exists(self.frags_9mers):
            raise RuntimeError(
                "Error making fragments - could not find fragment files:\n{}\n{}\nPlease check logfile {} for details.".format(
                    self.frags_3mers, self.frags_9mers, logfile
                )
            )

        logger.info('Fragments Done\n3mers at: ' + self.frags_3mers + '\n9mers at: ' + self.frags_9mers + '\n\n')

        psipred_ss2 = os.path.join(self.fragments_directory, self.name + '.psipred_ss2')
        if os.path.exists(psipred_ss2):
            psipred_parser.PsipredSs2Parser(psipred_ss2).check_content()
            self.psipred_ss2 = psipred_ss2

        return

    def get_version(self):
        """ Return the Rosetta version as a string"""
        # Get version
        version = None
        version_file = os.path.join(self.rosetta_dir, 'README.version')
        if os.path.exists(version_file):
            try:
                for line in open(version_file, 'r'):
                    line.strip()
                    if line.startswith('Rosetta'):
                        tversion = line.split()[1].strip()
                        # version can be 3 digits - e.g. 3.2.4 - we only care about 2
                        version = float(".".join(tversion.split(".")[0:2]))
                # logger.info( 'Your Rosetta version is: {0}'.format( version ) )
            except Exception as e:
                logger.critical("Error determining rosetta version from file: {0}\n{1}".format(version_file, e))
                return False
        else:
            # Version file is absent in 3.5, so we need to use the directory name
            logger.debug('Version file for Rosetta not found - checking to see if its 3.5 or 3.6')
            self.rosetta_dir = self.rosetta_dir.rstrip(os.sep)
            if self.rosetta_dir.endswith("3.5"):
                version = 3.5
            # 3.6 bundles seem to look like: rosetta_2014.30.57114_bundle
            # elif re.search("rosetta_\d{4}\.\d{2}\.\d{5}_bundle",dirname):
            # Ignore as people change the directory names - just check for the presence of the folders:
            elif self._chk36(self.rosetta_dir):
                version = 3.6
            else:
                logger.debug("Cannot determine rosetta version in directory: %s", self.rosetta_dir)
                return False
        logger.info('Rosetta version is: {0}'.format(version))
        return version

    def _chk36(self, rosetta_dir):
        # Make sure all the expected directories are found
        # expected=frozenset(['demos','documentation','main','tools'])
        expected = frozenset(['demos', 'main', 'tools'])
        found = frozenset(
            [os.path.basename(d) for d in os.listdir(rosetta_dir) if os.path.isdir(os.path.join(rosetta_dir, d))]
        )
        if len(expected.intersection(found)) == len(expected):
            return True
        else:
            return False

    def get_bin_dir(self):
        """Determine the binary directory for the version"""
        assert (
            self.rosetta_version and type(self.rosetta_version) is float
        ), "self.rosetta_version needs to be set before calling get_bin_dir"
        assert os.path.isdir(self.rosetta_dir), "self.rosetta_dir needs to have been set before calling get_bin_dir"
        bin_dir = os.path.join(self.rosetta_dir, 'rosetta_source', 'bin')
        if self.rosetta_version >= 3.6:
            bin_dir = os.path.join(self.rosetta_dir, 'main', 'source', 'bin')
        return bin_dir

    def idealize_cmd(self, pdbin):
        """Return command to idealize pdbin"""
        return [self.rosetta_idealize_jd2, "-database", self.rosetta_db, "-s", pdbin]

    def idealize_models(self, models):
        # Loop through each model, idealise them and get an alignment
        owd = os.getcwd()
        idealise_dir = os.path.join(self.work_dir, 'idealised_models')
        os.mkdir(idealise_dir)
        os.chdir(idealise_dir)
        logger.info("Idealising {0} models in directory: {1}".format(len(models), idealise_dir))
        collector = ScriptCollector(None)
        id_pdbs = []
        for model in models:
            name = os.path.splitext(os.path.basename(model))[0]
            script = Script(directory=idealise_dir, prefix=name, stem="_idealize")

            script.append(" ".join(self.idealize_cmd(pdbin=model)))
            # Get the name of the pdb that will be output
            id_pdbs.append(self.idealize_pdbout(pdbin=model, directory=idealise_dir))
            collector.add(script)

        # Run the jobs
        success = self.run_scripts(collector, run_dir=idealise_dir, job_time=7200, job_name='idealize', monitor=None)
        if not success:
            raise RuntimeError(
                "Error running ROSETTA in directory: {0}\nPlease check the log files for more information.".format(
                    idealise_dir
                )
            )
        os.chdir(owd)
        return id_pdbs

    def idealize_pdbout(self, pdbin, directory=None):
        """Return the path to the pdb generated by idealize for pdbin"""
        pdir, fname = os.path.split(pdbin)
        if not directory:
            directory = pdir
        name, ext = os.path.splitext(fname)
        return os.path.join(directory, "{0}_0001{1}".format(name, ext))

    def model_from_flagsfile(self, flagsfile, rosetta_executable=None, job_time=43200, consolidate_pdbs=True):
        """Run ROSETTA modelling from a flagsfile"""
        assert os.path.isfile(flagsfile), "Cannot find ROSETTA flagsfile: {}".format(flagsfile)
        # Remember starting directory
        owd = os.getcwd()
        os.chdir(self.work_dir)
        if not os.path.isdir(self.models_dir):
            os.mkdir(self.models_dir)

        # Split jobs onto separate processors - 1 for cluster, as many as will fit for desktop
        if self.submit_qtype != 'local':
            jobs_per_proc = [1] * self.nmodels
        else:
            jobs_per_proc = self.split_jobs(self.nmodels, self.nproc)
        if rosetta_executable is None:
            rosetta_executable = self.rosetta_relax_exe
        seeds = self.generate_seeds(len(jobs_per_proc))
        collector = ScriptCollector(None)
        dir_list = []
        for i, njobs in enumerate(jobs_per_proc):
            if njobs < 1:
                continue
            d = os.path.join(self.work_dir, "job_{0}".format(i))
            os.mkdir(d)
            dir_list.append(d)
            script = Script(directory=d, prefix="job_", stem=str(i))
            script.append("cd {}".format(d))
            script.append("{} \\".format(rosetta_executable))
            script.append("-database {} \\".format(self.rosetta_db))
            script.append("@{} \\".format(flagsfile))
            script.append("-out:nstruct {} \\".format(njobs))
            script.append("-run:constant_seed \\")
            script.append("-run:jran {}".format(seeds[i]))
            script.append("cd {}".format(self.work_dir))
            collector.add(script)

        success = self.run_scripts(collector, run_dir=self.work_dir, job_time=job_time, job_name='abinitio')
        if not success:
            raise RuntimeError(
                "Error running ROSETTA in directory: {0}\nPlease check the log files for more information.".format(
                    self.work_dir
                )
            )

        # Copy the models into the models directory - need to rename them accordingly
        _pdbs = []
        for d in dir_list:
            ps = glob.glob(os.path.join(d, "*.pdb"))
            _pdbs += ps

        if consolidate_pdbs:
            # Copy files into the models directory
            pdbs = []
            for i, pdbin in enumerate(_pdbs):
                pdbout = os.path.join(self.models_dir, "model_{0}.pdb".format(i))
                shutil.copyfile(pdbin, pdbout)
                pdbs.append(pdbout)
        else:
            pdbs = _pdbs

        if not pdbs:
            msg = (
                "No models created after modelling!"
                + os.linesep
                + "Please check the log files in the directory {0} "
                + "for more information."
            )
            raise RuntimeError(msg.format(self.work_dir))

        if len(pdbs) != self.nmodels:
            msg = (
                "Expected to create {0} models but found {1}."
                + os.linesep
                + "Please check the log files in the directory {2} "
                "for more information."
            )
            raise RuntimeError(msg.format(self.nmodels, len(pdbs), self.work_dir))

        os.chdir(owd)  # Go back to where we came from
        return pdbs

    def mr_cmd(self, template, alignment, nstruct, seed):
        cmd = [
            self.rosetta_mr_protocols,
            '-database ',
            self.rosetta_db,
            '-MR:mode',
            'cm',
            '-in:file:extended_pose',
            '1',
            '-in:file:fasta',
            self.fasta,
            '-in:file:alignment',
            alignment,
            '-in:file:template_pdb',
            template,
            '-loops:frag_sizes',
            '9 3 1',
            '-loops:frag_files',
            self.frags_9mers,
            self.frags_3mers,
            'none',
            '-loops:random_order',
            '-loops:random_grow_loops_by',
            '5',
            '-loops:extended',
            '-loops:remodel',
            'quick_ccd',
            '-loops:relax',
            'relax',
            '-relax:default_repeats',
            '4',
            '-relax:jump_move',
            'true',
            '-cm:aln_format',
            'grishin',
            '-MR:max_gaplength_to_model',
            '8',
            '-nstruct',
            str(nstruct),
            '-ignore_unrecognized_res',
            '-overwrite',
            '-run:constant_seed',
            '-run:jran',
            str(seed),
        ]
        cmd = self.cmd_add_restraints(cmd)
        return cmd

    def do_nmr_remodel(self, models):
        if self.mnr_remodel_fasta and not os.path.isfile(self.mnr_remodel_fasta):
            raise RuntimeError("Cannot find remodel_fasta: {0}".format(self.mnr_remodel_fasta))
        if self.nmr_process_ntimes and not isinstance(self.nmr_process_ntimes, int):
            raise RuntimeError("ntimes is not an int: {0}".format(self.nmr_process_ntimes))
        num_nmr_models = len(models)
        if not self.nmr_process_ntimes:
            self.nmr_process_ntimes = 1000 / num_nmr_models
        num_models = self.nmr_process_ntimes * num_nmr_models
        logger.info('Processing each model %d times. %d models will be made', self.nmr_process_ntimes, num_models)

        # Idealize all the nmr models to have standard bond lengths, angles etc
        id_pdbs = self.idealize_models(models)
        logger.info('{0} models successfully idealized'.format(len(id_pdbs)))

        owd = os.getcwd()
        remodel_dir = os.path.join(self.work_dir, 'remodelling')
        os.mkdir(remodel_dir)
        os.chdir(remodel_dir)

        # Sequence object for idealized models
        id_seq = sequence_util.Sequence(pdb=id_pdbs[0])

        # Get the alignment for the structure - assumes all models have the same sequence
        if not self.nmr_alignment_file:
            # fasta sequence of first model
            remodel_seq = sequence_util.Sequence(fasta=self.mnr_remodel_fasta)
            alignment_file = align_mafft(remodel_seq, id_seq, logger)
        else:
            alignment_file = self.nmr_alignment_file

        # Remodel each idealized model nmr_process times
        pdbs_to_return = self.remodel(id_pdbs, self.nmr_process_ntimes, alignment_file)

        os.chdir(owd)
        return pdbs_to_return

    @staticmethod
    def process_cmd_list(cmds):
        """Create a string from a list of commands"""
        cmd_str = ""
        for i, c in enumerate(cmds):
            if c.startswith('-') and i > 0:
                cmd_str += os.linesep
            cmd_str += c + " "
        return cmd_str

    def remodel(self, id_pdbs, ntimes, alignment_file):
        remodel_dir = os.getcwd()
        proc_map = self.remodel_proc_map(id_pdbs, ntimes)
        seeds = self.generate_seeds(len(proc_map))
        dir_list = []
        job_time = 7200
        collector = ScriptCollector(None)
        for i, (id_model, nstruct) in enumerate(proc_map):
            name = "job_{0}".format(i)
            d = os.path.join(remodel_dir, name)
            os.mkdir(d)
            dir_list.append(d)  # job script
            script = Script(directory=d, stem=name)
            cmd = self.mr_cmd(template=id_model, alignment=alignment_file, nstruct=nstruct, seed=seeds[i])
            script.append("cd {}".format(d))
            script.append(" \\\n".join(cmd))
            script.append("cd {}".format(remodel_dir))
            collector.add(script)

        success = self.run_scripts(collector, run_dir=remodel_dir, job_time=job_time, job_name='remodel', monitor=None)
        if not success:
            msg = (
                "Error running ROSETTA in directory: {0}."
                + os.linesep
                + "Please check the log files for more information."
            )
            raise RuntimeError(msg.format(remodel_dir))

        # Copy the models into the models directory - need to rename them accordingly
        pdbs = []
        for d in dir_list:
            ps = glob.glob(os.path.join(d, "*.pdb"))
            pdbs += ps

        if not len(pdbs):
            msg = (
                "No pdbs after remodelling in directory: {0}."
                + os.linesep
                + "Please check the log files for more information."
            )
            raise RuntimeError(msg.format(remodel_dir))

        # We also need to strip the H atoms as the remodelling occasionally fails to model an H atom,
        # which leads to pdbs with differing numbers of atoms, which causes theseus to fail. However
        # H-atoms are unlikely to be useful for our purposes.
        pdbs_moved = []
        for i, remodelled_pdb in enumerate(pdbs):
            final_pdb = os.path.join(self.models_dir, "model_{0}.pdb".format(i))
            logger.debug(
                "Stripping H atoms from remodelled pdb {0} to create final model: {1}".format(remodelled_pdb, final_pdb)
            )
            pdb_edit.strip(remodelled_pdb, final_pdb, hydrogen=True)
            pdbs_moved.append(final_pdb)
        return pdbs_moved

    def remodel_proc_map(self, id_pdbs, ntimes):
        if self.submit_qtype != 'local':
            # For clusters we saturate the queue with single model jobs (ideally in batch mode) so that cluster
            # can manage the usage for us
            proc_map = [(pdb, 1) for pdb in id_pdbs for _ in range(ntimes)]
        else:
            len_id_pdbs = len(id_pdbs)
            if self.nproc < len_id_pdbs:
                # if we have fewer processors then pdbs to remodel, each job is an idealised pdb that will be processed
                # on a processor nmr_process times
                proc_map = [(pdb, ntimes) for pdb in id_pdbs]
            else:
                proc_per_pdb = self.nproc / len(
                    id_pdbs
                )  # ignore remainder - we're dealing with a shed-load of cpus here
                if proc_per_pdb == 0:  # More cpus then jobs, so we just assign one to each as in parallel case
                    proc_map = [(pdb, 1) for pdb in id_pdbs for _ in range(ntimes)]
                else:
                    jobs_per_proc = self.split_jobs(ntimes, proc_per_pdb)
                    proc_map = []
                    for pdb in id_pdbs:
                        for count in jobs_per_proc:
                            proc_map.append(
                                (pdb, count)
                            )  # We need to split things so that each processor does a chunk of the work
        return proc_map

    def run_scripts(self, collector, run_dir=None, job_time=None, job_name=None, monitor=None):
        with TaskFactory(
                self.submit_qtype,
                collector,
                directory=run_dir,
                environment=self.submit_pe,
                run_time=job_time,
                name=job_name,
                processes=self.nprocesses,
                max_array_size=self.submit_max_array,
                queue=self.submit_queue,
                shell="/bin/bash",
        ) as task:
            task.run()
            task.wait(interval=5, monitor_f=monitor)

        return task.completed

    def setup_domain_restraints(self):
        """Create the file for restricting the domain termini and return the path to the file."""
        logger.info('restricting termini distance: {0}'.format(self.domain_termini_distance))
        restraints_file = os.path.join(self.work_dir, 'domain_constraints.txt')
        optd = {
            'atom1': 'CA',
            'res1_seq': 1,
            'atom2': 'CA',
            'res2_seq': self.sequence_length,
            'mean': self.domain_termini_distance,
            'stddev': 5.0,
        }
        restraint = energy_functions.RosettaFunctionConstructs().GAUSSIAN.format(**optd)
        with open(restraints_file, "w") as w:
            w.write(restraint + os.linesep)
        return restraints_file

    def set_from_dict(self, optd):
        """Set the values from a dictionary"""

        # Common variables
        self.fasta = optd['fasta']
        self.sequence_length = optd['fasta_length']
        self.name = optd['name']
        self.nmodels = optd['nmodels']
        self.nproc = optd['nproc']

        # Directories
        self.ample_dir = optd['work_dir']
        self.work_dir = os.path.join(self.ample_dir, 'modelling')
        if not optd['models_dir']:
            self.models_dir = os.path.join(self.ample_dir, "models")
        else:
            self.models_dir = optd['models_dir']

        # psipred secondary structure prediction
        if optd['psipred_ss2']:
            if not os.path.isfile(optd['psipred_ss2']):
                raise RuntimeError("Cannot find psipred_ss2 file: {}".format(optd['psipred_ss2']))
            self.psipred_ss2 = optd['psipred_ss2']

        if not optd['make_frags']:
            self.frags_3mers = optd['frags_3mers']
            self.frags_9mers = optd['frags_9mers']
            if not (os.path.exists(self.frags_3mers) and os.path.exists(self.frags_9mers)):
                raise RuntimeError(
                    "Cannot find both fragment files:\n{0}\n{1}\n".format(self.frags_3mers, self.frags_9mers)
                )
        else:
            # Fragment variables
            self.use_homs = optd['use_homs']
            self.fragments_directory = os.path.join(self.work_dir, 'rosetta_fragments')

        # Extra modelling options
        self.all_atom = optd['all_atom']
        self.domain_termini_distance = optd['domain_termini_distance']
        self.rad_gyr_reweight = optd['rg_reweight']

        if optd['improve_template']:
            if not os.path.exists(optd['improve_template']):
                raise RuntimeError('cant find template to improve')
            self.improve_template = optd['improve_template']
        if optd['restraints_file']:
            if not os.path.exists(optd['restraints_file']):
                raise RuntimeError("Cannot find restraints file: {0}".format(optd['restraints_file']))
            self.restraints_file = optd['restraints_file']
        self.restraints_weight = optd['restraints_weight']
        if optd['disulfide_constraints_file']:
            if not os.path.exists(optd['disulfide_constraints_file']):
                raise RuntimeError(
                    "Cannot find disulfide constraints file: {0}".format(optd['disulfide_constraints_file'])
                )
            self.disulfide_constraints_file = optd['disulfide_constraints_file']
        if optd['rosetta_flagsfile']:
            self.rosetta_flagsfile = optd['rosetta_flagsfile']

        # NMR options
        self.nmr_remodel = optd['nmr_remodel']
        self.nmr_process_ntimes = optd['nmr_process']
        self.nmr_alignment_file = optd['alignment_file']
        self.mnr_remodel_fasta = optd['nmr_remodel_fasta']

        # Multimer modelling
        self.multimer_modelling = optd['multimer_modelling']
        self.num_chains = optd['nmasu']

        # Runtime options
        self.submit_qtype = optd['submit_qtype']
        self.nprocesses = optd['nproc']
        self.submit_max_array = optd['submit_max_array']
        self.submit_queue = optd['submit_queue']
        self.submit_array = optd['submit_array']
        self.submit_pe = optd['submit_pe']

        if optd['transmembrane_old']:
            self.transmembrane_old = True
            if optd['blast_dir']:
                blastpgp = os.path.join(optd['blast_dir'], "bin/blastpgp")
                self.blastpgp = ample_util.find_exe(blastpgp)
                if self.blastpgp:
                    logger.debug("Using user-supplied blast_dir for blastpgp executable: {0}".format(self.blastpgp))

            # nr database
            if optd['nr']:
                if not os.path.exists(optd['nr'] + ".pal"):
                    raise RuntimeError(
                        "Cannot find the nr database: {0}\n"
                        "Please give the location with the nr argument to the script.".format(optd['nr'])
                    )
                else:
                    self.nr = optd['nr']
                    if self.nr:
                        logger.debug("Using user-supplied nr database: {0}".format(self.nr))

            self.spanfile = optd['transmembrane_spanfile']
            self.lipofile = optd['transmembrane_lipofile']
            self.octopusTopology = optd['transmembrane_octopusfile']

            # Check if we've been given files
            if self.octopusTopology and not (os.path.isfile(self.octopusTopology)):
                msg = "Cannot find provided transmembrane octopus topology prediction: {0}".format(self.octopusTopology)
                raise RuntimeError(msg)

            if self.spanfile and not os.path.isfile(self.spanfile):
                msg = "Cannot find provided transmembrane spanfile: {0}".format(self.spanfile)
                raise RuntimeError(msg)

            if self.lipofile and not (os.path.isfile(self.lipofile)):
                msg = "Cannot find provided transmembrane lipofile: {0}".format(self.lipofile)
                raise RuntimeError(msg)

            if (self.spanfile and not self.lipofile) or (self.lipofile and not self.spanfile):
                msg = "You need to provide both a spanfile and a lipofile"
                raise RuntimeError(msg)
        elif optd['transmembrane']:
            self.transmembrane = True
        # End transmembrane checks

        return

    def set_paths(self, optd=None, rosetta_dir=None):
        if rosetta_dir and os.path.isdir(rosetta_dir):
            self.rosetta_dir = rosetta_dir
        elif 'rosetta_dir' not in optd or not optd['rosetta_dir']:
            raise RuntimeError(
                "rosetta_dir not set - please use the -rosetta_dir flag to point at the directory where ROSETTA is installed"
            )
        elif not os.path.isdir(optd['rosetta_dir']):
            raise RuntimeError(
                "Cannot find rosetta_dir directory: {0}\nPlease set the correct rosetta_dir variable to point at the top Rosetta directory.".format(
                    optd['rosetta_dir']
                )
            )
        else:
            self.rosetta_dir = optd['rosetta_dir']

        # Determine version
        if optd and 'rosetta_version' in optd and optd['rosetta_version'] is not None:
            logger.debug('Using user-supplied Rosetta version: {0}'.format(optd['rosetta_version']))
            version = optd['rosetta_version']
        else:
            version = self.get_version()
            if not version:
                raise RuntimeError('Cannot determine Rosetta version in directory: {0}'.format(self.rosetta_dir))

        self.rosetta_version = version

        # Find the path to the binary directory
        self.rosetta_bin = self.get_bin_dir()

        # Now set all relevant paths

        # Rosetta db
        if optd and optd['rosetta_db'] and os.path.isfile(optd['rosetta_db']):
            self.rosetta_db = optd['rosetta_db']
        else:
            if self.rosetta_version < 3.6:
                self.rosetta_db = os.path.join(self.rosetta_dir, 'rosetta_database')
            else:
                self.rosetta_db = os.path.join(self.rosetta_dir, 'main', 'database')

        if not os.path.exists(self.rosetta_db):
            raise RuntimeError('cannot find Rosetta DB: {0}'.format(self.rosetta_db))

        self.rosetta_relax_exe = self.find_binary('AbinitioRelax')
        if not self.rosetta_relax_exe:
            raise RuntimeError("Cannot find ROSETTA AbinitioRelax binary in: {0}".format(self.rosetta_bin))

        # Set path to script
        if optd and optd['rosetta_fragments_exe'] and os.path.isfile(optd['rosetta_fragments_exe']):
            self.fragments_exe = optd['rosetta_fragments_exe']
        else:
            if self.rosetta_version == 3.3:
                self.fragments_exe = os.path.join(self.rosetta_dir, 'rosetta_fragments', 'make_fragments.pl')
            elif self.rosetta_version == 3.4 or self.rosetta_version == 3.5:
                self.fragments_exe = os.path.join(
                    self.rosetta_dir, 'rosetta_tools', 'fragment_tools', 'make_fragments.pl'
                )
            elif self.rosetta_version == 3.6:
                self.fragments_exe = os.path.join(self.rosetta_dir, 'tools', 'fragment_tools', 'make_fragments.pl')

        # for nmr
        self.rosetta_mr_protocols = self.find_binary('mr_protocols')
        self.rosetta_idealize_jd2 = self.find_binary('idealize_jd2')
        self.rosetta_minirosetta_exe = self.find_binary('minirosetta')

        if optd and optd['transmembrane_old']:
            self.tm_set_paths(optd)

    def split_jobs(self, njobs, nproc):
        """
        Return a list of number of jobs to run on each processor
        """
        if njobs < nproc:
            return [1] * njobs
        split_jobs = njobs / nproc  # split jobs between processors
        remainder = njobs % nproc
        jobs = []
        for _ in range(nproc):
            njobs = split_jobs
            # Separate out remainder over jobs
            if remainder > 0:
                njobs += 1
                remainder -= 1
            jobs.append(njobs)
        return jobs

    def tm_make_files(self):
        """
        Generate the various files needed for modelling transmembrane proteins

        REM the fasta as it needs to reside in this directory or the script may fail
        due to problems with parsing directory names with 'funny' characters
        """

        # Files have already been created
        if os.path.isfile(str(self.spanfile)) and os.path.isfile(str(self.lipofile)):
            logger.debug("Using given span file: {0}\n and given lipo file: {1}".format(self.spanfile, self.lipofile))
            return

        owd = os.getcwd()  # Remember where we started
        tm_dir = os.path.join(self.work_dir, 'tm_predict')
        os.mkdir(tm_dir)
        os.chdir(tm_dir)

        # It seems that the script can't tolerate "-" in the directory name leading to the fasta file,
        # so we need to copy the fasta file into the fragments directory
        fasta = os.path.split(self.fasta)[1]
        shutil.copy2(self.fasta, os.path.join(tm_dir, fasta))

        # See if we need to query the octopus server
        if os.path.isfile(str(self.octopusTopology)):
            logger.info("Using user-supplied topology prediction file: {0}".format(self.octopusTopology))
        else:
            # Query octopus server for prediction
            octo = octopus_predict.OctopusPredict()
            logger.info(
                "Generating predictions for transmembrane regions using octopus server: {0}".format(octo.octopus_url)
            )
            # fastaseq = octo.getFasta(self.fasta)
            # Problem with 3LBW prediction when remove X
            fastaseq = octo.getFasta(self.fasta)
            octo.getPredict(self.name, fastaseq, directory=tm_dir)
            self.octopusTopology = octo.topo
            logger.debug("Got topology prediction file: {0}".format(self.octopusTopology))

        # Generate span file from predict
        self.spanfile = os.path.join(tm_dir, self.name + ".span")
        logger.debug('Generating span file {0}'.format(self.spanfile))
        cmd = [self.octopus2span, self.octopusTopology]
        retcode = ample_util.run_command(cmd, logfile=self.spanfile, directory=tm_dir)
        if retcode != 0:
            raise RuntimeError("Error generating span file. Please check the log in {0}".format(self.spanfile))

        # Now generate lips file
        logger.debug('Generating lips file from span')
        collector = ScriptCollector(None)
        lips_script = Script(directory=tm_dir, stem="run_lips")
        cmd = [self.run_lips, fasta, self.spanfile, self.blastpgp, self.nr, self.align_blast]
        lips_script.append(" ".join(cmd))
        collector.add(lips_script)

        success = self.run_scripts(collector, run_dir=tm_dir, job_time=21600, monitor=None)
        # Script only uses first 4 chars to name files
        lipofile = os.path.join(tm_dir, self.name[0:4] + ".lips4")
        if not success or not os.path.exists(lipofile):
            logfile = "run_lips.log"
            raise RuntimeError("Error generating lips file {0}. Please check the logfile {1}".format(lipofile, logfile))

        self.lipofile = lipofile
        os.chdir(owd)

    def tm2_make_patch(self, work_dir):
        wts_str = """fa_atr = 0.8
fa_rep = 0.8
fa_sol = 0.0
"""
        self.tm_patch_file = os.path.abspath(os.path.join(work_dir, 'membrane.wts'))
        with open(self.tm_patch_file, 'w') as w:
            w.write(wts_str)
        return

    def tm_set_paths(self, optd):
        # Transmembrane stuff
        mem_bin = 'membrane_abinitio2'
        self.transmembrane_exe = self.find_binary(mem_bin)
        if not self.transmembrane_exe:
            raise RuntimeError("Cannot find ROSETTA {0} binary in: {1}".format(mem_bin, self.rosetta_bin))

        if self.rosetta_version < 3.6:
            tm_script_dir = os.path.join(
                self.rosetta_dir, 'rosetta_source', 'src', 'apps', 'public', 'membrane_abinitio'
            )
        else:
            tm_script_dir = os.path.join(self.rosetta_dir, 'tools', 'membrane_tools')

        self.octopus2span = os.path.join(tm_script_dir, "octopus2span.pl")
        self.run_lips = os.path.join(tm_script_dir, "run_lips.pl")
        self.align_blast = os.path.join(tm_script_dir, "alignblast.pl")

        self.octopusTopology = optd['transmembrane_octopusfile']
        self.lipofile = optd['transmembrane_lipofile']
        self.spanfile = optd['transmembrane_spanfile']

        # Check if we've been given files
        if self.octopusTopology and not (os.path.isfile(self.octopusTopology)):
            raise RuntimeError(
                "Cannot find provided transmembrane octopus topology prediction: {0}".format(self.octopusTopology)
            )

        if self.spanfile and not (os.path.isfile(self.spanfile)):
            raise RuntimeError("Cannot find provided transmembrane spanfile: {0}".format(self.spanfile))

        if self.lipofile and not (os.path.isfile(self.lipofile)):
            raise RuntimeError("Cannot find provided transmembrane lipofile: {0}".format(self.lipofile))

        if (self.spanfile and not self.lipofile) or (self.lipofile and not self.spanfile):
            raise RuntimeError("You need to provide both a spanfile and a lipofile")

        if not (self.spanfile and self.lipofile):
            # We only need to check these if we haven't been given a spanfile and lipofile

            if (
                not os.path.exists(self.octopus2span)
                or not os.path.exists(self.run_lips)
                or not os.path.exists(self.align_blast)
            ):
                raise RuntimeError(
                    "Cannot find the required executables: octopus2span.pl ,run_lips.pl and align_blast.pl in the directory\n"
                    + "{0}\nPlease check these files are in place".format(tm_script_dir)
                )

            if optd['blast_dir']:
                blastpgp = os.path.join(optd['blast_dir'], "bin/blastpgp")
                self.blastpgp = ample_util.find_exe(blastpgp)
                if self.blastpgp:
                    logger.debug("Using user-supplied blast_dir for blastpgp executable: {0}".format(self.blastpgp))

            # nr database
            if optd['nr']:
                if not os.path.exists(optd['nr'] + ".pal"):
                    raise RuntimeError(
                        "Cannot find the nr database: {0}\nPlease give the location with the nr argument to the script.".format(
                            optd['nr']
                        )
                    )
                else:
                    self.nr = optd['nr']
                    if self.nr:
                        logger.debug("Using user-supplied nr database: {0}".format(self.nr))

            if not (self.nr and self.blastpgp):
                if self.rosetta_version > 3.5:
                    blastpgp = os.path.join(self.rosetta_dir, 'tools', 'fragment_tools', 'blast', 'bin', 'blastpgp')
                    nr = os.path.join(self.rosetta_dir, 'tools', 'fragment_tools', 'databases', 'nr')
                    if not (os.path.exists(blastpgp) and os.path.exists(nr + '.pal')):
                        msg = (
                            "Cannot find blastpgp executable and nr database requried for transmembrane modelling\n"
                            + "Please try running the {0} script to download these for rosetta.".format(
                                os.path.join(self.rosetta_dir, 'tools', 'fragment_tools', 'install_dependencies.pl')
                            )
                        )
                        raise RuntimeError(msg)
                    else:
                        self.blastpgp = blastpgp
                        self.nr = nr
                        if self.blastpgp:
                            logger.debug("Using blastpgp executable: {0}".format(self.blastpgp))
                        if self.nr:
                            logger.debug("Using nr database: {0}".format(self.nr))
                else:
                    msg = (
                        "Cannot find blastpgp executable and nr database requried for transmembrane modelling\n"
                        + "Please install blast and the nr database and use the -blast_dir and -nr_dir options to ample."
                    )
                    raise RuntimeError(msg)
        return
