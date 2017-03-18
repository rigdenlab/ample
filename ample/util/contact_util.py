"""Wrapper module for the ConKit package"""

__author__ = "Felix Simkovic"
__date__ = "02 Dec 2016"
__version__ = "2.0"

import conkit
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import os

from ample.modelling import energy_functions
from ample.util import ample_util
from ample.util import workers_util

logger = logging.getLogger(__name__)


class ContactUtil(object):
    """Utility class to handle the contact information

    Attributes
    ----------
    bbcontacts_file : str
       The path to the bbcontacts contact file
    bbcontacts_format : str
       The format of ``bbcontacts_file``
    contact_file : str
       The path to the contact file
    contact_format : str
       The format of ``contact_file``
    distance_to_neighbour : int
       The minimum distance between contacting residues
    energy_function : str
       The energy function to use
    native_cutoff : int
       The distance cutoff for structure file contact extraction
    plot_file : str
       The path to the contact map plot
    plot_format : str
       The format of ``plot_file`` [default: pdf]
    precision : float
       The precision score
    restraint_factor : float
       The factor use to define the number of contacts
    restraint_file : str
       The path to the restraint file
    restraint_format : str
       The format of ``restraint_file``
    sequence_file : str
       The path to the sequence file
    sequence_format : str
       The format of ``sequence_file`` [default: fasta]
    structure_file : str
       The path to the protein structure file
    structure_format : str
       The format of ``structure_file`` [default: pdb]

    """
    def __init__(self, optd):
        """Create a new instance of the :obj:`ContactUtil`"""

        # Define user-hidden attributes
        self._bbcontacts_format = None
        self._current_contact_map = None
        self._contact_format = None
        self._energy_function = None
        self._plot_format = None
        self._restraint_format = None
        self._sequence_format = None
        self._structure_file = None
        self._structure_format = None

        # Assign the values to the user-hidden variables
        # Do not do this above - checks will not be performed otherwise
        self.bbcontacts_format = 'bbcontacts'
        self.contact_format = optd['contact_format']
        self.energy_function = optd['energy_function']
        self.plot_format = 'pdf'
        self.restraint_format = optd['restraints_format']
        self.sequence_format = 'fasta'
        self.structure_format = 'pdb'

        # Define user-exposed attributes
        self.bbcontacts_file = optd['bbcontacts_file']
        self.contact_file = optd['contact_file']
        self.plot_file = os.path.join(optd['work_dir'], optd['name'] + ".cm." + self.plot_format)
        self.restraint_file = os.path.join(optd['work_dir'], optd['name'] + ".cst")
        self.sequence_file = optd['fasta']

        self.distance_to_neighbour = optd['distance_to_neighbour']
        self.native_cutoff = 8
        self.precision = 0.0
        self.restraint_factor = optd['restraints_factor']

        ## Define the native pdb to be the right structure file.
        ## Can easily be changed via the structure_file property.
        #
        # Check for some further optional attributes
        if optd['native_pdb'] and optd['native_pdb_std']:
            self.structure_file = optd['native_pdb_std']
        elif optd['native_pdb']:
            self.structure_file = optd['native_std']
        self.structure_map = None

        # Determine the native cutoff
        if optd['native_cutoff']:
            self.native_cutoff = optd['native_cutoff']

    @property
    def bbcontacts_format(self):
        """The format of ``bbcontacts_file``"""
        return self._bbcontacts_format

    @bbcontacts_format.setter
    def bbcontacts_format(self, value):
        """Define the format of ``bbcontacts_file``

        Raises
        ------
        ValueError
           Unknown contact file format

        """
        if value != 'bbcontacts':
            raise ValueError('Unknown contact file format: {0}'.format(value))
        self._bbcontacts_format = value

    @property
    def contact_format(self):
        """The format of ``contact_file``"""
        return self._contact_format

    @contact_format.setter
    def contact_format(self, value):
        """Define the format of ``contact_file``

        Raises
        ------
        ValueError
           Unknown contact file format

        """
        if value not in conkit.io.CONTACT_FILE_PARSERS.keys():
            raise ValueError('Unknown contact file format: {0}'.format(value))
        self._contact_format = value

    @property
    def energy_function(self):
        """The energy function to use"""
        return self._energy_function

    @energy_function.setter
    def energy_function(self, value):
        """Define the value of ``energy_function``

        Raises
        ------
        ValueError
           Rosetta energy function not defined
        ValueError
           SAINT2 energy function not defined

        """
        if self.restraint_format == 'rosetta' and not hasattr(energy_functions.RosettaFunctionConstructs, value):
            raise ValueError('Rosetta energy function not defined: {0} for {1}'.format(value, self.restraint_format))
        elif self.restraint_format == 'saint2' and not hasattr(energy_functions.Saint2FunctionConstructs, value):
            raise ValueError('SAINT2 energy function not defined: {0} for {1}'.format(value, self.restraint_format))

        self._energy_function = value

    @property
    def plot_format(self):
        """The format of ``plot_file`` [default: pdf]"""
        return self._plot_format

    @plot_format.setter
    def plot_format(self, value):
        """Define the format of ``plot_file``

        Raises
        ------
        ValueError
           Unknown plot format

        """
        if value not in matplotlib.pyplot.gcf().canvas.get_supported_filetypes():
            raise ValueError('Unknown plot format: {0}'.format(value))
        self._plot_format = value

    @property
    def restraint_format(self):
        """The format of ``restraint_file``"""
        return self._restraint_format

    @restraint_format.setter
    def restraint_format(self, value):
        """Define the restraint format

        Raises
        ------
        ValueError
           Restraint format not defined

        """
        if value not in ['rosetta', 'saint2']:
            raise ValueError('Restraint format not defined: {0}'.format(value))
        self._restraint_format = value

    @property
    def sequence_format(self):
        """The format of ``sequence_file`` [default: fasta]"""
        return self._sequence_format

    @sequence_format.setter
    def sequence_format(self, value):
        """Define the format of ``sequence_file``

        Raises
        ------
        ValueError
           Unknown sequence file format

        """
        if value != 'fasta':
            raise ValueError('Unknown structure file format: {0}'.format(value))
        self._sequence_format = value

    @property
    def structure_file(self):
        """The path to the protein structure file"""
        return self._structure_file

    @structure_file.setter
    def structure_file(self, value):
        """Define the path to the structure file"""
        self._structure_file = value

    @property
    def structure_format(self):
        """The format of ``structure_file``"""
        return self._structure_format

    @structure_format.setter
    def structure_format(self, value):
        """Define the format of ``structure_file``

        Raises
        ------
        ValueError
           Unknown structure file format

        """
        if value != 'pdb':
            raise ValueError('Unknown structure file format: {0}'.format(value))
        self._structure_format = value

    @staticmethod
    def check_options(optd):
        """Function to check that all contact files are available

        Raises
        ------
        ValueError
           You must provide ``-contact_file`` when using ``-bbcontacts_file`` or use as ``-contact_file`` instead
        ValueError
           Cannot find contact file
        ValueError
           Rosetta energy function unavailable

        """
        # Make sure contact file is provided with bbcontacts_file
        if not optd['contact_file'] and optd['bbcontacts_file']:
            msg = "You must provide -contact_file when using -bbcontacts_file or use as -contact_file instead"
            logger.critical(msg)
            raise ValueError(msg)

        # Check the existence of the contact file
        if optd['contact_file'] and not os.path.isfile(optd['contact_file']):
            msg = "Cannot find contact file:\n{0}".format(optd['contact_file'])
            logger.critical(msg)
            raise ValueError(msg)

        # Check the existence of the contact file
        if optd['bbcontacts_file'] and not os.path.isfile(optd['bbcontacts_file']):
            msg = "Cannot find contact file:\n{0}".format(optd['contact_file'])
            logger.critical(msg)
            raise ValueError(msg)

        # Check that the contact file format was provided
        if optd['contact_file'] and not optd['contact_format']:
            msg = "You must define the contact file format via -contact_format"
            logger.critical(msg)
            raise ValueError(msg)

        # Check that the contact file format is defined in ConKit
        if optd['contact_format'] not in conkit.io.CONTACT_FILE_PARSERS:
            msg = "The provided contact file format is not yet implemented"
            logger.critical(msg)
            raise ValueError(msg)

        # Make sure user selected energy function is pre-defined
        if optd['restraints_format'] == 'rosetta' and optd['energy_function']:
            if not hasattr(energy_functions.RosettaFunctionConstructs, optd['energy_function']):
                msg = "Rosetta energy function {0} unavailable".format(optd['energy_function'])
                logger.critical(msg)
                raise ValueError(msg)

        if optd['restraints_format'] == 'saint2' and optd['energy_function']:
            if not hasattr(energy_functions.Saint2FunctionConstructs, optd['energy_function']):
                msg = "SAINT2 energy function {0} unavailable".format(optd['energy_function'])
                logger.critical(msg)
                raise ValueError(msg)

        if optd['subselect_mode'] and optd['subselect_mode'].lower() not in ['linear', 'scaled']:
            msg = "Subselection mode not valid"
            logger.critical(msg)
            raise ValueError(msg)

    def _preprocess(self, match=False):
        """Pre-process the data according to the data provided

        Parameters
        ----------
        match : bool
           Match the contact maps

        Returns
        -------
        :obj:`conkit.core.ContactMap`
           The modified and processed contact map

        """
        logger.info('Provided contact file and format are: {0} - {1}'.format(self.contact_file, self.contact_format))
        contact_map = conkit.io.read(self.contact_file, self.contact_format)[0]

        logger.info('Provided sequence file and format are: {0} - {1}'.format(self.sequence_file, self.sequence_format))
        sequence = conkit.io.read(self.sequence_file, self.sequence_format)[0]
        contact_map.sequence = sequence
        contact_map.assign_sequence_register()

        logger.info('Calculating the scalar score')
        contact_map.calculate_scalar_score()

        dtn = self.distance_to_neighbour
        logger.info('Removing neighboring residues to distance of {0} residues'.format(dtn))
        contact_map.remove_neighbors(min_distance=dtn, inplace=True)

        sort_key = 'raw_score'
        logger.info('Sorting the contact map based on {0}'.format(sort_key))
        contact_map.sort(sort_key, reverse=True, inplace=True)

        ncontacts = int(contact_map.sequence.seq_len * self.restraint_factor)
        logger.info('Slicing contact map to contain top {0} contacts only'.format(ncontacts))
        contact_map = contact_map[:ncontacts]

        if self.bbcontacts_file:
            logger.info(
                'Provided contact file and format are: {0} - {1}'.format(self.bbcontacts_file, self.bbcontacts_format)
            )
            bbcontact_map = conkit.io.read(self.bbcontacts_file, self.bbcontacts_format)[0]
            bbcontact_map.sequence = sequence
            bbcontact_map.assign_sequence_register()
            bbcontact_map.rescale(inplace=True)
            bbcontact_map.calculate_scalar_score()
            bbcontact_map.sort(sort_key, reverse=True, inplace=True)

            for bbcontact in bbcontact_map:
                if bbcontact.id in contact_map:
                    contact_map[bbcontact.id].weight = 2
                    for d in (1, 2):
                        alternate_positions = [
                            (bbcontact.res1_seq, bbcontact.res2_seq + d),
                            (bbcontact.res1_seq, bbcontact.res2_seq - d),
                            (bbcontact.res1_seq + d, bbcontact.res2_seq),
                            (bbcontact.res1_seq - d, bbcontact.res2_seq),
                        ]
                        for alt_pos in alternate_positions:
                            if alt_pos in contact_map:
                                contact_map[alt_pos].weight = 2
                else:
                    contact_map.add(bbcontact)

            contact_map.sort(sort_key, reverse=True, inplace=True)

        if self.structure_file and match:
            logger.info(
                'Provided structure file and format are: {0} - {1}'.format(self.structure_file, self.structure_format)
            )
            self.structure_map = conkit.io.read(self.structure_file, self.structure_format)[0]
            contact_map.match(self.structure_map, inplace=True)

        return contact_map

    def create_restraints(self):
        """Write a list of restraints"""

        # Process the contact map according to the parameters defined here
        contact_map = self._preprocess(match=False)

        with open(self.restraint_file, 'w') as f_out:

            if self.restraint_format == 'rosetta':
                construct = getattr(
                    energy_functions.RosettaFunctionConstructs, self.energy_function
                ).fget(energy_functions.RosettaFunctionConstructs)

                for contact in contact_map:
                    contact_dict = contact._to_dict()
                    contact_dict['atom1'] = 'CA' if contact.res1 == 'G' else 'CB'
                    contact_dict['atom2'] = 'CA' if contact.res2 == 'G' else 'CB'
                    contact_dict['energy_bonus'] = contact.weight * 15.00
                    contact_dict['scalar_score'] = contact.scalar_score * contact.weight
                    contact_dict['sigmoid_cutoff'] = energy_functions.DynamicDistances.cutoff(contact.res1, contact.res2)
                    contact_dict['sigmoid_slope'] = energy_functions.DynamicDistances.percentile(contact.res1, contact.res2)
                    f_out.write(construct.format(**contact_dict) + os.linesep)

            elif self.restraint_format == 'saint2':
                construct = getattr(
                    energy_functions.Saint2FunctionConstructs, 'DEFAULT'
                ).fget(energy_functions.Saint2FunctionConstructs)

                for contact in contact_map:
                    contact_dict = contact._to_dict()
                    f_out.write(construct.format(**contact_dict) + os.linesep)

            else:
                msg = 'Restraint format unknown: {0}'.format(self.restraint_format)
                logger.critical(msg)
                raise ValueError(msg)

    def subselect_decoys(self, decoys, mode='linear', **kwargs):
        """Subselect decoys excluding those not satisfying long-distance restraints

        Parameters
        ----------
        decoy_dir : list, tuple
           A list containing paths to decoy files
        mode : str, optional
           The subselection mode to use
            * scaled: keep the decoys with scaled scores of >= 0.5
            * linear: keep the top half of decoys

        Returns
        -------
        list
           A list of paths to the sub-selected decoys

        """
        import tempfile

        # Backup original data
        org_structure_file = self.structure_file
        org_structure_format = self.structure_format
        org_structure_map = self.structure_map

        # Compute the long range contact satisfaction on a per-decoy basis
        logger.info('Long-range contacts are defined with sequence separation of 25+')

        # Hack a custom copy of the contact map together that we can use with the script
        # All decoys should be sequence identical and thus we can just match it to the top
        self.structure_file = decoys[0]
        contact_map = self._preprocess(match=True)
        tmp_contact_file = tempfile.NamedTemporaryFile(delete=False)
        conkit.io.write(tmp_contact_file.name, 'casprr', contact_map)

        # Construct the job scripts
        job_scripts = []    # Hold job scripts
        log_files = []      # Hold paths to log files
        executable = 'conkit-precision.bat' if sys.platform.startswith('win') else 'conkit-precision')
        for decoy in decoys:
            # Some file names
            decoy_name = os.path.splitext(os.path.basename(decoy))[0]
            contact_name = os.path.splitext(os.path.basename(self.contact_file))[0]
            prefix = '{0}_{1}_'.format(contact_name, decoy_name)
            # Create the run scripts
            script = tempfile.NamedTemporaryFile(prefix=prefix, suffix=ample_util.SCRIPT_EXT, delete=False)

            # Construct the command
            # TODO: Get the log file business working properly
            cmd = [executable, '-d', 25, decoy,
                   self.sequence_file, self.sequence_format,
                   tmp_contact_file.name, 'casprr']

            script.write(
                ample_util.SCRIPT_HEADER + os.linesep
                + " ".join(map(str, cmd)) + os.linesep
            )
            script.close()
            os.chmod(script.name, 0o777)
            job_scripts.append(script.name)
            # Save some more information
            log_files.append(os.path.splitext(script.name)[0] + ".log")

        # Execute the scripts
        success = workers_util.run_scripts(
            job_scripts=job_scripts,
            monitor=None,
            check_success=None,
            early_terminate=None,
            nproc=kwargs['nproc'],
            job_time=7200,          # Might be too long/short, taken from Rosetta modelling
            job_name='subselect',
            submit_cluster=kwargs['submit_cluster'],
            submit_qtype=kwargs['submit_qtype'],
            submit_queue=kwargs['submit_queue'],
            submit_array=kwargs['submit_array'],
            submit_max_array=kwargs['submit_max_array'],
        )

        if not success:
            msg = "Error running decoy subselection"
            raise RuntimeError(msg)

        # Evaluate the scores
        scores = numpy.zeros(len(decoys))
        for i, (decoy, log, script) in enumerate(zip(decoys, log_files, job_scripts)):
            for line in open(log, 'r'):
                if line.startswith('Precision score'):
                    scores[i] = float(line.strip().split()[-1])
            os.unlink(log)
            os.unlink(script)

        # Subselect the decoys
        logger.info('Model selection mode: {0}'.format(mode))
        if mode == 'scaled':
            scores_scaled = scores / numpy.mean(scores)
            to_keep = numpy.where(scores_scaled >= 0.5)[0].tolist()
            to_throw = numpy.where(scores_scaled < 0.5)[0].tolist()
        elif mode == 'linear':
            scores_sorted = sorted(list(enumerate(scores)), key=lambda x: x[1], reverse=True)
            to_keep = zip(*scores_sorted[:len(scores_sorted) / 2])[0]
            to_throw = zip(*scores_sorted[len(scores_sorted) / 2:])[0]
        else:
            msg = "Unknown sub-selection mode: {0}".format(mode)
            logger.critical(msg)
            raise ValueError(msg)

        # Some checks
        if len(to_keep) < 1:
            msg = "Number of decoys to keep is 0"
            raise RuntimeError(msg)

        logger.info('Excluding {0} decoy(s) from ensembling'.format(len(to_throw)))

        # Restore object to original state
        self.structure_file = org_structure_file
        self.structure_format = org_structure_format
        self.structure_map = org_structure_map

        # TODO: return the scores so we can store them in AMPLE dict
        # Return the list of decoys to keep
        return tuple([decoys[i] for i in to_keep])

    def summarise(self):
        """Process the contact file etc"""
        # Process the contact map according to the parameters defined here
        contact_map = self._preprocess(match=True)

        # Calculate the precision score
        if self.structure_file:
            self.precision = contact_map.precision

        # Draw a contact map plot
        contact_map.plot_map(reference=self.structure_map, file_format=self.plot_format, file_name=self.plot_file)

