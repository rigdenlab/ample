"""Wrapper module for the ConKit package"""

from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "18 Mar 2017"
__version__ = "2.1"

import conkit
import logging
import numpy
import os
import sys
import tempfile

from ample.modelling import energy_functions
from ample.util import ample_util
from ample.util import workers_util

logger = logging.getLogger(__name__)


class ContactUtil(object):
    """

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
    cutoff_factor : float
       The contact list truncation factor
    distance_to_neighbor : int
       The minimum distance between contacting residues
    sequence_file : str
       The path to the sequence file
    sequence_format : str
       The format of the ``sequence_file``

    """

    def __init__(self, contact_file, contact_format, sequence_file, sequence_format, bbcontacts_file=None,
                 cutoff_factor=1.0, distance_to_neighbor=5):
        """Initialise a new :obj:`ContactUtil` instance

        Parameters
        ----------
        contact_file : str
           The path to the contact file
        contact_format : str
           The format of ``contact_file``
        sequence_file : str
           The path to the sequence file
        sequence_format : str
           The format of the ``sequence_file``
        bbcontacts_file : str, optional
           The path to the bbcontacts contact file
        cutoff_factor : float, optional
           The contact list truncation factor [default: 1.0]
        distance_to_neighbor : int, optional
           The minimum distance between contacting residues [default: 5]

        """
        self._bbcontacts_file = None
        self._bbcontacts_format = None
        self._contact_file = None
        self._contact_format = None
        self._cutoff_factor = None
        self._distance_to_neighbor = None
        self._sequence_file = None
        self._sequence_format = None

        self.bbcontacts_format = 'bbcontacts'
        self.contact_format = contact_format
        self.sequence_format = sequence_format

        self.bbcontacts_file = bbcontacts_file
        self.contact_file = contact_file
        self.sequence_file = sequence_file

        self.cutoff_factor = cutoff_factor
        self.distance_to_neighbor = distance_to_neighbor

    @property
    def bbcontacts_file(self):
        """The path to ``bbcontacts_file``"""
        return self._bbcontacts_file

    @bbcontacts_file.setter
    def bbcontacts_file(self, file):
        """Define the path to the ``bbcontacts_file``"""
        self._bbcontacts_file = file

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
    def contacts_file(self):
        """The path to ``contacts_file``"""
        return self._contact_file

    @contacts_file.setter
    def contacts_file(self, file):
        """Define the path to the ``contacts_file``"""
        if not os.path.isfile(file):
            msg = "contact file does not exist: {0}".format(file)
            raise ValueError(msg)
        self._contact_file = file

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
    def cutoff_factor(self):
        """The contact list truncation factor"""
        return self._cutoff_factor

    @cutoff_factor.setter
    def cutoff_factor(self, value):
        """Define the contact list truncation factor"""
        if value < 0.0:
            msg = "cutoff factor needs to be positive: {0}".format(value)
            raise ValueError(msg)
        self._cutoff_factor = float(value)

    @property
    def distance_to_neighbor(self):
        """The minimum distance between neighboring contacts"""
        return self._distance_to_neighbor

    @distance_to_neighbor.setter
    def distance_to_neighbor(self, value):
        """"Define the minimum distance between neighboring contacts"""
        if value < 0:
            msg = "cutoff factor needs to be positive: {0}".format(value)
            raise ValueError(msg)
        self._distance_to_neighbor = int(value)

    @property
    def sequence_file(self):
        """The path to ``sequence_file``"""
        return self._sequence_file

    @sequence_file.setter
    def sequence_file(self, file):
        """Define the path to the ``sequence_file``"""
        if not os.path.isfile(file):
            msg = "sequence file does not exist: {0}".format(file)
            raise ValueError(msg)
        self._sequence_file = file

    @property
    def sequence_format(self):
        """The format of ``sequence_file``"""
        return self._sequence_format

    @sequence_format.setter
    def sequence_format(self, value):
        """Define the format of ``sequence_format``

        Raises
        ------
        ValueError
           Unknown sequence file format

        """
        if value not in conkit.io.SEQUENCE_FILE_PARSERS.keys():
            raise ValueError('Unknown sequence file format: {0}'.format(value))
        self._sequence_format = value

    def _preprocess(self):
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
        contact_map = conkit.io.read(self.contact_file, self.contact_format).top_map

        logger.info('Provided sequence file and format are: {0} - {1}'.format(self.sequence_file, self.sequence_format))
        sequence = conkit.io.read(self.sequence_file, self.sequence_format).top_sequence
        contact_map.sequence = sequence
        contact_map.assign_sequence_register()

        logger.info('Calculating the scalar score')
        contact_map.calculate_scalar_score()

        dtn = self.distance_to_neighbor
        logger.info('Removing neighboring residues to distance of {0} residues'.format(dtn))
        contact_map.remove_neighbors(min_distance=dtn, inplace=True)

        sort_key = 'raw_score'
        logger.info('Sorting the contact map based on {0}'.format(sort_key))
        contact_map.sort(sort_key, reverse=True, inplace=True)

        ncontacts = int(contact_map.sequence.seq_len * self.cutoff_factor)
        logger.info('Slicing contact map to contain top {0} contacts only'.format(ncontacts))
        contact_map = contact_map[:ncontacts]

        if self.bbcontacts_file:
            logger.info(
                'Provided contact file and format are: {0} - {1}'.format(self.bbcontacts_file, self.bbcontacts_format)
            )
            bbcontact_map = conkit.io.read(self.bbcontacts_file, self.bbcontacts_format).top_map
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

        return contact_map

    def subselect_decoys(self, decoys, decoy_format, mode='linear', subdistance_to_neighbor=25, **kwargs):
        """Subselect decoys excluding those not satisfying long-distance restraints

        Parameters
        ----------
        decoys : list, tuple
           A list containing paths to decoy files
        decoy_format : str
           The file format of ``decoys``
        mode : str, optional
           The subselection mode to use
            * scaled: keep the decoys with scaled scores of >= 0.5
            * linear: keep the top half of decoys
        subdistance_to_neighbor : int, optional
           The minimum distance between neighboring residues in the subselection [default: 25]
        **kwargs
           Job submission related keyword arguments

        Returns
        -------
        list
           A list of paths to the sub-selected decoys

        """

        # Compute the long range contact satisfaction on a per-decoy basis
        logger.info('Long-range contacts are defined with sequence separation of 25+')

        # Hack a custom copy of the contact map together that we can use with the script
        # All decoys should be sequence identical and thus we can just match it to the top
        contact_map = self._preprocess()
        contact_map.match(conkit.io.read(decoys[0], decoy_format).top_map, inplace=True)
        tmp_contact_file = tempfile.NamedTemporaryFile(delete=False)
        conkit.io.write(tmp_contact_file.name, 'casprr', contact_map)

        # Construct the job scripts
        job_scripts = []    # Hold job scripts
        log_files = []      # Hold paths to log files
        executable = 'conkit-precision.bat' if sys.platform.startswith('win') else 'conkit-precision'
        for decoy in decoys:
            # Some file names
            decoy_name = os.path.splitext(os.path.basename(decoy))[0]
            contact_name = os.path.splitext(os.path.basename(self.contact_file))[0]
            prefix = '{0}_{1}_'.format(contact_name, decoy_name)
            # Create the run scripts
            script = tempfile.NamedTemporaryFile(prefix=prefix, suffix=ample_util.SCRIPT_EXT, delete=False)

            # Construct the command
            # TODO: Get the log file business working properly
            cmd = [executable, '-d', subdistance_to_neighbor, decoy,
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

        # Collate the scores
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
            keep, throw = ContactUtil._select_scaled(scores)
        elif mode == 'linear':
            keep, throw = ContactUtil._select_linear(scores)
        else:
            msg = "Unknown sub-selection mode: {0}".format(mode)
            logger.critical(msg)
            raise ValueError(msg)

        # Some checks
        if len(keep) < 1:
            msg = "Number of decoys to keep is 0"
            raise RuntimeError(msg)

        logger.info('Excluding {0} decoy(s) from ensembling'.format(len(throw)))

        # TODO: return the scores so we can store them in AMPLE dict
        # Return the list of decoys to keep
        return tuple([decoys[i] for i in keep])

    def summarize(self, plot_file, structure_file=None, structure_format=None, native_cutoff=8):
        """Process the contact file etc

        Parameters
        ----------
        plot_file : str
           The path to the contact map plot
        structure_file : str
           A reference structure file
        structure_format : str
           The format of ``structure_file``
        native_cutoff : int
           The distance cutoff for contact extraction from ``structure_file``

        Returns
        -------
        str
           The path to the contact map plot
        float
           The precision score, if calculated, else 0.0

        Raises
        ------
        ValueError
           A structure file also needs a structure format
        ValueError
           A structure format also needs structure file
        ValueError
           Unknown structure format

        """
        # Process the contact map according to the parameters defined here
        contact_map = self._preprocess()

        if structure_file and not structure_format:
            msg = "A structure file also needs a structure format"
            raise ValueError(msg)
        elif structure_format and not structure_file:
            msg = "A structure format also needs structure file"
            raise ValueError(msg)
        elif structure_file and structure_format and structure_format not in conkit.io.CONTACT_FILE_PARSERS.keys():
            msg = "Unknown structure format"
            raise ValueError(msg)
        elif structure_file and structure_format:
            logger.info(
                'Provided structure file and format are: {0} - {1}'.format(structure_file, structure_format)
            )
            structure_map = conkit.io.read(structure_file, structure_format).top_map
            contact_map.match(structure_map, inplace=True)

            # Calculate the precision score
            precision = contact_map.precision
        else:
            structure_map = None
            precision = 0.0

        # Draw a contact map plot
        contact_map.plot_map(reference=structure_map, file_name=plot_file)

        return plot_file, precision

    def write_restraints(self, restraint_file, restraint_format, energy_function):
        """Write a list of restraints

        Parameters
        ----------
        restraint_file : str
           The file to write the restraints to
        restraint_format : str
           The restraints format, depends primarily on the program for which the restraints will be used
        energy_function : str
           The energy function

        Raises
        ------
        ValueError
           Unknown restraint format
        ValueError
           Unknown Rosetta energy function
        ValueError
           Unknown SAINT2 energy function
        """
        # Process the contact map according to the parameters defined here
        contact_map = self._preprocess()

        if restraint_format not in ['rosetta', 'saint2']:
            msg = 'Unknown restraint format: {0}'.format(restraint_format)
            logger.critical(msg)
            raise ValueError(msg)
        elif restraint_format == 'rosetta' and not hasattr(energy_functions.RosettaFunctionConstructs, energy_function):
            msg = 'Unknown Rosetta energy function: {0} for {1}'.format(energy_function, restraint_format)
            logger.critical(msg)
            raise ValueError(msg)
        elif restraint_format == 'saint2' and not hasattr(energy_functions.Saint2FunctionConstructs, energy_function):
            msg = 'Unknown SAINT2 energy function: {0} for {1}'.format(energy_function, restraint_format)
            logger.critical(msg)
            raise ValueError(msg)

        with open(restraint_file, 'w') as f_out:

            if format == 'rosetta':
                construct = getattr(
                    energy_functions.RosettaFunctionConstructs, energy_function
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

            elif format == 'saint2':
                construct = getattr(
                    energy_functions.Saint2FunctionConstructs, 'DEFAULT'
                ).fget(energy_functions.Saint2FunctionConstructs)

                for contact in contact_map:
                    contact_dict = contact._to_dict()
                    f_out.write(construct.format(**contact_dict) + os.linesep)

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

    @staticmethod
    def _select_linear(data):
        data_sorted = sorted(list(enumerate(data)), key=lambda x: x[1], reverse=True)
        midpoint = int(round(len(data_sorted) / 2.))
        keep = zip(*data_sorted[:midpoint])[0]
        throw = zip(*data_sorted[midpoint:])[0]
        return keep, throw

    @staticmethod
    def _select_scaled(data):
        data_scaled = data / numpy.mean(data)
        keep = numpy.where(data_scaled >= 0.5)[0]
        throw = numpy.where(data_scaled < 0.5)[0]
        return tuple(keep), tuple(throw)
