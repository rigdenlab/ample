"""Wrapper module for the ConKit package"""

from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "18 Mar 2017"
__version__ = "2.1"

from distutils.version import StrictVersion

import inspect
import logging
import numpy
import os
import sys
import tempfile

from ample.modelling import energy_functions
from ample.util import ample_util

import conkit
import conkit.io
import conkit.plot

logger = logging.getLogger(__name__)


class SubselectionAlgorithm(object):
    """A class to collect all subselection algorithms"""
    @staticmethod
    def _numpify(data):
        """Convert a Python array to a Numpy array"""
        if type(data).__module__ == numpy.__name__:
            return data
        else:
            return numpy.asarray(data)

    @staticmethod
    def cutoff(data, cutoff=0.287):
        """A cutoff-defined subselection algorithm

        Description
        -----------
        This algorithm removes a decoy, if its score is l
        ess than the cutoff.

        Parameters
        ----------
        data : list, tuple
           A 1D array of scores
        cutoff : float, optional
           The cutoff of keeping decoys

        Returns
        -------
        list
           The decoy indices to keep
        list
           The decoy indices to throw

        """
        data = SubselectionAlgorithm._numpify(data)
        keep = numpy.where(data >= cutoff)[0]
        throw = numpy.where(data < cutoff)[0]
        return keep.tolist(), throw.tolist()

    @staticmethod
    def linear(data, cutoff=0.5):
        """A linearly-defined subselection algorithm

        Description
        -----------
        This algorithm removes the worst 500 decoys.

        Parameters
        ----------
        data : list, tuple
           A 1D array of scores
        cutoff : float, optional
           The porportion of the total number of decoys to keep

        Returns
        -------
        list
           The decoy indices to keep
        list
           The decoy indices to throw

        """
        sorted_indices = SubselectionAlgorithm._numpify(data).argsort()[::-1]
        point = numpy.ceil(sorted_indices.shape[0] * cutoff).astype(numpy.int)
        keep = sorted_indices[:point]
        throw = sorted_indices[point:]
        return keep.tolist(), throw.tolist()

    @staticmethod
    def scaled(data, cutoff=0.5):
        """A scaling-defined subselection algorithm

        Description
        -----------
        This algorithm removes a decoy, if its scaled score
        is less than 0.5. The scaled score is calculated by
        dividing the satisfaction score by the average of the
        set.

        Parameters
        ----------
        data : list, tuple
           A 1D array of scores
        cutoff : float, optional
           The cutoff of keeping decoys

        Returns
        -------
        list
           The decoy indices to keep
        list
           The decoy indices to throw

        """
        data = SubselectionAlgorithm._numpify(data)
        data_scaled = data / numpy.mean(data)
        keep = numpy.where(data_scaled >= cutoff)[0]
        throw = numpy.where(data_scaled < cutoff)[0]
        return keep.tolist(), throw.tolist()

    @staticmethod
    def ignore(data):
        """"A subselection algorithm to keep all

        Description
        -----------
        This algorithm doesn't do anything except mimic others.

        It will not discard any decoys and keep all!!

        Parameters
        ----------
        data : list, tuple
           A 1D array of scores

        Returns
        -------
        list
           The decoy indices to keep
        list
           The decoy indices to throw

        """
        data = SubselectionAlgorithm._numpify(data)
        keep = numpy.where(data >= 0)[0]
        throw = numpy.where(data < 0)[0]
        return keep.tolist(), throw.tolist()


# Populate the available subselection modes into a list
SUBSELECTION_MODES = [func_name for func_name, _ in inspect.getmembers(SubselectionAlgorithm)
                      if not func_name.startswith('_')]


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

    def __init__(self, sequence_file, sequence_format, contact_file=None, contact_format=None, bbcontacts_file=None, bbcontacts_format="bbcontacts", cutoff_factor=1.0, distance_to_neighbor=5):
        """Initialise a new :obj:`ContactUtil` instance

        Parameters
        ----------
        sequence_file : str
           The path to the sequence file
        sequence_format : str
           The format of the ``sequence_file``
        contact_file : str, optional
           The path to the contact file
        contact_format : str, optional
           The format of ``contact_file``
        bbcontacts_file : str, optional
           The path to the bbcontacts contact file
        bbcontacts_format : str
           The format of the ``bbcontacts_file``
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

        self.sequence_file = sequence_file
        self.sequence_format = sequence_format
        self.contact_file = contact_file
        self.contact_format = contact_format
        self.bbcontacts_file = bbcontacts_file
        self.bbcontacts_format = bbcontacts_format

        self.cutoff_factor = cutoff_factor
        self.distance_to_neighbor = distance_to_neighbor

    def predict_contacts_from_sequence(self, wdir=".", min_neff=200):
        from conkit.applications import CCMpredCommandline, HHblitsCommandline
        tag = os.path.basename(self.sequence_file).rsplit(".", 1)[0]
        a3m = os.path.join(wdir, tag + ".a3m")
        jon = os.path.join(wdir, tag + ".jon")
        hhblitsdb = os.environ["HHBLITSDB"]
        conkit.applications.HHblitsCommandline(
            cmd="hhblits", database=hhblitsdb, input=self.sequence_file,
            output=os.devnull, oa3m=a3m, niterations=3, id=99,
            show_all=True, cov=60, diff='inf', maxfilt=500000
        )()
        if conkit.io.read(a3m, "a3m").neff >= min_neff:
            logger.info("More than 200 effective sequences in alignment, "
                        + "predicting contacts ...")
            conkit.io.convert(a3m, "a3m", jon, "jones")
            mat = os.path.join(wdir, tag + ".mat")
            rr = os.path.join(wdir, tag + ".rr")
            conkit.applications.CCMpredCommandline(
                cmd="ccmpred", alnfile=jon, matfile=mat, renormalize=True
            )()
            conkit.io.convert(mat, "ccmpred", rr, "casprr")
            self.contact_file = rr
            self.contact_format = "casprr"
        else:
            logger.info("Less than 200 effective sequences in alignment, "
                        + "not predicting contacts ...")

    def subselect_decoys(self, decoys, decoy_format, mode='linear', subdistance_to_neighbor=24, **kwargs):
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
            * cutoff: keep all decoys with satisfaction scores of >= 0.287
            * ignore: keep all decoys
        subdistance_to_neighbor : int, optional
           The minimum distance between neighboring residues in the subselection [default: 24]
        **kwargs
           Job submission related keyword arguments

        Returns
        -------
        list
           A 2-D list of paths and scores of all sub-selected decoys

        """
        from ample.util import ample_util
        from ample.util import workers_util

        # Compute the long range contact satisfaction on a per-decoy basis
        logger.info(
            'Long-range contacts are defined with sequence separation of 24+')

        # Hack a custom copy of the contact map together that we can use with the script
        # All decoys should be sequence identical and thus we can just match it to the top
        contact_map = self.contact_map

        contact_map.match(
            conkit.io.read(
                decoys[0], decoy_format
            ).top_map, inplace=True
        )
        tmp_contact_file = tempfile.NamedTemporaryFile(delete=False)
        conkit.io.write(tmp_contact_file.name, 'casprr', contact_map)

        executable = 'conkit-precision.bat' \
            if sys.platform.startswith('win') \
            else 'conkit-precision'

        job_scripts, log_files = [], []
        for decoy in decoys:
            decoy_name = os.path.splitext(os.path.basename(decoy))[0]
            contact_name = os.path.splitext(
                os.path.basename(self.contact_file)
            )[0]

            # TODO: Get the log file business working properly
            cmd = [executable, '-d', subdistance_to_neighbor]
            if StrictVersion(conkit.__version__) <= StrictVersion('0.6.3'):
                cmd += [decoy]
            else:
                cmd += [decoy, decoy_format]
            cmd += [self.sequence_file, self.sequence_format]
            cmd += [tmp_contact_file.name, 'casprr']

            prefix = '{0}_{1}_'.format(contact_name, decoy_name)
            script = tempfile.NamedTemporaryFile(prefix=prefix, suffix=ample_util.SCRIPT_EXT,
                                                 delete=False)
            script.write(
                ample_util.SCRIPT_HEADER + os.linesep
                + " ".join(map(str, cmd)) + os.linesep
            )
            script.close()

            os.chmod(script.name, 0o777)
            job_scripts.append(script.name)
            log_files.append(os.path.splitext(script.name)[0] + ".log")

        success = workers_util.run_scripts(
            job_scripts=job_scripts,
            monitor=None,
            check_success=None,
            early_terminate=None,
            nproc=kwargs['nproc'] if 'nproc' in kwargs else 1,
            job_time=7200,          # Might be too long/short, taken from Rosetta modelling
            job_name='subselect',
            submit_cluster=kwargs['submit_cluster'] if 'submit_cluster' in kwargs else False,
            submit_qtype=kwargs['submit_qtype'] if 'submit_qtype' in kwargs else None,
            submit_queue=kwargs['submit_queue'] if 'submit_queue' in kwargs else False,
            submit_array=kwargs['submit_array'] if 'submit_array' in kwargs else None,
            submit_max_array=kwargs['submit_max_array'] if 'submit_max_array' in kwargs else None,
        )

        if not success:
            msg = "Error running decoy subselection"
            raise RuntimeError(msg)

        scores = numpy.zeros(len(decoys))
        data = zip(decoys, log_files, job_scripts)
        for i, (decoy, log, script) in enumerate(data):
            for line in open(log, 'r'):
                if line.startswith('Precision score'):
                    scores[i] = float(line.strip().split()[-1])
            map(os.remove, [script, log])

        logger.info('Model selection mode: %s', mode)
        if mode == 'scaled':
            keep, throw = SubselectionAlgorithm.scaled(scores)
        elif mode == 'linear':
            keep, throw = SubselectionAlgorithm.linear(scores)
        elif mode == 'cutoff':
            keep, throw = SubselectionAlgorithm.cutoff(scores)
        elif mode == 'ignore':
            keep, throw = SubselectionAlgorithm.cutoff(scores)
        else:
            msg = "Unknown sub-selection mode: {0}".format(mode)
            logger.critical(msg)
            raise ValueError(msg)

        if len(keep) < 1:
            msg = "Number of decoys to keep is 0 - defaulting to keeping all"
            logger.warning(msg)
            keep, throw = range(len(decoys)), []

        logger.info('Excluding %d decoy(s) from ensembling', len(throw))

        return [(decoys[i], scores[i]) for i in keep]

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
        contact_map = self.contact_map

        logger.debug(structure_file)
        logger.debug(structure_format)
        if structure_file and not structure_format:
            msg = "A structure file also needs a structure format"
            logger.critical(msg)
            raise ValueError(msg)
        elif structure_file and structure_format and structure_format not in conkit.io.CONTACT_FILE_PARSERS.keys():
            msg = "Unknown structure format"
            logger.critical(msg)
            raise ValueError(msg)
        elif structure_file and structure_format:
            logger.info(
                'Provided structure file and format are: {0} - {1}'.format(
                    structure_file, structure_format)
            )
            structure_map = conkit.io.read(
                structure_file, structure_format).top_map
            contact_map.match(structure_map, inplace=True)

            # Calculate the precision score
            precision = contact_map.precision
        else:
            structure_map = None
            precision = 0.0

        conkit.plot.ContactMapFigure(contact_map, reference=structure_map,
                                     file_name=plot_file)

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

        Returns
        -------
        str
           The file the restraints were written to

        Raises
        ------
        ValueError
           Unknown restraint format
        ValueError
           Unknown Rosetta energy function
        ValueError
           Unknown SAINT2 energy function
        """
        contact_map = self.contact_map

        if restraint_format not in ['rosetta', 'saint2']:
            msg = 'Unknown restraint format: {0}'.format(restraint_format)
            logger.critical(msg)
            raise ValueError(msg)
        elif restraint_format == 'rosetta' and not hasattr(energy_functions.RosettaFunctionConstructs, energy_function):
            msg = 'Unknown Rosetta energy function: {0} for {1}'.format(
                energy_function, restraint_format)
            logger.critical(msg)
            raise ValueError(msg)
        elif restraint_format == 'saint2' and not hasattr(energy_functions.Saint2FunctionConstructs, energy_function):
            msg = 'Unknown SAINT2 energy function: {0} for {1}'.format(
                energy_function, restraint_format)
            logger.critical(msg)
            raise ValueError(msg)

        with open(restraint_file, 'w') as f_out:

            if restraint_format == 'rosetta':
                construct = getattr(
                    energy_functions.RosettaFunctionConstructs, energy_function
                ).fget(energy_functions.RosettaFunctionConstructs)

                for contact in contact_map:
                    contact_dict = contact._to_dict()
                    contact_dict['atom1'] = 'CA' if contact.res1 == 'G' else 'CB'
                    contact_dict['atom2'] = 'CA' if contact.res2 == 'G' else 'CB'
                    contact_dict['energy_bonus'] = contact.weight * 15.00
                    contact_dict['scalar_score'] = contact.scalar_score * \
                        contact.weight
                    contact_dict['sigmoid_cutoff'] = energy_functions.DynamicDistances.cutoff(
                        contact.res1, contact.res2)
                    contact_dict['sigmoid_slope'] = 1.0 / energy_functions.DynamicDistances.percentile(
                        contact.res1, contact.res2)
                    f_out.write(construct.format(**contact_dict) + os.linesep)

            elif restraint_format == 'saint2':
                construct = getattr(
                    energy_functions.Saint2FunctionConstructs, 'DEFAULT'
                ).fget(energy_functions.Saint2FunctionConstructs)

                for contact in contact_map:
                    contact_dict = contact._to_dict()
                    f_out.write(construct.format(**contact_dict) + os.linesep)

        return restraint_file

    @staticmethod
    def check_options(optd):

        if not optd['contact_file'] and optd['bbcontacts_file']:
            msg = "You must provide -contact_file when using -bbcontacts_file or use as -contact_file instead"
            logger.critical(msg)
            raise ValueError(msg)

        if optd['contact_file'] and not os.path.isfile(optd['contact_file']):
            msg = "Cannot find contact file:\n{0}".format(optd['contact_file'])
            logger.critical(msg)
            raise ValueError(msg)

        if optd['bbcontacts_file'] and not os.path.isfile(optd['bbcontacts_file']):
            msg = "Cannot find contact file:\n{0}".format(optd['contact_file'])
            logger.critical(msg)
            raise ValueError(msg)

        if optd['contact_file'] and not optd['contact_format']:
            msg = "You must define the contact file format via -contact_format"
            logger.critical(msg)
            raise ValueError(msg)

        if optd["contact_file"] and optd['contact_format'] not in conkit.io.CONTACT_FILE_PARSERS:
            msg = "The provided contact file format is not yet implemented"
            logger.critical(msg)
            raise ValueError(msg)

        if optd['restraints_format'] == 'rosetta' and optd['energy_function']:
            if not hasattr(energy_functions.RosettaFunctionConstructs, optd['energy_function']):
                msg = "Rosetta energy function {0} unavailable".format(
                    optd['energy_function'])
                logger.critical(msg)
                raise ValueError(msg)

        if optd['restraints_format'] == 'saint2' and optd['energy_function']:
            if not hasattr(energy_functions.Saint2FunctionConstructs, optd['energy_function']):
                msg = "SAINT2 energy function {0} unavailable".format(
                    optd['energy_function'])
                logger.critical(msg)
                raise ValueError(msg)

        if optd['subselect_mode'] and optd['subselect_mode'].lower() not in SUBSELECTION_MODES:
            msg = "Subselection mode not valid"
            logger.critical(msg)
            raise ValueError(msg)

    @property
    def bbcontacts_file(self):
        return self._bbcontacts_file

    @bbcontacts_file.setter
    def bbcontacts_file(self, fname):
        self._bbcontacts_file = fname

    @property
    def bbcontacts_format(self):
        return self._bbcontacts_format

    @bbcontacts_format.setter
    def bbcontacts_format(self, value):
        if value and value not in conkit.io.CONTACT_FILE_PARSERS.keys():
            raise ValueError('Unknown contact file format: {0}'.format(value))
        self._bbcontacts_format = value

    @property
    def contacts_file(self):
        return self._contact_file

    @contacts_file.setter
    def contacts_file(self, fname):
        self._contact_file = fname

    @property
    def contact_format(self):
        return self._contact_format

    @contact_format.setter
    def contact_format(self, value):
        if value and value not in conkit.io.CONTACT_FILE_PARSERS.keys():
            raise ValueError('Unknown contact file format: {0}'.format(value))
        self._contact_format = value

    @property
    def cutoff_factor(self):
        return self._cutoff_factor

    @cutoff_factor.setter
    def cutoff_factor(self, value):
        if value < 0.0:
            msg = "cutoff factor needs to be positive: {0}".format(value)
            raise ValueError(msg)
        self._cutoff_factor = float(value)

    @property
    def distance_to_neighbor(self):
        return self._distance_to_neighbor

    @distance_to_neighbor.setter
    def distance_to_neighbor(self, value):
        if value < 0:
            msg = "cutoff factor needs to be positive: {0}".format(value)
            raise ValueError(msg)
        self._distance_to_neighbor = int(value)

    @property
    def sequence_file(self):
        return self._sequence_file

    @sequence_file.setter
    def sequence_file(self, fname):
        self._sequence_file = fname

    @property
    def sequence_format(self):
        return self._sequence_format

    @sequence_format.setter
    def sequence_format(self, value):
        if value not in conkit.io.SEQUENCE_FILE_PARSERS.keys():
            raise ValueError('Unknown sequence file format: {0}'.format(value))
        self._sequence_format = value

    @property
    def contact_map(self):
        if self.contact_file is None:
            raise ValueError("No contact file provded!")

        logger.info('Provided contact file and format are: %s - %s',
                    self.contact_file, self.contact_format)
        contact_map = conkit.io.read(
            self.contact_file, self.contact_format).top_map

        logger.info('Provided sequence file and format are: %s - %s',
                    self.sequence_file, self.sequence_format)
        sequence = conkit.io.read(
            self.sequence_file, self.sequence_format).top_sequence
        contact_map.sequence = sequence
        contact_map.assign_sequence_register()

        logger.info('Calculating the scalar score')
        contact_map.calculate_scalar_score()

        dtn = self.distance_to_neighbor
        logger.info(
            'Removing neighboring residues to distance of %d residues', dtn)
        contact_map.remove_neighbors(min_distance=dtn, inplace=True)

        sort_key = 'raw_score'
        logger.info('Sorting the contact map based on %s', sort_key)
        contact_map.sort(sort_key, reverse=True, inplace=True)

        ncontacts = int(contact_map.sequence.seq_len * self.cutoff_factor)
        logger.info(
            'Slicing contact map to contain top %d contacts only', ncontacts)
        contact_map = contact_map[:ncontacts]

        if self.bbcontacts_file:
            logger.info('Provided contact file and format are: %s - %s',
                        self.bbcontacts_file, self.bbcontacts_format)
            bbcontact_map = conkit.io.read(
                self.bbcontacts_file, self.bbcontacts_format).top_map
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

    @property
    def found_ccmpred_contact_prediction_deps(self):
        from ample.util.ample_util import FileNotFoundError
        found_ccmpred_exe = False
        found_hhblits_exe = False
        try:
            ample_util.find_exe("ccmpred")
            found_ccmpred_exe = True
            ample_util.find_exe("hhblits")
            found_hhblits_exe = True
        except FileNotFoundError:
            pass
        found_hhblitsdb = "HHBLITSDB" in os.environ
        return found_ccmpred_exe and found_hhblits_exe and found_hhblitsdb

    @property
    def require_contact_prediction(self):
        return self.contact_file is None

    @property
    def do_contact_analysis(self):
        return not self.require_contact_prediction


def _create_parsers():
    DISABLE = ["casprr", "flib", "genericstructure"]
    parsers = conkit.io.CONTACT_FILE_PARSERS.keys()
    for d in DISABLE:
        if d in parsers:
            parsers.pop(parsers.index(d))
    return sorted(parsers)


CONTACT_FILE_PARSERS = _create_parsers()
