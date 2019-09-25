"""Tests for util.config_util"""

import os
import tempfile
import unittest

from ample.util.config_util import AMPLEConfigOptions
from ample.util import argparse_util
from ample.util import options_processor
from ample.util.mrbump_util import SHELXE_MAX_PERMITTED_RESOLUTION

__author__ = "Jens Thomas"


class Test(unittest.TestCase):
    def test_process_mr_options(self):
        """Make sure we handle resolution limits around SHELXE use properly
        """
        options = AMPLEConfigOptions()
        argso = argparse_util.process_command_line(args=['-mtz', 'foo', '-fasta', 'bar'])
        options.populate(argso)

        # Make sure we turn off shelxe if resolution is outside limit
        self.assertTrue(options.d['use_shelxe'])  # default so should be true
        options.d['shelxe_rebuild'] = True
        # Set resolution below limit
        options.d['mtz_min_resolution'] = SHELXE_MAX_PERMITTED_RESOLUTION + 1
        options_processor.process_mr_options(options.d)
        self.assertFalse(options.d['use_shelxe'])
        self.assertFalse(options.d['shelxe_rebuild'])

        # Check we don't accidently turn it off if above limit
        options.d['use_shelxe'] = True
        options.d['shelxe_rebuild'] = True
        options.d['mtz_min_resolution'] = SHELXE_MAX_PERMITTED_RESOLUTION - 1
        options_processor.process_mr_options(options.d)
        self.assertTrue(options.d['use_shelxe'])
        self.assertTrue(options.d['shelxe_rebuild'])


if __name__ == "__main__":
    unittest.main()
