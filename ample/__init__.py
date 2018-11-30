import os
if 'CCP4' not in os.environ:
    raise ImportError('CCP4 needs to be installed and its environment sourced!')

from ample.util import version
__version__ = version.__version__

