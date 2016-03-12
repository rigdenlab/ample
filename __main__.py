"""main route into Ample"""

import sys

from ample import main
from ample.util import exit_util

try:
    main.Ample().run()
except Exception as e:
    msg = "Error running main AMPLE program: {0}".format(e.message)
    exit_util.exit_error(msg, sys.exc_info()[2])
