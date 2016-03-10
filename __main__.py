"""main route into Ample"""

import sys

from ample import main
from ample.python import ample_exit

try:
    main.Ample().run()
except Exception as e:
    msg = "Error running main AMPLE program: {0}".format(e.message)
    ample_exit.exit_error(msg, sys.exc_info()[2])
