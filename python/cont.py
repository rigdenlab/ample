#!/usr/bin/python2.6

#

import signal, time
import subprocess
import os
from signal import alarm, signal, SIGALRM, SIGKILL
from subprocess import PIPE, Popen

x = 1
while x <2:
 print  'running'
 time.sleep(2)
