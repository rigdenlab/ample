#!/usr/bin/python2.6

#

import signal, time
import subprocess
import os
from signal import alarm, signal, SIGALRM, SIGKILL
from subprocess import PIPE, Popen



has_worked = False
while has_worked == False:

   p = subprocess.Popen('ls', shell = True)
   print p.pid
   time.sleep(10)
   print p.poll()
   if p.poll() is None: #still running      
      p.kill()
      print 'timed out'
   else:
     print p.communicate()
     has_worked = True
