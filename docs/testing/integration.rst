.. _integration_framework:

=====================
Integration Framework
=====================

AMPLE contains it's own Integration framework to test its functionality.

If you want to test AMPLE post-installation to see whether it behaves normally, execute the following command:

**UNIX (Linux|Mac)**:

.. code-block:: bash 

   $ ccp4-python -m ample.testing integration

**Windows**:

.. code-block:: batch

   $ ccp4-python -m %CCP4%\lib\py2\ample\testing\run_tests.py integration

Regardless of the operating system, you can provide extra flags:

* ``-c`` - clean up all test files/directories
* ``-d`` - don't actually run the jobs
* ``-h`` - show this help message and exit
* ``-n`` - number of processors to run on (1 per job)
* ``-r`` - location of rosetta installation directory [*UNIX only*]
* ``-s`` - submit to a cluster queueing system
* ``-w`` - directory to run jobs in
