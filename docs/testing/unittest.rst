.. _unittest_framework:

==================
Unittest Framework
==================

AMPLE contains it's own Unittesting framework to test its functionality.

If you want to test AMPLE post-installation to see whether it behaves normally, execute the following command:

**UNIX (Linux|Mac)**:

.. code-block:: bash 

   $ ccp4-python -m ample.testing unittest

**Windows**:

.. code-block:: batch

   $ ccp4-python -m %CCP4%\lib\py2\ample\testing\run_tests.py unittest

Regardless of the operating system, you can provide extra flags:

* ``-h`` - show this help message and exit
* ``-b`` - debugging by printing print messages
* ``-v`` - level of verbosity [default: 2]


