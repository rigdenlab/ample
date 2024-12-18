.. _example_abinitio:

========================
Using *ab initio* models
========================

Rosetta Installation
====================
ROSETTA must be installed on the local system before AMPLE can be used to generate *ab inito* models. AMPLE needs to know the location of the ROSETTA installation in order to find all of the various tools it needs for creating decoy structures. Please ensure that the path to the ROSETTA top level directory is specified in the interface or with the ``-rosetta_dir`` flag if using a script. 

.. note:: 
   ROSETTA is a comprehensive package and requires compilation from source code. For a detailed explanation about how to install ROSETTA see `ROSETTA installation`_.

------------------------------------------------------------------

Running AMPLE
=============
Input Files
-----------
AMPLE requires a FASTA file and an MTZ file in order to run. There are some other files required, which will be described below.

.. note::
   You can download all the data files `here <https://github.com/rigdenlab/ample-examples/archive/master.zip>`_.

Rosetta Input Files
-------------------
For *ab initio* modelling ROSETTA requires Robetta fragment files (3 and 9 residues), (fragment files can be generated using the `Robetta online server`_. Note that registration is required for this service). For this example these fragment files have already been calculated and part of the downloaded files.

AMPLE Setup
-----------
System-dependent example scripts to run AMPLE are shown below:

UNIX (Linux|Mac)
^^^^^^^^^^^^^^^^

.. literalinclude:: /../examples/toxd-example/run.sh
   :language: bash
   :lines: 15-25

Windows
^^^^^^^

.. literalinclude:: /../examples/toxd-example/run.bat
   :language: batch
   :lines: 3-10

First we set the path to the location where ROSETTA is installed. This is then input into ample using the ``-rosetta_dir`` flag.

Next we need to provide the locations of our input files, this is done using the following flags:

* ``-fasta_input`` – location of the FASTA file.
* ``-mtz_input`` – location of the MTZ file.

*UNIX (Linux|Mac) only:*

* ``-frags_3mers`` – location of the 3 residue fragment from the Robetta server.
* ``-frags_9mers`` – location of the 9 residue fragment from the Robetta server.
* ``-nmodels`` – (optional, default 1000) flag to specify the number of models we want to make with ROSETTA.

.. note::
   In this test case only 30 models are created as we known Rosetta models this structure well . However, in a typical AMPLE run we will generate at least 500 models to increase the likelihood that the correct fold is modelled.

.. note:: 
   Models for AMPLE can also be generated using the `QUARK online server`_. This is particularly useful, if you are running Windows or do not have access to a local Rosetta installation. 
   Go to the website, submit your sequence and download the tarball. Provide the tarball containing your models to AMPLE using the ``-models`` flag.

Finally we can specify some options about how AMPLE will run. Here we use:

* ``-percent`` – flag which specifies by what percentage to truncate the models.
* ``-use_shelxe`` – specifies whether we should use shelxe or not.
* ``-nproc`` – lets you specify how many processors you want to use.

For a full list possible options see :ref:`AMPLE options <cl_options>`.

------------------------------------------------------------------

AMPLE Output
============
On starting a separate window will appear summarising the progress of AMPLE and any results. The window will contain up to four tabs, the contents of which are explained below:

Summary
-------
The summary tab contains different sections. Below you can find information about each:

Ensembling Results
^^^^^^^^^^^^^^^^^^
There is a brief summary of the type of truncation that was undertaken and then a table listing each ensemble. The columns of the table are:

* **Name:** the name of the ensemble. This is used to name the pdb file and the directory where mrbump carries out molecular replacement.
* **Truncation Level:** the percentage of the model remaining after the varying residues were pruned away.
* **Variance Threshold:** AMPLE constructs ensembles by pruning back the most variable residues based on their variance as calculated by THESEUS. The variance threshold is the THESEUS variance score for the most variable residue that remains in this ensemble.
* **No. Residues:** the number of residues for each model in the ensemble.
* **Radius Threshold:** the truncated models are sub-clustered after truncation under 3 different radius thresholds to create the ensemble, and this is the radius used for this sub-cluster.
* **No. Decoys:** the number of models within this ensemble.
* **Number of Atoms:** the number of atoms for each model in the ensemble.
* **Sidechain Treatment:** 

  * *allatom* – all sidechains were retained
  * *reliable* – MET, ASP, PRO, GLN, LYS, ARG, GLU, SER were retained
  * *polyAla* – all sidechains were stripped back to polyalanine

MrBUMP Results
^^^^^^^^^^^^^^
This section displays a table with the results of running MrBUMP on each of the ensembles, for this example you will have information for the following headings.

.. figure:: ../images/summary_toxd.png
   :align: center

* **ensemble_name:** this matches the name from the ensemble section.
* **MR_program:** the program used for Molecular Replacement.
* **Solution type:** the MrBUMP categorisation of the solution

  * *GOOD* - final Rfree <=0.35
  * *MARGINAL* - final Rfree <= 0.48 OR final Rfree <= 0.5 and the ratio between the initial and final Rfree is <= 0.8, OR final Rfree <= 0.55 and the ratio between the initial and final Rfree is <= 0.95
  * *POOR* - anything else
  * *no_job_directory* - a script has been prepared, but the job hasn’t been run yet
  * *unfinished* - the job is running or has stopped without generating any results

* **PHASER_LLG:** the PHASER log-likelihood gain for the Molecular Replacement solution.
* **PHASER_TFZ:** PHASER Translation Function Z-score for the Molecular Replacement solution.
* **REFMAC_Rfact:** Rfact score for REFMAC refinement of the Molecular Replacement solution.
* **REFMAC_Rfree:** Rfree score for REFMAC refinement of the Molecular Replacement solution.
* **SHELXE_CC:** SHELXE Correlation Coefficient score after C-alpha trace.
* **SHELXE_ACL:** Average Chain Length of the fragments of the SHELXE C-alpha trace.

Results
-------
The Results tab displays the final results of AMPLE after running MrBUMP on the ensembles.

.. figure:: ../images/results_toxd.png
   :align: center

The tab is split into two sections. The upper section shows the top three results as ranked by their SHELXE CC score. The lower section shows the top three results as ranked by their PHASER TFZ score. These may or may not be different. Within each section, the left-hand menu displays a list of ensemble names – these match the names from the Ensembles section in the Summary tab. Clicking on any item will display the results for that ensemble in the central pane. At the top is a table that matches the MrBUMP entry from the Summary tab, and there are then sections for the files output by each program run by MrBUMP. The files can either be displayed directly or opened directly with COOT or CCP4MG using the displayed buttons.

Typically a result with a SHELXE CC score of 25 or higher **and** a SHELXE ACL of 10 or higher will indicate a correct solution. 

.. note:: 
   The results you obtain may be slightly different to those presented above as you are generating a new slightly different set of *ab initio* models.


Log File
--------
This displays the text output by AMPLE as it is running. Any problems or errors will be displayed here.

.. figure:: ../images/log_toxd.png
   :align: center

Citations
---------
This section lists the programs and algoriths that are using in the AMPLE job and gives a list of references to be cited should AMPLE find a solution.

.. figure:: ../images/citation_toxd.png
   :align: center

------------------------------------------------------------------


.. _QUARK online server: http://zhanglab.ccmb.med.umich.edu/QUARK
.. _Robetta online server: http://robetta.bakerlab.org/fragmentsubmit.jsp
.. _Rosetta installation: http://ccp4wiki.org/~ccp4wiki/wiki/index.php?title=Installing_Rosetta
