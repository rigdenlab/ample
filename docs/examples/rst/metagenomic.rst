.. _example_metagenomic:

===================================
Using models derived from databases
===================================

Covariance-assisted *ab initio* models representing structurally uncharacterised Pfam families are now available on a large scale in databases such as the `GREMLIN database <https://gremlin2.bakerlab.org/meta_struct.php>`_ and the `PconsFam database <http://pconsfam.bioinfo.se/>`_. If models are available for your target protein, these models can be used in AMPLE.

------------------------------------------------------------------

Running AMPLE
=============
Input Files
-----------
AMPLE requires a FASTA file and an MTZ file in order to run. There are some other files required, which will be described below.

.. note::
   You can download all the data files `here <hhttps://github.com/rigdenlab/ample-examples/archive/master.zip>`_.

AMPLE Setup
-----------
System-dependent example scripts to run AMPLE are shown below:

UNIX (Linux|Mac)
^^^^^^^^^^^^^^^^

.. literalinclude:: /../examples/metagenomic/run.sh
   :language: bash
   :lines: 12-18

Windows
^^^^^^^

.. literalinclude:: /../examples/metagenomic-helices/run.bat
   :language: batch
   :lines: 6-13

First we provide the location of the input files. This is done with the following flags:

* ``-fasta_input`` – location of the FASTA file.
* ``-mtz_input`` – location of the MTZ file.

Next we provide the flag that tells AMPLE to run in ideal helices mode:

* ``-models`` - specifies the path to the directory of the covariance-assisted *ab initio* models

Finally we provide some options about how AMPLE will run:

* ``-nproc`` – lets you specify how many processors you want to use.
* ``-show_gui`` - Flag to display the AMPLE gui. This is set by default when running through CCP4i or CCP4i2 but must be manually set on the command line to generate the output GUI shown below.

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

.. figure:: ../images/summary_ensembling_metagenomic.png
   :align: center

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

.. figure:: ../images/summary_mrbump_metagenomic.png
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


Results
-------
The Results tab displays the final results of AMPLE after running MrBUMP on the ensembles.

.. figure:: ../images/results_metagenomic.png
   :align: center

The tab is split into two sections. The upper section shows the top three results as ranked by their SHELXE CC score. The lower section shows the top three results as ranked by their PHASER TFZ score. These may or may not be different. Within each section, the left-hand menu displays a list of ensemble names. Clicking on any item will display the results for that ensemble in the central pane. At the top is a table that matches the MrBUMP entry from the Summary tab, and there are then sections for the files output by each program run by MrBUMP. The files can either be displayed directly or opened directly with COOT or CCP4MG using the displayed buttons.

Typically a result with a SHELXE CC score of 25 or higher **and** a SHELXE ACL of 10 or higher will indicate a correct solution.

.. note::
   The results you obtain may be slightly different to those presented above as you are generating a new slightly different set of *ab initio* models.


Log File
--------
This displays the text output by AMPLE as it is running. Any problems or errors will be displayed here.

.. figure:: ../images/log_metagenomic.png
   :align: center

Citations
---------
This section lists the programs and algoriths that are using in the AMPLE job and gives a list of references to be cited should AMPLE find a solution.

.. figure:: ../images/citation_metagenomic.png
   :align: center