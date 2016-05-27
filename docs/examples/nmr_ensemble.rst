.. _nmr_ensemble:

=====================
Using an NMR ensemble
=====================

AMPLE can attempt molecular replacement with ensembles created from NMR ensembles. In the simplest case, AMPLE will split an NMR ensemble into its constituent models and carry out its standard truncation/clustering algorithm to generate the ensembles. When used in this way, ROSETTA is not required to be installed. 

.. note::
   To use AMPLE in this way, just supply the NMR model with the -nmr_model_in flag and set the -nmr_remodel flag to False.

Running AMPLE
=============
Input Files
-----------
AMPLE requires a FASTA file and an MTZ file in order to run. For molecular replacement with ensembles created from NMR ensembles you must also supply the NMR model. 

.. note::
   For this example the NMR model is available `here`_

AMPLE Setup
-----------
System-dependent example scripts to run AMPLE are shown below:

UNIX (Linux|Mac)
^^^^^^^^^^^^^^^^

.. code-block:: bash

   #!/bin/bash

   # NMR ensembling example – 2LC9 is an ensemble model of a minor and transiently formed state of
   # a T4 lysozyme mutant. The target 102l is X-ray data

   $CCP4/bin/ample \
   -mtz input/102l.mtz \
   -fasta input/102L.fasta \
   -name 102l \
   -nmr_model_in input/2LC9.pdb \
   -nmr_remodel False \
   -quick_mode True

Windows
^^^^^^^

.. code-block:: batch

   REM NMR ensembling example – 2LC9 is an ensemble model of a minor and transiently formed state of
   REM a T4 lysozyme mutant. The target 102l is X-ray data

   %CCP4%\bin\ample.bat ^
   -mtz input\102l.mtz ^
   -fasta input\102L.fasta ^
   -name 102l ^
   -nmr_model_in input\2LC9.pdb ^
   -nmr_remodel False ^
   -quick_mode True

We need to provide the locations of our input files, this is done using the following flag:

* ``-fasta_input`` – for our FASTA file.
* ``-mtz_input``  – for our MTZ file.
* ``-nmr_model``  – for our NMR model.

Next we can specify a few run options for AMPLE:

* ``-name`` – specifies the job name.
* ``-nmr_remodel`` – specifies whether to remodel the NMR structures.
* ``-quick_mode`` – Preset options to run quickly, but less thoroughly.

For a full list of options see :ref:`AMPLE options <ample_options>`.

AMPLE Output
============
On starting a separate window will appear summarising the progress of AMPLE and any results. The window will contain up to three tabs, the contents of which are explained below:

Summary tab
-----------
This is divided into two sections that display a summary of the results of the ensembling process and the results of the Molecular Replacement with MrBUMP respectively:

.. image:: ../_static/summary_nmr.png

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

.. image:: ../_static/results_nmr.png

The tab is split into two sections. The upper section shows the top three results as ranked by their SHELXE CC score. The lower section shows the top three results as ranked by their PHASER TFZ score. These may or may not be different. Within each section, the left-hand menu displays a list of ensemble names – these match the names from the Ensembles section in the Summary tab. Clicking on any item will display the results for that ensemble in the central pane. At the top is a table that matches the MrBUMP entry from the Summary tab, and there are then sections for the files output by each program run by MrBUMP. The files can either be displayed directly or opened directly with COOT or CCP4MG using the displayed buttons.

Typically a result with a SHELXE CC score of 25 or higher **and** a SHELXE ACL of 10 or higher will indicate a correct solution..


Log File
--------
This displays the text output by AMPLE as it is running. Any problems or errors will be displayed here.

.. image:: ../_static/log_nmr.png





.. _here: https://drive.google.com/file/d/0B3NdI1poe0RhSVFyRjRHSER1Y0k/view.



