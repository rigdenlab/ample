Changelog
=========

1.5.0
------

Added
~~~~~
- logging can be configured with a JSON file.
- modelling can be run as a standalone module.
- added a Coiled-coil mode.
- added the ability to model multimers (developed with Owen Davies at Newcastle).
- ``ample_into_ccp4.sh`` script for linking AMPLE into the current CCP4 installation.

Changed
~~~~~~~
- updated the processing of supplied models so multiple chains are supported.
- automatic determination of the type of models supplied (homologs, NMR ensemble, etc.)
- boolean arguments (e.g. -ideal_helices) now don't require an argument to be True.
- removed all MAXCLUSTER code.
- ROSETTA modelling now runs using a flagsfile.
- turned off all rebuilding for resolutions > 4.0
- updated integration tests to be less fragile.



1.4.5
------

Added
~~~~~
- code should be now both Python 2 and 3 compliant
- added citation tab to jscofe gui and updated citation tracking and references.
- '-percent_fixed_intervals' option for truncating so that we can specify a list of percentage intervals for truncating models.
- '-purge' option now accepts a level of 0, 1 or 2 to select how much data is purged.

Changed
~~~~~~~
- SHELXE is only run depending on the xtal data resolution (as specified by mrbump_util.SHELXE_MAX_PERMITTED_RESOLUTION - currently 3.0A)
- don't exit on error code from MRBUMP job submission - just warn about errors.
- numerous small improvements and bug fixes
