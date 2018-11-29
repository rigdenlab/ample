Changelog
=========

1.4.5
------

Added
~~~~~
- code should be now both Python 2 and 3 compliant
- '-percent_fixed_intervals' option for truncating so that we can specify a list of percentage intervals for truncating models.
- '-purge' option now accepts a level of 0, 1 or 2 to select how much data is purged.

Changed
~~~~~~~
- SHELXE is only run depending on the xtal data resolution (as specified by mrbump_util.SHELXE_MAX_PERMITTED_RESOLUTION - currently 3.0A)
- don't exit on error code from MRBUMP job submission - just warn about errors.
- numerous small improvements and bug fixes
