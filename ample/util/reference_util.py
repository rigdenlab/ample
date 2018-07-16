"""Various miscellaneous functions"""

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "16 Jul 2018"

import os

from ample.util.ample_util import is_file
from ample.constants import SHARE_DIR

def construct_references(optd, write_file=True):
    """Construct the reference list

    Description
    -----------
    This is somewhat a very basic BibTex parser. It is under no circumstances foolproof but rather
    the bare minimum for what is required in AMPLE.

    If a more sophisticated parser is required, refer to one of the python packages available online.

    Parameters
    ----------
    optd : dict
       A dictionary containing AMPLE options

    Returns
    -------
    str
        A string containing the references
    """
    # ========================================
    # Get the filename and check we can use it
    # ========================================
    ref_fname = os.path.join(SHARE_DIR, "include", "ample.bib")
    if not is_file(ref_fname):
        msg = "Cannot find BibTex file containing references. " \
              "Please determine them yourself and cite AMPLE."
        return msg

    articles = []
    article = {}
    entry = False
    with open(ref_fname, "r") as fhin:
        for line in fhin.readlines():

            line = line.strip()

            if not line:                                # Make sure line is not empty
                continue

            elif line.startswith("@"):                  # Beginning of all BibTex entry blocks
                entry = True                            # Notify that we have an article block
                unique_id = line.replace("@article{", "").replace(",", "")
                article = {'unique_id': unique_id}      # Reset the article dictionary

            elif line == "}":                           # End of all BibTex entry blocks
                entry = False                           # Notify that we article block is over
                articles.append(article)                # Save the article dictionary

            elif entry:                                 # BibTex entry block
                # Some dirty line handling.
                # Not very bulletproof but should do for now
                line = line.replace("{", "").replace("}", "")
                key, value = [l.strip() for l in line.split("=")]
                value = value.rstrip(",").replace("\"", "")
                # Save the data to the article entry
                article[key] = value

    applications = used_apps(optd)
    to_remove = [i for i, art in enumerate(articles) if art['label'] not in applications]
    # reverse order so that we don't throw off the subsequent indexes
    for index in sorted(to_remove, reverse=True):
        del articles[index]
    # =========================================================================
    # Somewhat a template of how we want to write each article in BibTex format
    # =========================================================================
    template_bib = "@article{{{unique_id},{sep}author = {{{author}}},{sep}doi = {{{doi}}},{sep}" \
                   "journal = {{{journal}}},{sep}number = {{{number}}},{sep}pages = {{{pages}}},{sep}" \
                   "title = {{{{{title}}}}},{sep}volume = {{{volume}}},{sep}year = {{{year}}},{sep}}}{sep}"
    references_bib = [template_bib.format(sep=os.linesep, **a) for a in articles]

    if write_file:
        ref_fname = os.path.join(optd['work_dir'], optd['name']+".bib")
        with open(ref_fname, "w") as fhout:
            fhout.write(os.linesep.join(references_bib))

    # ==========================================================
    # Somewhat a template of how we want to display each article
    # ==========================================================
    for article in articles:
        # Shorten the author list for the displayed message
        article['author'] = article['author'].split(" and ")[0].split(",")[0] + " et al."
        # Display page separator as single dash
        article['pages'] = article['pages'].replace("--", "-")
    template_msg = "* {author} ({year}). {title}. {journal} {volume}({number}), {pages}. [doi:{doi}]"
    references_msg = [template_msg.format(**article) for article in articles]

    return (os.linesep*2).join(references_msg)

def used_apps(optd):
    """determine which applications were used and return their labels"""
    from ample.ensembler.constants import SPICKER_RMSD, SPICKER_TM

    # Papers that are always listed
    labels = ['AMPLE', 'AMPLE_QUARK', 'AMPLE_COILS', 'AMPLE_CONTACTS', 'CCP4', 'GESAMT']

    # This conditional attempts to recognise which programs have run and therefore
    # prints only the relevant references. It is under no circumstance perfect.
    if optd['nmr_model_in']:
        labels += ['AMPLE_NMR']
    if optd['make_models']:
        labels += ['ROSETTA']
    if optd['transmembrane']:
        labels += ['AMPLE_TM']
    if optd['quark_models']:
        labels += ['QUARK']
        
    # Flags related to cluster-and-truncate approach
    elif not optd['import_ensembles']:
        labels += ['THESEUS']
        if optd['use_scwrl']:
            labels += ['SCWRL4']
        elif optd['cluster_method'] in [SPICKER_RMSD, SPICKER_TM]:
            labels = ['SPICKER']
        elif optd['cluster_method'] in ['fast_protein_cluster']:
            labels += ['FPC']

    # Flags related to Molecular Replacement
    elif optd['do_mr']:
        labels += ['MrBUMP']
        if optd['refine_rebuild_arpwarp'] or optd['shelxe_rebuild_arpwarp']:
            labels += ['Arp_Warp']
        elif optd['refine_rebuild_buccaneer'] or optd['shelxe_rebuild_buccaneer']:
            labels += ['Buccaneer']
        elif optd['use_shelxe']:
            labels += ['SHELXE']
        elif 'molrep' in optd['mrbump_programs']:
            labels += ['Molrep', 'REFMAC']
        elif 'phaser' in optd['mrbump_programs']:
            labels += ['Phaser', 'REFMAC']

    return labels

# ======================================================================
# Some default string messages that we need during the program to inform
# the user of certain information
# ======================================================================

header = """
#########################################################################
#########################################################################
#########################################################################
# CCP4: AMPLE - Ab Initio Modelling Molecular Replacement               #
#########################################################################

"""

# ======================================================================
# ======================================================================

reference_str_log = """
#########################################################################

The authors of specific programs should be referenced where applicable:

{refs}

"""

reference_str_gui = """
AMPLE is a pipeline that 

{refs}

"""

# ======================================================================
# ======================================================================

survey_url = "http://goo.gl/forms/7xP9M4P81O"
footer = """

#########################################################################
#***********************************************************************#
#*                          How did we do?                             *#
#*                                                                     *#
#*           Please follow this link and leave some feedback!          *#
#*                                                                     *#
#*                 {url}                      *#
#***********************************************************************#
#########################################################################
""".format(url=survey_url)
# ======================================================================
# ======================================================================

