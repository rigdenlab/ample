"""Various miscellaneous functions"""

__author__ = "Jens Thomas, and Felix Simkovic"
__date__ = "16 Jul 2018"

from collections import OrderedDict
import copy
from enum import Enum
import os

from ample.util.ample_util import is_file
from ample.constants import SHARE_DIR
from ample.ensembler.constants import SPICKER_RMSD, SPICKER_TM

 
class ReferenceManager():
     
    # Section Names
    class SECTIONS(Enum):
        __order__ = 'GENERAL MODELLING MODEL_PREP MR REFINEMENT DM AUTOBUILD'
        GENERAL = 'General'
        MODELLING = 'Modelling'
        MODEL_PREP = 'Search model preparation'
        MR = 'Molecular Replacement'
        REFINEMENT = 'Refinement'
        DM = 'Main-chain tracing and density modification'
        AUTOBUILD = 'Autobuilding'
        
    SEC_TAG = 'h3'
     
    def __init__(self, optd):
        self.references = {}
        self.ordered_labels = []
        self.section_labels = OrderedDict()
        self.setup_references(optd)
        self.setup_sections(optd)
        
    def setup_references(self, optd):
        ref_fname = os.path.join(SHARE_DIR, "include", "ample.bib")
        if not is_file(ref_fname):
            msg = "Cannot find BibTex file containing references. " \
                  "Please determine them yourself and cite AMPLE."
            return msg
        article = {}
        entry = False
        with open(ref_fname, "r") as fhin:
            for line in fhin.readlines():
                line = line.strip()
                if not line:
                    continue
                elif line.startswith("@"):
                    # Beginning of all BibTex entry blocks                  
                    entry = True
                    unique_id = line.replace("@article{", "").replace(",", "")
                    article = {'unique_id': unique_id}      # Reset the article dictionary
                elif line == "}":
                    # End of all BibTex entry blocks                      
                    entry = False
                    self.references[article['label']] = article
                elif entry:
                    # BibTex entry block                               
                    # Some dirty line handling.
                    # Not very bulletproof but should do for now
                    line = line.replace("{", "").replace("}", "")
                    key, value = [l.strip() for l in line.split("=")]
                    value = value.rstrip(",").replace("\"", "")
                    # Save the data to the article entry
                    article[key] = value
        return


    def setup_sections(self, optd):
        def key_set(key):
            return key in optd and optd[key]
        # Build up list of program reference labels, ordered by setions
        for section in self.SECTIONS:
            if section == self.SECTIONS.GENERAL:
                self.section_labels[self.SECTIONS.GENERAL] = ['AMPLE', 'AMPLE_QUARK', 'AMPLE_COILS', 'AMPLE_CONTACTS', 'CCP4']
            elif section == self.SECTIONS.MODELLING:
                labels = []
                if key_set('nmr_model_in'):
                    labels.append('AMPLE_NMR')
                if key_set('make_models'):
                    labels.append('ROSETTA')
                if key_set('transmembrane'):
                    labels.append('AMPLE_TM')
                if key_set('quark_models'):
                    labels.append('QUARK')
                self.section_labels[self.SECTIONS.MODELLING] = labels
            elif section == self.SECTIONS.MODEL_PREP:
                labels = []
                if not key_set('import_ensembles'):
                    labels.append('THESEUS')
                    if key_set('use_scwrl'):
                        labels.append('SCWRL4')
                    elif optd['cluster_method'] in [SPICKER_RMSD, SPICKER_TM]:
                        labels.append('SPICKER')
                    elif optd['cluster_method'] in ['fast_protein_cluster']:
                        labels.append('FPC')
                self.section_labels[self.SECTIONS.MODEL_PREP] = labels
            if optd['do_mr']:
                if section == self.SECTIONS.MR:
                    labels = ['MRBUMP']
                    if 'molrep' in optd['mrbump_programs']:
                        labels.append('MOLREP')
                    if 'phaser' in optd['mrbump_programs']:
                        labels.append('PHASER')
                    self.section_labels[self.SECTIONS.MR] = labels
                elif section == self.SECTIONS.REFINEMENT:
                    self.section_labels[self.SECTIONS.REFINEMENT] = ['REFMAC']
                elif section == self.SECTIONS.DM:
                    self.section_labels[self.SECTIONS.DM] = []
                    if key_set('use_shelxe'):
                        self.section_labels[self.SECTIONS.DM].append('SHELXE')
                elif section == self.SECTIONS.AUTOBUILD:
                    labels = []
                    if key_set('refine_rebuild_arpwarp') or key_set('shelxe_rebuild_arpwarp'):
                        labels += ['ARPWARP']
                    elif key_set('refine_rebuild_buccaneer') or key_set('shelxe_rebuild_buccaneer'):
                        labels += ['BUCCANEER']        
                    self.section_labels[self.SECTIONS.AUTOBUILD] = labels
                
        # Generate ordered list of all relevant reference labels
        for section in self.SECTIONS:
            if section in self.section_labels:
                self.ordered_labels += self.section_labels[section]
        return
        
    def methods_as_html(self):
        html = "<p>The following lists the different programs that were used and general references that should be cited." + \
               "Numbers in superscript next to program names refer to the number of the program reference in the overall list of references.</p>"
        for section in self.SECTIONS:
            if section == self.SECTIONS.GENERAL:
                html += '<p>The first {} references are general AMPLE references and should be cited in all cases</p>'.format(len(self.section_labels[self.SECTIONS.GENERAL]))
            elif section == self.SECTIONS.MODELLING and len(self.section_labels[self.SECTIONS.MODELLING]):
                html += '<{}>{}</{}>'.format(self.SEC_TAG, self.SECTIONS.MODELLING.value, self.SEC_TAG)
                html += "<p>The following programs were used for model building:</p>"
                html += '<ul>'
                for label in self.section_labels[self.SECTIONS.MODELLING]:
                    html += "<li>{}<sup>{}</sup></li>".format(label, self.ordered_labels.index(label) + 1)
                html += "</ul>"
            elif section == self.SECTIONS.MODEL_PREP and len(self.section_labels[self.SECTIONS.MODEL_PREP]):
                html += '<{}>{}</{}>'.format(self.SEC_TAG, self.SECTIONS.MODEL_PREP.value, self.SEC_TAG)
                html += '<p>Model analysis and search model preparation was carried out with the following programs:</p>'
                html += '<ul>'
                for label in self.section_labels[self.SECTIONS.MODEL_PREP]:
                    html += "<li>{}<sup>{}</sup></li>".format(label, self.ordered_labels.index(label) + 1)
                html += "</ul>"
            elif section == self.SECTIONS.MR and len(self.section_labels[self.SECTIONS.MR]):
                html += '<{}>{}</{}>'.format(self.SEC_TAG, self.SECTIONS.MR.value, self.SEC_TAG)
                html += "<p>Molecular Replacement was undertaken with the following programs:</p>"
                html += '<ul>'
                for label in self.section_labels[self.SECTIONS.MR]:
                    html += "<li>{}<sup>{}</sup></li>".format(label, self.ordered_labels.index(label) + 1)
                html += "</ul>"
            elif section == self.SECTIONS.DM and len(self.section_labels[self.SECTIONS.DM]):
                html += '<{}>{}</{}>'.format(self.SEC_TAG, self.SECTIONS.DM.value, self.SEC_TAG)
                html += "<p>Density modification and main-chain tracing was carried out with the following programs:</p>"
                html += '<ul>'
                for label in self.section_labels[self.SECTIONS.DM]:
                    html += "<li>{}<sup>{}</sup></li>".format(label, self.ordered_labels.index(label) + 1)
                html += "</ul>"
            elif section == self.SECTIONS.AUTOBUILD and len(self.section_labels[self.SECTIONS.AUTOBUILD]):
                html += '<{}>{}</{}>'.format(self.SEC_TAG, self.SECTIONS.AUTOBUILD.value, self.SEC_TAG)
                html += "<p>Autobuilding of the final structure was carried out with the following programs:</p>"
                html += '<ul>'
                for label in self.section_labels[self.SECTIONS.AUTOBUILD]:
                    html += "<li>{}<sup>{}</sup></li>".format(label, self.ordered_labels.index(label) + 1)
                html += "</ul>"
        return html
    
    def references_as_html(self):
        html = '<{}>References</{}>'.format(self.SEC_TAG, self.SEC_TAG)
        html += '<ol>'
        template_txt = "<li> {author} ({year}). {title}. {journal} {volume}({number}), {pages}. [doi:{doi}]</li>"
        for label in self.ordered_labels:
            ref = copy.copy(self.references[label])
            ref['author'] = ref['author'].split(" and ")[0].split(",")[0] + " et al."
            ref['pages'] = ref['pages'].replace("--", "-")
            html += template_txt.format(**ref)
        html += '</ol>'
        return html
    
    def references_as_text(self):
        template_txt = "* {author} ({year}). {title}. {journal} {volume}({number}), {pages}. [doi:{doi}]"
        text = ""
        for label in self.ordered_labels:
            ref = copy.copy(self.references[label])
            ref['author'] = ref['author'].split(" and ")[0].split(",")[0] + " et al."
            ref['pages'] = ref['pages'].replace("--", "-")
            text += template_txt.format(**ref) + os.linesep*2
        return text
    
    def save_references_to_file(self, optd=None):
        # =========================================================================
        # Somewhat a template of how we want to write each article in BibTex format
        # =========================================================================
        template_bib = "@article{{{unique_id},{sep}author = {{{author}}},{sep}doi = {{{doi}}},{sep}" \
                       "journal = {{{journal}}},{sep}number = {{{number}}},{sep}pages = {{{pages}}},{sep}" \
                       "title = {{{{{title}}}}},{sep}volume = {{{volume}}},{sep}year = {{{year}}},{sep}}}{sep}"
        references_bib = [template_bib.format(sep=os.linesep, **self.references[l]) for l in self.ordered_labels]
        #ref_fname = os.path.join(optd['work_dir'], optd['name']+".bib")
        ref_fname = 'foo.txt'
        with open(ref_fname, "w") as fhout:
            fhout.write(os.linesep.join(references_bib))
        return ref_fname

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

