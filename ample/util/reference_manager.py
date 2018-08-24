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
        self.citation_file_path = None
        self.section_labels = OrderedDict()
        self.setup_references()
        self.setup_sections(optd)
        
    def setup_references(self):
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
        # Create default lists
        for section in self.SECTIONS:
            self.section_labels[section] = []
        
        # Build up list of program reference labels, ordered by sections
        for section in self.SECTIONS:
            if section == self.SECTIONS.GENERAL:
                self.section_labels[section] = ['AMPLE', 'CCP4', 'AMPLE_COILED-COILS', 'AMPLE_CONTACTS']
            elif section == self.SECTIONS.MODELLING:
                labels = []
                if optd.get('ideal_helices'):
                    labels.append('AMPLE_COILED-COILS')
                if optd.get('make_models'):
                    labels.append('ROSETTA')
                if optd.get('nmr_model_in'):
                    labels.append('AMPLE_NMR')
                if optd.get('quark_models'):
                    labels.append('QUARK')
                    labels.append('AMPLE_QUARK')
                if optd.get('transmembrane'):
                    labels.append('AMPLE_TRANSMEMBRANE')
                self.section_labels[section] = labels
            elif section == self.SECTIONS.MODEL_PREP:
                labels = []
                if not optd.get('import_ensembles'):
                    labels += ['CCTBX', 'THESEUS', 'GESAMT']
                    if optd.get('use_scwrl'):
                        labels.append('SCWRL4')
                    elif optd['cluster_method'] in [SPICKER_RMSD, SPICKER_TM]:
                        labels.append('SPICKER')
                    elif optd['cluster_method'] in ['fast_protein_cluster']:
                        labels.append('FPC')
                    self.section_labels[section] = labels
            if optd.get('do_mr'):
                if section == self.SECTIONS.MR:
                    labels = ['MRBUMP']
                    if optd.get('mrbump_programs'):
                        if 'molrep' in optd['mrbump_programs']:
                            labels.append('MOLREP')
                        if 'phaser' in optd['mrbump_programs']:
                            labels.append('PHASER')
                    self.section_labels[section] = labels
                elif section == self.SECTIONS.REFINEMENT:
                    self.section_labels[self.SECTIONS.REFINEMENT] = ['REFMAC']
                elif section == self.SECTIONS.DM:
                    if optd.get('use_shelxe'):
                        self.section_labels[self.SECTIONS.DM].append('SHELXE')
                elif section == self.SECTIONS.AUTOBUILD:
                    labels = []
                    if optd.get('refine_rebuild_arpwarp') or optd.get('shelxe_rebuild_arpwarp'):
                        labels += ['ARPWARP']
                    elif optd.get('refine_rebuild_buccaneer') or optd.get('shelxe_rebuild_buccaneer'):
                        labels += ['BUCCANEER']        
                    self.section_labels[section] = labels
                
        # Generate ordered list of all relevant reference labels
        for section in self.SECTIONS:
            if section in self.section_labels:
                self.ordered_labels += self.section_labels[section]
        return
    
    @property
    def methods_as_html(self):
        html = "<p>This section lists the programs and algorithms that were used in this job and the references that should be cited. " + \
               "Numbers in superscript next to program/reference names refer to the number of the program reference in the overall list of references.</p>"
        for section in self.SECTIONS:
            if section == self.SECTIONS.GENERAL:
                html += '<p>The first 2 references should be cited in all cases.</p>' + \
                        '<p>If your protein was a coiled-coil protein, please cite reference number 3.</p>' + \
                        '<p>If coevolutionary contact information was used in the generation of your models, please cite reference number 4.</p>'
            elif section == self.SECTIONS.MODELLING and len(self.section_labels[self.SECTIONS.MODELLING]):
                standfirst = "<p>The following programs or algorithims were used for model building:</p>"
                html += self._methods_section_html(self.SECTIONS.MODELLING, standfirst)
            elif section == self.SECTIONS.MODEL_PREP and len(self.section_labels[self.SECTIONS.MODEL_PREP]):
                standfirst ='<p>Model analysis and search model preparation was carried out with the following programs:</p>'
                html += self._methods_section_html(self.SECTIONS.MODEL_PREP, standfirst)
            elif section == self.SECTIONS.MR and len(self.section_labels[self.SECTIONS.MR]):
                standfirst ='<p>Molecular Replacement was carried out with the following programs:</p>'
                html += self._methods_section_html(self.SECTIONS.MR, standfirst)
            elif section == self.SECTIONS.REFINEMENT and len(self.section_labels[self.SECTIONS.REFINEMENT]):
                standfirst ='<pRefinement of the MR solutions carried out with the following programs:</p>'
                html += self._methods_section_html(self.SECTIONS.REFINEMENT, standfirst)
            elif section == self.SECTIONS.DM and len(self.section_labels[self.SECTIONS.DM]):
                standfirst ='<p>Density modification and main-chain tracing was carried out with the following programs:</p>'
                html += self._methods_section_html(self.SECTIONS.DM, standfirst)
            elif section == self.SECTIONS.AUTOBUILD and len(self.section_labels[self.SECTIONS.AUTOBUILD]):
                standfirst ='Autobuilding of the final structure was carried out with the following programs:</p>'
                html += self._methods_section_html(self.SECTIONS.AUTOBUILD, standfirst)
        return html
    
    def _methods_section_html(self, section, standfirst):
        mysec = self.SECTIONS(section)
        html = '<{}>{}</{}>'.format(self.SEC_TAG, mysec.value, self.SEC_TAG)
        html += standfirst
        html += '<ul>'
        for label in self.section_labels[mysec]:
            html += "<li>{}<sup>{}</sup></li>".format(label, self.ordered_labels.index(label) + 1)
        html += "</ul>"
        return html

    @property
    def citations_as_html(self):
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
    
    @property
    def citations_as_text(self):
        txt = """A number of programs and algorithms were used within the this run of AMPLE.

The following is a list of citations for this run:

{0}
""".format(self.citation_list_as_text)
        if self.citation_file_path:
            txt += """
A bibtex file with these references has been saved to the following file:

{0}

""".format(self.citation_file_path)
            return txt
    
    @property
    def citation_list_as_text(self):
        template_txt = "* {author} ({year}). {title}. {journal} {volume}({number}), {pages}. [doi:{doi}]"
        text = ""
        for label in self.ordered_labels:
            ref = copy.copy(self.references[label])
            ref['author'] = ref['author'].split(" and ")[0].split(",")[0] + " et al."
            ref['pages'] = ref['pages'].replace("--", "-")
            text += template_txt.format(**ref) + os.linesep*2
        return text
    
    def save_citations_to_file(self, optd):
        # =========================================================================
        # Somewhat a template of how we want to write each article in BibTex format
        # =========================================================================
        template_bib = "@article{{{unique_id},{sep}author = {{{author}}},{sep}doi = {{{doi}}},{sep}" \
                       "journal = {{{journal}}},{sep}number = {{{number}}},{sep}pages = {{{pages}}},{sep}" \
                       "title = {{{{{title}}}}},{sep}volume = {{{volume}}},{sep}year = {{{year}}},{sep}}}{sep}"
        references_bib = [template_bib.format(sep=os.linesep, **self.references[l]) for l in self.ordered_labels]
        ref_fname = os.path.join(optd['work_dir'], optd['name']+".bib")
        with open(ref_fname, "w") as fhout:
            fhout.write(os.linesep.join(references_bib))
        self.citation_file_path = ref_fname
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

