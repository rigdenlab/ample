#!/usr/bin/env ccp4-python
__author__ = "Jens Thomas"
__date__ = "03 Mar 2015"
__version__ = "1.0"

import logging
import json
import os
import subprocess
import traceback
import urlparse
import uuid

from ample import ensembler
from ample.util import mrbump_util
from ample.util import reference_manager

try:
    import pyrvapi
except:
    pyrvapi = None

# Hack to use Andre's pyrvapi API
try:
    import pyrvapi_ext as API
except ImportError:
    API = None

logger = logging.getLogger(__name__)

# Utility functions
def have_files(results_dict, *keylist, **kwargs):
    """Check if files in dictionary exist.
    
    if kwarg check is 'all' (default) return True if all exist,
    if check is 'any' return True of any exist
    """
    check = 'all'
    if 'check' in kwargs:
        assert kwargs['check'] in ['all', 'any']
        check = kwargs['check']
    found = 0
    for k in keylist:
        if k in results_dict and os.path.isfile(str(results_dict[k])):
            if check == 'any':
                return True
            found += 1
    if check == 'all' and len(keylist) == found:
        return True
    return False


class AmpleOutput(object):
    """Display the output of an AMPLE job."""

    _ensemble_tooltips = {
        "Name": "Ensemble name - used to name the pdb file and the directory where mrbump carries out molecular replacement.",
        "Cluster": "The SPICKER cluster that this ensemble was derived from.",
        "Truncation Level": "Percentage of the model remaining after the varying residues were pruned away",
        "Variance Threshold (A^2)": "THESEUS variance score for the most variable residue that remains in this ensemble",
        "No. Residues": "Number of residues for each model in the ensemble",
        "Radius Threshold (A)": "Radius threshold (1,2 or 3 A) used for subclustering the models in a truncation level",
        "No. Decoys": "Number of models within this ensemble",
        "Number of Atoms": "Number of atoms for each model in the ensemble",
        "Sidechain Treatment": "allatom - all sidechains were retained, reliable - MET, ASP, PRO, GLN, LYS, ARG, GLU, SER were retained, polyAla - all sidechains were stripped back to polyalanine",
    }

    _mrbump_tooltips = {
        "ensemble_name": "The identifier of the AMPLE ensemble search model",
        "MR_program": "Molecular replacement program",
        "Solution_Type": "MRBUMP categorisation of the solution",
        "PHASER_LLG": "PHASER Log-likelihood gain for the Molecular Replacement solution",
        "PHASER_TFZ": "PHASER Translation Function Z-score for the Molecular Replacement solution",
        "REFMAC_Rfact": "Rfact score for REFMAC refinement of the Molecular Replacement solution",
        "REFMAC_Rfree": "Rfree score for REFMAC refinement of the Molecular Replacement solution",
        "BUCC_final_Rfact": "Rfact score for BUCCANEER rebuild of the Molecular Replacement solution",
        "BUCC_final_Rfree": "Rfree score for BUCCANEER rebuild of the Molecular Replacement solution",
        "ARP_final_Rfact": "Rfact score for ARPWARP rebuild of the Molecular Replacement solution",
        "ARP_final_Rfree": "Rfree score for ARPWARP rebuild of the Molecular Replacement solution",
        "SHELXE_CC": "SHELXE Correlation Coefficient score after C-alpha trace",
        "SHELXE_ACL": "Average Chain Length of the fragments of the SHELXE C-alpha trace",
        "SXRBUCC_final_Rfact": "Rfact score for BUCCANEER rebuild of the SHELXE C-alpha trace",
        "SXRBUCC_final_Rfree": "Rfree score for BUCCANEER rebuild of the SHELXE C-alpha trace",
        "SXRARP_final_Rfact": "Rfact score for ARPWARP rebuild of the SHELXE C-alpha trace",
        "SXRAP_final_Rfree": "Rfree score for ARPWARP rebuild of the SHELXE C-alpha trace",
    }

    def __init__(self, amopt):
        self.header = False
        self.jsrview_dir = None
        self.citation_tab_id = None
        self.log_tab_id = None
        self.old_mrbump_results = None
        self.results_tab_id = None
        self.results_tab_sections = []
        self.summary_tab_id = None
        self.summary_tab_ensemble_sec_id = None
        self.summary_tab_results_sec_id = None
        self.summary_tab_survey_sec_id = None
        self.summary_tab_results_sec_table_id = None
        self.summary_tab_pending_sec_id = None

        # Process variables from amopt
        ccp4i2_xml = amopt.get('ccp4i2_xml')
        rvapi_document = amopt.get('rvapi_document')
        work_dir = amopt['work_dir']
        run_dir = amopt['run_dir']
        show_gui = bool(amopt['show_gui'])
        webserver_uri = amopt['webserver_uri']

        # Process modes and set flags
        self.ccp4i2 = bool(ccp4i2_xml)  # Indicate we are running under CCP4I2
        self.jscofe = bool(rvapi_document)  # Indicate we are running under jscofe
        # Show the gui if we are under ccp4i2, jscofe or show_gui has been specified (e.g. ccp4i)
        self.generate_output = self.ccp4i2 | self.jscofe | show_gui
        # No log tab with jscofe or ccp4i2
        self.own_log_tab = not (self.ccp4i2 or self.jscofe)

        # For running under old CCP4online
        if webserver_uri:
            # don't start browser and setup variables for the path on the webserver
            self._webserver_start = len(run_dir) + 1
            self.webserver_uri = webserver_uri
        else:
            self.webserver_uri = None
            self.wbeserver_start = None

        self.setup(work_dir=work_dir, ccp4i2_xml=ccp4i2_xml, rvapi_document=rvapi_document, show_gui=show_gui)
        return

    def setup(self, work_dir=None, ccp4i2_xml=None, rvapi_document=None, show_gui=False):
        if not pyrvapi or not self.generate_output:
            return
        title = "AMPLE Results"
        logger.debug("Using Andre's Pyrvapi" if API else "COULD NOT FIND Andre's API!")
        if rvapi_document:
            logger.debug("Restoring document: %s", rvapi_document)
            pyrvapi.rvapi_restore_document2(rvapi_document)
            self.jsrview_dir = os.path.dirname(rvapi_document)
        else:
            # Quick hack to init with Andre's stuff - can switch out for Felix's API when done
            logger.debug("Starting with xml %s", ccp4i2_xml)
            self.jsrview_dir = os.path.join(work_dir, "jsrview")
            if not os.path.isdir(self.jsrview_dir):
                os.mkdir(self.jsrview_dir)
            kwargs = dict(
                wintitle=title,
                reportdir=self.jsrview_dir,
                xml=ccp4i2_xml,
                abspaths=False,
                # bug in jsrview:
                # layout = 4 if i1 else 7,
            )
            API.document.newdoc(**kwargs)
        if not self.own_log_tab:
            self.log_tab_id = pyrvapi.rvapi_get_meta()
        if show_gui:
            # We start our own browser
            jsrview = os.path.join(os.environ["CCP4"], "libexec", "jsrview")
            subprocess.Popen([jsrview, os.path.join(self.jsrview_dir, "index.html")])
        return

    def create_citation_tab(self, ample_dict):
        if self.citation_tab_id:
            return
        self.citation_tab_id = "citation_tab"
        pyrvapi.rvapi_insert_tab(self.citation_tab_id, "Citation", self.log_tab_id, False)
        refMgr = reference_manager.ReferenceManager(ample_dict)
        bibtex_file = refMgr.save_citations_to_file(ample_dict)
        if self.ccp4i2:
            # The horror of ccp4i2 means that this all gets dumped into xml so we can't use any markup tags
            tdata = refMgr.citations_as_text
        else:
            tdata = refMgr.methods_as_html
            tdata += refMgr.citations_as_html
            tdata += '<hr><p>A bibtex file with the relevant citations has been saved to: {}</p>'.format(bibtex_file)
        pyrvapi.rvapi_add_text(tdata, self.citation_tab_id, 0, 0, 1, 1)
        if not self.ccp4i2:
            pyrvapi.rvapi_add_data(
                "bibtex_file",
                "Citations as BIBTEX",
                self.fix_path(bibtex_file),
                "text",
                self.citation_tab_id,
                2,
                0,
                1,
                1,
                True,
            )
        return self.citation_tab_id

    def create_log_tab(self, ample_dict):
        if not self.own_log_tab or self.log_tab_id:
            return
        logfile = ample_dict['ample_log']
        if not os.path.isfile(logfile):
            return False
        self.log_tab_id = "log_tab"
        logurl = self.fix_path(logfile)
        pyrvapi.rvapi_add_tab(self.log_tab_id, "Log file", True)
        pyrvapi.rvapi_append_content(logurl, True, self.log_tab_id)
        return self.log_tab_id

    def create_results_tab(self, ample_dict):
        if self.ccp4i2 or not self.summary_tab_id or not self._got_mrbump_results(ample_dict):
            return
        mrb_results = ample_dict.get('mrbump_results')
        if mrb_results == self.old_mrbump_results:
            return
        self.old_mrbump_results = mrb_results
        if not self.results_tab_id:
            self.results_tab_id = "results_tab"
            pyrvapi.rvapi_insert_tab(self.results_tab_id, "Results", self.summary_tab_id, False)
        # Delete old sections:
        pyrvapi.rvapi_flush()
        for section_id in self.results_tab_sections:
            pyrvapi.rvapi_remove_widget(section_id)
        pyrvapi.rvapi_flush()
        self.results_tab_sections = []
        ensemble_results = ample_dict['ensembles_data'] if 'ensembles_data' in ample_dict['ensembles_data'] else None
        mrbsum = mrbump_util.ResultsSummary(results=mrb_results[0 : min(len(mrb_results), mrbump_util.TOP_KEEP)])
        mrbsum.sortResults(prioritise="SHELXE_CC")
        self.results_section(
            self.results_tab_id, mrbsum.results, ensemble_results, "Top {0} SHELXE Results".format(mrbump_util.TOP_KEEP)
        )
        mrbsum.sortResults(prioritise="PHASER_TFZ")
        # Add seperator between results - doesn't work as not deleted on refresh
        # pyrvapi.rvapi_add_text("<br/><hr/><br/>", self.results_tab_id, 0, 0, 1, 1)
        self.results_section(
            self.results_tab_id, mrbsum.results, ensemble_results, "Top {0} PHASER Results".format(mrbump_util.TOP_KEEP)
        )
        return self.results_tab_id

    def results_section(self, results_tab_id, mrb_results, ensemble_results, section_title):
        """Results Tab"""
        if not mrb_results:
            return
        # Create unique identifier for this section by using the id
        # All ids will have this appended to avoid clashes
        uid = str(uuid.uuid4())
        section_id = section_title.replace(" ", "_") + uid
        self.results_tab_sections.append(section_id)  # Add to list so we can remove if we update
        pyrvapi.rvapi_add_panel(section_id, results_tab_id, 0, 0, 1, 1)
        pyrvapi.rvapi_add_text("<h3>{0}</h3>".format(section_title), section_id, 0, 0, 1, 1)
        results_tree = "results_tree" + section_id
        pyrvapi.rvapi_add_tree_widget(results_tree, section_title, section_id, 0, 0, 1, 1)
        for r in mrb_results:
            ensemble_name = r['ensemble_name']
            container_id = "sec_{0}".format(ensemble_name) + uid
            pyrvapi.rvapi_add_panel(container_id, results_tree, 0, 0, 1, 1)

            header = "<h3>Results for ensemble: {0}</h3>".format(ensemble_name)
            pyrvapi.rvapi_add_text(header, container_id, 0, 0, 1, 1)

            sec_table = "sec_table_{0}".format(ensemble_name) + uid
            title = "Results table: {0}".format(ensemble_name)
            title = "Summary"
            pyrvapi.rvapi_add_section(sec_table, title, container_id, 0, 0, 1, 1, True)
            table_id = "table_{0}".format(ensemble_name) + uid
            pyrvapi.rvapi_add_table(table_id, "", sec_table, 1, 0, 1, 1, False)
            tdata = mrbump_util.ResultsSummary().results_table([r])
            self.fill_table(table_id, tdata, tooltips=self._mrbump_tooltips)
            # Ensemble
            if ensemble_results:
                epdb = self.ensemble_pdb(r, ensemble_results)
                if epdb:
                    sec_ensemble = "sec_ensemble_{0}".format(ensemble_name) + uid
                    pyrvapi.rvapi_add_section(sec_ensemble, "Ensemble Search Model", container_id, 0, 0, 1, 1, False)
                    data_ensemble = "data_ensemble_{0}".format(ensemble_name) + uid
                    pyrvapi.rvapi_add_data(
                        data_ensemble, "Ensemble PDB", self.fix_path(epdb), "XYZOUT", sec_ensemble, 2, 0, 1, 1, True
                    )
            # PHASER
            self.add_results_section(
                result_dict=r,
                ensemble_name=ensemble_name,
                program_name='PHASER',
                logfile_key='PHASER_logfile',
                pdb_key='PHASER_pdbout',
                mtz_key='PHASER_mtzout',
                uid=uid,
                container_id=container_id,
            )
            # REFMAC
            self.add_results_section(
                result_dict=r,
                ensemble_name=ensemble_name,
                program_name='Refmac',
                logfile_key='REFMAC_logfile',
                pdb_key='REFMAC_pdbout',
                mtz_key='REFMAC_mtzout',
                uid=uid,
                container_id=container_id,
            )
            # Buccaner
            self.add_results_section(
                result_dict=r,
                ensemble_name=ensemble_name,
                program_name='BUCCANEER',
                logfile_key='BUCC_logfile',
                pdb_key='BUCC_pdbout',
                mtz_key='BUCC_mtzout',
                uid=uid,
                container_id=container_id,
            )
            # Arpwarp
            self.add_results_section(
                result_dict=r,
                ensemble_name=ensemble_name,
                program_name='ArpWarp',
                logfile_key='ARP_logfile',
                pdb_key='ARP_pdbout',
                mtz_key='ARP_mtzout',
                uid=uid,
                container_id=container_id,
            )
            # SHELXE
            self.add_results_section(
                result_dict=r,
                ensemble_name=ensemble_name,
                program_name='SHELXE',
                logfile_key='SHELXE_logfile',
                pdb_key='SHELXE_pdbout',
                mtz_key='SHELXE_mtzout',
                uid=uid,
                container_id=container_id,
            )
            # Buccaner Rebuild
            self.add_results_section(
                result_dict=r,
                ensemble_name=ensemble_name,
                program_name='BUCCANEER SHELXE Trace Rebuild',
                logfile_key='SXRBUCC_logfile',
                pdb_key='SXRBUCC_pdbout',
                mtz_key='SXRBUCC_mtzout',
                uid=uid,
                container_id=container_id,
            )
            # Arpwarp Rebuild
            self.add_results_section(
                result_dict=r,
                ensemble_name=ensemble_name,
                program_name='ARPWARP SHELXE Trace Rebuild',
                logfile_key='SXRARP_logfile',
                pdb_key='SXRARP_pdbout',
                mtz_key='SXRARP_mtzout',
                uid=uid,
                container_id=container_id,
            )
            pyrvapi.rvapi_set_tree_node(results_tree, container_id, "{0}".format(ensemble_name), "auto", "")
        return

    def add_results_section(
        self,
        result_dict=None,
        ensemble_name=None,
        program_name=None,
        logfile_key=None,
        pdb_key=None,
        mtz_key=None,
        uid=None,
        container_id=None,
    ):
        assert (
            result_dict
            and ensemble_name
            and program_name
            and logfile_key
            and pdb_key
            and mtz_key
            and uid
            and container_id
        )
        have_logfile = have_files(result_dict, logfile_key)
        have_pdb_and_mtz = have_files(result_dict, pdb_key, mtz_key)
        if not (have_logfile or have_pdb_and_mtz):
            return
        program_id = program_name.lower().replace(' ', '_')
        this_sec_id = "sec_{0}_{1}".format(program_id, ensemble_name) + uid
        pyrvapi.rvapi_add_section(this_sec_id, "{} Outputs".format(program_name), container_id, 0, 0, 1, 1, False)
        if have_pdb_and_mtz:
            data_id = "o{0}{1}".format(program_id, ensemble_name) + uid
            pyrvapi.rvapi_add_data(
                data_id,
                "{} OUTPUTS".format(program_name),
                self.fix_path(result_dict[pdb_key]),
                "xyz",
                this_sec_id,
                2,
                0,
                1,
                1,
                True,
            )
            pyrvapi.rvapi_append_to_data(data_id, self.fix_path(result_dict[mtz_key]), "hkl:map")
        if have_logfile:
            data_id = "l{0}{1}".format(program_id, ensemble_name) + uid
            pyrvapi.rvapi_add_data(
                data_id,
                "{} Logfile".format(program_name),
                self.fix_path(result_dict[logfile_key]),
                # "summary",
                "text",
                this_sec_id,
                2,
                0,
                1,
                1,
                True,
            )

    def create_summary_tab(self, ample_dict):
        self._create_summary_tab()
        if self.do_create_ensembles_section(ample_dict):
            self.create_ensembles_section(ample_dict)
        if not self._got_mrbump_results(ample_dict):
            return self.summary_tab_id
        if not self.summary_tab_results_sec_id:
            self.rm_pending_section()
            # Only create the table once
            self.summary_tab_results_sec_id = "mrbump"
            pyrvapi.rvapi_add_section(self.summary_tab_results_sec_id, "MRBUMP", self.summary_tab_id, 0, 0, 1, 1, True)
            self.summary_tab_results_sec_table_id = "mrbump_table"
            pyrvapi.rvapi_add_table1(
                self.summary_tab_results_sec_id + "/" + self.summary_tab_results_sec_table_id,
                "MRBUMP Results",
                1,
                0,
                1,
                1,
                True,
            )
        mrb_results = ample_dict.get('mrbump_results')
        if not mrb_results == self.old_mrbump_results:
            # We set old_mrbump_results when we create the results_tab
            self.fill_table(
                self.summary_tab_results_sec_table_id,
                mrbump_util.ResultsSummary().results_table(mrb_results),
                tooltips=self._mrbump_tooltips,
            )
        if not self.summary_tab_survey_sec_id and not self.ccp4i2:
            # Only create the table once
            self.summary_tab_survey_sec_id = "survey"
            pyrvapi.rvapi_add_section(self.summary_tab_survey_sec_id, "Feedback", self.summary_tab_id, 0, 0, 1, 1, True)
            rstr = "<h2>How did we do?</h2><h3>Please follow this link and leave some feedback:</h3><a href='{0}' style='color: blue'>{0}</a>".format(
                reference_manager.survey_url
            )
            pyrvapi.rvapi_add_text(rstr, self.summary_tab_survey_sec_id, 0, 0, 1, 1)
        return self.summary_tab_id

    def _create_summary_tab(self):
        if not self.summary_tab_id:
            self.summary_tab_id = "summary_tab"
            title = "Summary"
            pyrvapi.rvapi_insert_tab(self.summary_tab_id, title, self.citation_tab_id, False)
            # Create pending section until we have data to show
            self.summary_tab_pending_sec_id = 'summary_tab_pending'
            pyrvapi.rvapi_add_section(
                self.summary_tab_pending_sec_id, "Processing...", self.summary_tab_id, 0, 0, 1, 1, True
            )
            rstr = "<p>No results are currently available. Please check back later.</p>"
            pyrvapi.rvapi_add_text(rstr, self.summary_tab_pending_sec_id, 0, 0, 1, 1)
        return

    def do_create_ensembles_section(self, ample_dict):
        return (
            not (ample_dict.get('single_model_mode') or ample_dict.get('homologs') or ample_dict.get('ideal_helices'))
            and bool(ample_dict.get('ensembles_data'))
            and not self.summary_tab_ensemble_sec_id
        )

    def create_ensembles_section(self, ample_dict):
        self.rm_pending_section()
        ensembles_data = ample_dict['ensembles_data']
        self.summary_tab_ensemble_sec_id = "ensembles"
        pyrvapi.rvapi_add_section(self.summary_tab_ensemble_sec_id, "Ensembles", self.summary_tab_id, 0, 0, 1, 1, True)
        if ample_dict['import_ensembles']:
            rstr = 'Imported {0} ensembles.'.format(len(ensembles_data))
            pyrvapi.rvapi_add_text(rstr, self.summary_tab_ensemble_sec_id, 0, 0, 1, 1)
        else:
            # Get the ensembling data
            d = ensembler.collate_cluster_data(ensembles_data)
            clusters = d['clusters']

            rstr = ""
            rstr += "Ensemble Results<br/>"
            rstr += "----------------<br/><br/>"
            rstr += "Cluster method: {0}<br/>".format(d['cluster_method'])
            rstr += "Cluster score type: {0}<br/>".format(d['cluster_score_type'])
            rstr += "Truncation method: {0}<br/>".format(d['truncation_method'])
            rstr += "Percent truncation: {0}<br/>".format(d['percent_truncation'])
            rstr += "Side-chain treatments: {0}<br/>".format(d['side_chain_treatments'])
            rstr += "Number of clusters: {0}<br/><br/>".format(len(clusters.keys()))
            rstr += "Generated {0} ensembles<br/><br/>".format(len(ensembles_data))
            pyrvapi.rvapi_add_text(rstr, self.summary_tab_ensemble_sec_id, 0, 0, 1, 1)
            ensemble_table = "ensemble_table"
            pyrvapi.rvapi_add_table1(
                self.summary_tab_ensemble_sec_id + "/" + ensemble_table, "Ensembling Results", 1, 0, 1, 1, True
            )
            tdata = []
            for i, cluster_num in enumerate(sorted(d['clusters'].keys())):
                header = True if i == 0 else False
                tdata += ensembler.cluster_table_data(clusters, cluster_num, d['side_chain_treatments'], header=header)
            self.fill_table(ensemble_table, tdata, tooltips=self._ensemble_tooltips)
        return

    def display_results(self, ample_dict):
        """Display the results of an AMPLE run using pyrvapi

        Parameters
        ----------
        ample_dict : dict
          An AMPLE job dictionary

        """
        if not (pyrvapi or self.generate_output):
            return
        try:
            if not self.header:
                pyrvapi.rvapi_add_header("AMPLE Results")
                self.header = True
            self.create_log_tab(ample_dict)
            self.create_citation_tab(ample_dict)
            self.create_summary_tab(ample_dict)
            self.create_results_tab(ample_dict)
            pyrvapi.rvapi_flush()
        except Exception as e:
            logger.critical("Error displaying results: %s\n%s", e, traceback.format_exc())
        return True

    def ensemble_pdb(self, mrbump_result, ensembles_data):
        try:
            ensemble_dict = None
            for e in ensembles_data:
                if e['name'] == mrbump_result.get('ensemble_name'):
                    ensemble_dict = e
                    break
            if os.path.isfile(ensemble_dict['ensemble_pdb']):
                return ensemble_dict['ensemble_pdb']
            else:
                return False
        except:
            return False

    def fix_path(self, path):
        """Ammend path so it's suitable for the webserver or jscofe/standalone"""
        if self.webserver_uri:
            return urlparse.urljoin(self.webserver_uri, path[self._webserver_start :])
        elif self.jscofe:
            return os.path.join("..", os.path.relpath(path, self.jsrview_dir))
        return path

    def fill_table(self, table_id, tdata, tooltips={}):
        # Make column headers
        for i in range(len(tdata[0])):  # Skip name as it's the row header
            h = tdata[0][i]
            tt = tooltips[h] if h in tooltips else ""
            pyrvapi.rvapi_put_horz_theader(table_id, h.encode('utf-8'), tt, i)  # Add table data
        for i in range(1, len(tdata)):
            for j in range(len(tdata[i])):
                pyrvapi.rvapi_put_table_string(table_id, str(tdata[i][j]), i - 1, j)
        # REM - can use pyrvapi.rvapi_shape_table_cell to format cells is required
        return

    def _got_mrbump_results(self, ample_dict):
        return ample_dict.get('mrbump_results') and len(ample_dict['mrbump_results'])

    def rm_pending_section(self):
        if self.summary_tab_pending_sec_id:
            pyrvapi.rvapi_flush()
            pyrvapi.rvapi_remove_widget(self.summary_tab_pending_sec_id)
            pyrvapi.rvapi_flush()
            self.summary_tab_pending_sec_id = None
        return

    def rvapi_shutdown(self, amopt):
        """Return any results to jscofe

        Parameters
        ----------
        amopt : dict
            AMPLE results dictionary with all information
        """
        rvdoc = amopt['rvapi_document']
        if not rvdoc:
            return

        # Create dictionary we're going to return
        meta = {'results': []}
        nresults = 0
        if bool(amopt.get('mrbump_results')):
            mrb_results = amopt['mrbump_results']
            nresults = min(3, len(mrb_results))
        if nresults > 0:
            for fdata in mrbump_util.ResultsSummary(mrb_results[:nresults]).topFiles(nresults):
                # Mangle paths. relpath assumes args are directories so need to add ..
                fdata['pdb'] = self.fix_path(fdata['pdb'])
                fdata['mtz'] = self.fix_path(fdata['mtz'])
                meta['results'].append(fdata)

        # Commit to file
        logger.debug("Exporting pyrvapi metadata:\n{0}".format(meta))
        pyrvapi.rvapi_put_meta(json.dumps(meta))
        pyrvapi.rvapi_store_document2(rvdoc)
        return


if __name__ == "__main__":
    import copy, sys, time
    from ample.util import ample_util

    logging.basicConfig(level=logging.DEBUG)
    pklfile = sys.argv[1]
    ample_dict = ample_util.read_amoptd(pklfile)
    ample_dict['show_gui'] = True
    ample_dict['ample_log'] = os.path.abspath(__file__)

    report_dir = os.path.abspath(os.path.join(os.curdir, "pyrvapi_tmp"))
    AR = AmpleOutput(ample_dict)
    AR.display_results(ample_dict)

    view1_dict = copy.copy(ample_dict)
    del view1_dict['ensembles_data']
    del view1_dict['mrbump_results']

    SLEEP = 5

    AR.display_results(view1_dict)
    time.sleep(SLEEP)

    # for i in range(10):
    view1_dict['ensembles_data'] = ample_dict['ensembles_data']
    AR.display_results(view1_dict)
    time.sleep(SLEEP)

    mrbump_results = []
    for r in ample_dict['mrbump_results'][0:3]:
        r['SHELXE_CC'] = None
        r['SHELXE_ACL'] = None
        mrbump_results.append(r)
    view1_dict['mrbump_results'] = mrbump_results
    AR.display_results(view1_dict)
    time.sleep(SLEEP)

    view1_dict['mrbump_results'] = ample_dict['mrbump_results'][0:5]
    AR.display_results(view1_dict)
    time.sleep(SLEEP)

    view1_dict['mrbump_results'] = ample_dict['mrbump_results']
    AR.display_results(view1_dict)
