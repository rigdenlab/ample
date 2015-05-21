'''
Created on 3 Mar 2015

@author: jmht

'''

import cPickle
import logging
import os
import subprocess
import sys
import urlparse

sys.path.insert(0,"/opt/ample-dev1/python")
import ensemble
import mrbump_results

try: import pyrvapi
except: pyrvapi=None

_running=None
_tabs=None
_webserver_uri=None
_wbeserver_start=None

_tooltips={
           "ensemble_name" : "The identifier of the AMPLE ensemble search model",
           "MR_program" : "Molecular replacement program",
           "Solution_Type" : "MRBUMP categorisation of the solution",
           "PHASER_LLG" : "PHASER Log-likelihood gain for the Molecular Replacment solution",
           "PHASER_TFZ" : "PHASER Translation Function Z-score for the Molecular Replacment solution",
           "REFMAC_Rfact" : "Rfact score for REFMAC refinement of the Molecular Replacement solution",
           "REFMAC_Rfree" : "Rfree score for REFMAC refinement of the Molecular Replacement solution",
           "SHELXE_CC" : "SHELXE Correlation Coefficient score after C-alpha trace",
           "SHELXE_ACL" : "Average Chain Length of the fragments of the SHELXE C-alpha trace",
           "SXRBUCC_final_Rfact" : "Rfact score for BUCCANEER rebuild of the SHELXE C-alpha trace",
           "SXRBUCC_final_Rfree" : "Rfact score for BUCCANEER rebuild of the SHELXE C-alpha trace",
           "SXRARP_final_Rfact" : "Rfact score for ARPWARP rebuild of the SHELXE C-alpha trace",
           "SXRAP_final_Rfree" : "Rfact score for ARPWARP rebuild of the SHELXE C-alpha trace",
           }

def ensemble_pdb(mrbump_result,results_dict):
    ed=None
    for e in results_dict['ensembles_data']:
        if e['name']==mrbump_result['ensemble_name']:
            ed=e
            break
    assert ed,"Could not find ensemble!"
    return ed['ensemble_pdb']

def fix_path(path):
    """Ammend path so it's suitable for the webserver"""
    global _webserver_uri,_webserver_start
    if _webserver_uri:
        return urlparse.urljoin(_webserver_uri,path[_webserver_start:])
    else: return path
    
def fill_table(table_id, tdata):
    global _tooltips
    # Make column headers
    for i in range(len(tdata[0])): # Skip name as it's the row header
        h=tdata[0][i]
        tt=_tooltips[h] if h in _tooltips else ""
        pyrvapi.rvapi_put_horz_theader(table_id, h.encode('utf-8'), tt, i) # Add table data
    
    ir=len(tdata)-1 if len(tdata) > 2 else 2
    for i in range(1, ir):
        for j in range(len(tdata[i])):
            pyrvapi.rvapi_put_table_string(table_id, str(tdata[i][j]), i - 1, j)
    
    # Now colour the ensemble name cells
    for i in range(len(tdata) - 1):
        pyrvapi.rvapi_shape_table_cell(table_id, #tableId
            i, #row
            0, #column
            "", #tooltip
            "", # cell_css
            "table-blue-vh", #cell_style
            1, # rowSpan
            1) # colSpan
    return

def summary_tab(results_dict):
    #
    # Summary Tab
    #
    if not('ensembles_data' in results_dict and len(results_dict['ensembles_data'])): return
    ensembles_data=results_dict['ensembles_data']
    
    summary_tab="summary_tab"
    pyrvapi.rvapi_add_tab(summary_tab,"Summary",True) # Last arg is "open" - i.e. show or hide
    ensemble_sec="ensembles"
    pyrvapi.rvapi_add_section(ensemble_sec,"Ensembles",summary_tab,0,0,1,1,False)
    
    # Get the ensembling data
    clusters, cluster_method, truncation_method, percent_truncation = ensemble.collate_cluster_data(ensembles_data)
    
    rstr = ""
    rstr += "Ensemble Results<br/>"
    rstr += "----------------<br/><br/>"
    rstr += "Cluster method: {0}<br/>".format(cluster_method)
    rstr += "Truncation method: {0}<br/>".format(truncation_method)
    rstr += "Percent truncation: {0}<br/>".format(percent_truncation)
    rstr += "Number of clusters: {0}<br/><br/>".format(len(clusters.keys()))
    rstr += "Generated {0} ensembles<br/><br/>".format(len(ensembles_data))
    pyrvapi.rvapi_add_text(rstr,ensemble_sec,0,0,1,1 )
    
    ensemble_table="ensemble_table"
    pyrvapi.rvapi_add_table1(ensemble_sec+"/"+ensemble_table,"Ensembling Results",1,0,1,1,True)
    # for cluster_num in sorted(clusters.keys()):
    #     rstr += "\n"
    #     rstr += "Cluster {0}\n".format(cluster_num)
    #     rstr += "Number of models: {0}\n".format(clusters[cluster_num]['cluster_num_models'])
    #     rstr += "Cluster centroid: {0}\n".format(clusters[cluster_num]['cluster_centroid'])
    #     rstr += "\n"
    #     tdata = cluster_table_data(clusters, cluster_num)
    #     rstr += tableFormat.pprint_table(tdata)        
    # 
    cluster_num=1
    tdata = ensemble.cluster_table_data(clusters, cluster_num)
    fill_table(ensemble_table, tdata)
    
    #
    # MRBUMP Results
    #
    if not( 'mrbump_results' in results_dict and len(results_dict['mrbump_results'])): return summary_tab
    mrb_results=results_dict['mrbump_results']
    mrbump_sec="mrbump"
    pyrvapi.rvapi_add_section(mrbump_sec,"MRBUMP",summary_tab,0,0,1,1,True)
    mrbump_table="mrbump_table"
    pyrvapi.rvapi_add_table1(mrbump_sec+"/"+mrbump_table,"MRBUMP Results",1,0,1,1,True)
    mrb_data = mrbump_results.ResultsSummary().results_table(mrb_results)
    fill_table(mrbump_table, mrb_data)
    return summary_tab

def results_tab(results_dict):
    #
    # Results Tab
    #
    if not ('mrbump_results' in results_dict and len(results_dict['mrbump_results'])): return
    mrb_results=results_dict['mrbump_results']
    
    results_tab="results_tab"
    pyrvapi.rvapi_add_tab(results_tab,"Results",True) # Last arg is "open" - i.e. show or hide
    results_tree="results_tree"
    pyrvapi.rvapi_add_tree_widget(results_tree,"Final Results",results_tab,0,0,1,1)
    
    for r in mrb_results:
        name=r['ensemble_name']
        #container_id="sec_{0}".format(name)
        #pyrvapi.rvapi_add_section(container_id,"Results for: {0}".format(name),results_tree,0,0,1,1,True)
        container_id="sec_{0}".format(name)
        pyrvapi.rvapi_add_panel(container_id,results_tree,0,0,1,1)
        
        header="<h3>Results for ensemble: {0}</h3>".format(name)
        pyrvapi.rvapi_add_text(header,container_id,0,0,1,1 )
        
        sec_table="sec_table_{0}".format(name)
        title="Results table: {0}".format(name)
        title="Summary"
        pyrvapi.rvapi_add_section(sec_table,title,container_id,0,0,1,1,True)
        table_id="table_{0}".format(name)
        pyrvapi.rvapi_add_table(table_id,"",sec_table,1,0,1,1,False)
        tdata=mrbump_results.ResultsSummary().results_table([r])
        fill_table(table_id,tdata)
        
        # Ensemble
        epdb=ensemble_pdb(r,results_dict)
        sec_ensemble="sec_ensemble_{0}".format(name)
        pyrvapi.rvapi_add_section(sec_ensemble,"Ensemble Search Model",container_id,0,0,1,1,False )
        data_ensemble="data_ensemble_{0}".format(name)
        pyrvapi.rvapi_add_data(data_ensemble,
                                "Ensemble PDB",
                                fix_path(epdb),
                                "XYZOUT",
                                sec_ensemble,
                                2,0,1,1,True)
        # PHASER
        if os.path.isfile(str(r['PHASER_logfile'])) or (os.path.isfile(str(r['PHASER_pdbout'])) and os.path.isfile(str(r['PHASER_mtzout']))):
            sec_phaser="sec_phaser_{0}".format(name)
            pyrvapi.rvapi_add_section(sec_phaser,"PHASER Outputs",container_id,0,0,1,1,False )
            if os.path.isfile(str(r['PHASER_pdbout'])) and os.path.isfile(str(r['PHASER_mtzout'])):
                data_phaser="data_phaser_out_{0}".format(name)
                pyrvapi.rvapi_add_data(data_phaser,
                                        "PHASER PDB",
                                        fix_path(r['PHASER_pdbout']),
                                        "xyz:map",
                                        sec_phaser,
                                        2,0,1,1,True)
                pyrvapi.rvapi_append_to_data(data_phaser,fix_path(r['PHASER_mtzout']),"xyz:map")
            if os.path.isfile(str(r['PHASER_logfile'])):
                pyrvapi.rvapi_add_data("data_phaser_logfile_{0}".format(name),
                                        "PHASER Logfile",
                                        fix_path(r['PHASER_logfile']),
                                        "text",
                                        sec_phaser,
                                        2,0,1,1,True)
                
        # REFMAC
        if os.path.isfile(str(r['REFMAC_logfile'])) or (os.path.isfile(str(r['REFMAC_pdbout'])) and os.path.isfile(str(r['REFMAC_mtzout']))):
            sec_refmac="sec_refmac_{0}".format(name)
            pyrvapi.rvapi_add_section(sec_refmac,"REFMAC Outputs",container_id,0,0,1,1,False)
            if os.path.isfile(str(r['REFMAC_pdbout'])) and os.path.isfile(str(r['REFMAC_mtzout'])):
                data_refmac="data_refmac_out_{0}".format(name)
                pyrvapi.rvapi_add_data(data_refmac,
                                        "REFMAC PDB",
                                        fix_path(r['REFMAC_pdbout']),
                                        "xyz:map",
                                        sec_refmac,
                                        2,0,1,1,True)
                pyrvapi.rvapi_append_to_data(data_refmac,fix_path(r['REFMAC_mtzout']),"xyz:map")
            if os.path.isfile(str(r['REFMAC_logfile'])):
                pyrvapi.rvapi_add_data("data_refmac_logfile_{0}".format(name),
                                        "REFMAC Logfile",
                                        fix_path(r['REFMAC_logfile']),
                                        "text",
                                        sec_refmac,
                                        2,0,1,1,True)
                
        # SHELXE
        if os.path.isfile(str(r['SHELXE_logfile'])) or (os.path.isfile(str(r['SHELXE_pdbout'])) and os.path.isfile(str(r['SHELXE_mtzout']))):
            sec_shelxe="sec_shelxe_{0}".format(name)
            pyrvapi.rvapi_add_section(sec_shelxe,"SHELXE Outputs",container_id,0,0,1,1,False)
            if os.path.isfile(str(r['SHELXE_pdbout'])) and os.path.isfile(str(r['SHELXE_mtzout'])):
                data_shelxe="data_shelxe_out_{0}".format(name)
                pyrvapi.rvapi_add_data(data_shelxe,
                                        "SHELXE PDB",
                                        fix_path(r['SHELXE_pdbout']),
                                        "xyz:map",
                                        sec_shelxe,
                                        2,0,1,1,True)
                pyrvapi.rvapi_append_to_data(data_shelxe,fix_path(r['SHELXE_mtzout']),"xyz:map")
            if os.path.isfile(str(r['SHELXE_logfile'])):
                pyrvapi.rvapi_add_data("data_shelxe_logfile_{0}".format(name),
                                        "SHELXE Logfile",
                                        fix_path(r['SHELXE_logfile']),
                                        "text",
                                        sec_shelxe,
                                        2,0,1,1,True)
        
        # Buccaner Rebuild
        if os.path.isfile(str(r['SXRBUCC_logfile'])) or (os.path.isfile(str(r['SXRBUCC_pdbout'])) and os.path.isfile(str(r['SXRBUCC_mtzout']))):
            sec_sxrbucc="sec_sxrbucc_{0}".format(name)
            pyrvapi.rvapi_add_section(sec_sxrbucc,"BUCCANEER SHELXE Trace Rebuild Outputs",container_id,0,0,1,1,False)
            if os.path.isfile(str(r['SXRBUCC_pdbout'])) and os.path.isfile(str(r['SXRBUCC_mtzout'])):
                data_sxrbucc="data_sxrbucc_out_{0}".format(name)
                pyrvapi.rvapi_add_data(data_sxrbucc,
                                        "SXRBUCC PDB",
                                        fix_path(r['SXRBUCC_pdbout']),
                                        "xyz:map",
                                        sec_sxrbucc,
                                        2,0,1,1,True)
                pyrvapi.rvapi_append_to_data(data_sxrbucc,fix_path(r['SXRBUCC_mtzout']),"xyz:map")
            if os.path.isfile(str(r['SXRBUCC_logfile'])):
                pyrvapi.rvapi_add_data("data_shelxe_logfile_{0}".format(name),
                                        "SXRBUCC Logfile",
                                        fix_path(r['SXRBUCC_logfile']),
                                        "text",
                                        sec_sxrbucc,
                                        2,0,1,1,True)
                
        # Arpwarp Rebuild
        if os.path.isfile(str(r['SXRARP_logfile'])) or (os.path.isfile(str(r['SXRARP_pdbout'])) and os.path.isfile(str(r['SXRARP_mtzout']))):
            sec_sxrarp="sec_sxrarp_{0}".format(name)
            pyrvapi.rvapi_add_section(sec_sxrarp,"ARPWARP SHELXE Trace Redbuild Outputs",container_id,0,0,1,1,False)
            if os.path.isfile(str(r['SXRARP_pdbout'])) and os.path.isfile(str(r['SXRARP_mtzout'])):
                data_sxrarp="data_sxrarp_out_{0}".format(name)
                pyrvapi.rvapi_add_data(data_sxrarp,
                                        "SXRARP PDB",
                                        fix_path(r['SXRARP_pdbout']),
                                        "xyz:map",
                                        sec_sxrarp,
                                        2,0,1,1,True)
                pyrvapi.rvapi_append_to_data(data_sxrarp,fix_path(r['SXRARP_mtzout']),"xyz:map")
            if os.path.isfile(str(r['SXRARP_logfile'])):
                pyrvapi.rvapi_add_data("data_sxrarp_logfile_{0}".format(name),
                                        "SXRARP Logfile",
                                        fix_path(r['SXRARP_logfile']),
                                        "text",
                                        sec_sxrarp,
                                        2,0,1,1,True)
        
        pyrvapi.rvapi_set_tree_node("results_tree",container_id,"{0}".format(name),"auto","")
    return results_tab
        
def log_tab(results_dict):
    log_tab="log_tab"
    logfile=results_dict['ample_log']
    logurl=fix_path(logfile)
    pyrvapi.rvapi_add_tab(log_tab,"Log file",True) # Last arg is "open" - i.e. show or hide
    # Add watched (updatable) content to the log tab. Note that the
    # log file does not exist yet.
    pyrvapi.rvapi_append_content(logurl,True,log_tab)
    #pyrvapi.rvapi_flush()
    return log_tab

def display_results(results_dict,run_dir=None):

    global _running,_tabs,_webserver_uri,_webserver_start
    logger=logging.getLogger()
    if not pyrvapi:
        msg="Cannot display results using pyrvapi!"
        logger.debug(msg)
        return False
    
    # Remove all tabs
    if _tabs:
        for t in _tabs: pyrvapi.rvapi_remove_tab(t)

    if not _running:        
        # Infrastructure to run
        ccp4 = os.environ["CCP4"]
        share_jsrview = os.path.join(ccp4, "share", "jsrview")
        if not run_dir: run_dir=os.path.join(results_dict['work_dir'],"jsrview")
        if not os.path.isdir(run_dir): os.mkdir(run_dir)
        pyrvapi.rvapi_init_document ("AMPLE_results",run_dir,"AMPLE Results",1,7,share_jsrview,None,None,None)
        if 'webserver_uri' in results_dict and results_dict['webserver_uri']:
            # don't start browser and setup variables for the path on the webserver
            _webserver_start=len(results_dict['run_dir'])+1
            _webserver_uri=results_dict['webserver_uri']
        else:
            # We start our own browser
            jsrview = os.path.join(ccp4, "libexec", "jsrview")
            subprocess.Popen([jsrview, os.path.join(run_dir,"index.html")])
        pyrvapi.rvapi_add_header("AMPLE Results")
        #print "RUNNING display_results ",len(results_dict['mrbump_results']) if 'mrbump_results' in results_dict else "NO RESULTS"
        _running=True
        
    _tabs=[]
    t=summary_tab(results_dict)
    if t: _tabs.append(t)
    t=results_tab(results_dict)
    if t: _tabs.append(t)
    t=log_tab(results_dict)
    if t: _tabs.append(t)
    pyrvapi.rvapi_flush()
    return True

if __name__=="__main__":
    import time
    pklfile="/opt/ample-dev1/examples/toxd-example/resultsd1.pkl"
    with open(pklfile) as f: results_dict=cPickle.load(f)
    results_dict['webserver_uri']="http:www.jensrules.co.uk/ample/stuff"
    results_dict['webserver_uri']=None
    display_results(results_dict)
    print "TAB 1"
    time.sleep(10)
    pklfile="/opt/ample-dev1/examples/toxd-example/resultsd2.pkl"
    with open(pklfile) as f: results_dict=cPickle.load(f)
    results_dict['webserver_uri']="http:www.jensrules.co.uk/ample/stuff"
    results_dict['webserver_uri']=None
    display_results(results_dict)
    print "TAB 2"
    time.sleep(10)
    pklfile="/opt/ample-dev1/examples/toxd-example/resultsd3.pkl"
    with open(pklfile) as f: results_dict=cPickle.load(f)
    results_dict['webserver_uri']="http:www.jensrules.co.uk/ample/stuff"
    results_dict['webserver_uri']=None
    display_results(results_dict)
    print "TAB 3"

