'''
Created on 3 Mar 2015

@author: jmht

'''

import cPickle
import os
import subprocess
import sys

sys.path.insert(0,"/opt/ample-dev1/python")
import ensemble
import mrbump_results

import pyrvapi
# try: import pyrvapi
# except: pyrvapi=None

def fill_table(table_id, tdata):
    # Make column headers
    for i in range(len(tdata[0])): # Skip name as it's the row header
        pyrvapi.rvapi_put_horz_theader(table_id, tdata[0][i], "", i) # Add table data
    
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
    
    pyrvapi.rvapi_flush()
    return

def summary_tab(results_dict):
    #
    # Summary Tab
    #
    pyrvapi.rvapi_add_tab("summary_tab","Summary",True) # Last arg is "open" - i.e. show or hide
    #
    # Ensemble Results
    #
    pyrvapi.rvapi_add_section("ensembles","Ensembles","summary_tab",0,0,1,1,False)
    
    # Get the ensembling data
    ensembles_data=results_dict['ensembles_data']
    clusters, cluster_method, truncation_method, percent_truncation = ensemble.collate_cluster_data(ensembles_data)
    
    rstr = ""
    rstr += "Ensemble Results<br/>"
    rstr += "----------------<br/><br/>"
    rstr += "Cluster method: {0}<br/>".format(cluster_method)
    rstr += "Truncation method: {0}<br/>".format(truncation_method)
    rstr += "Percent truncation: {0}<br/>".format(percent_truncation)
    rstr += "Number of clusters: {0}<br/><br/>".format(len(clusters.keys()))
    rstr += "Generated {0} ensembles<br/><br/>".format(len(ensembles_data))
    pyrvapi.rvapi_add_text(rstr,"ensembles",0,0,1,1 )
    #
    
    pyrvapi.rvapi_add_table1("ensembles/ensemble_table","Ensembling Results",1,0,1,1,True)
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
    fill_table("ensemble_table", tdata)
    
    #
    # MRBUMP Results
    #
    pyrvapi.rvapi_add_section("mrbump","MRBUMP","summary_tab",0,0,1,1,True)
    pyrvapi.rvapi_add_table1("mrbump/mrbump_table","MRBUMP Results",1,0,1,1,True)
    mrb_data = mrbump_results.ResultsSummary().results_table(results_dict['mrbump_results'])
    fill_table("mrbump_table", mrb_data)
    pyrvapi.rvapi_flush()
    return

def results_tab(results_dict):
    #
    # Results Tab
    #
    pyrvapi.rvapi_add_tab("results_tab","Results",True) # Last arg is "open" - i.e. show or hide
    pyrvapi.rvapi_add_tree_widget("results_tree","Final Results","results_tab",0,0,1,1)
    
    mrb_results=results_dict['mrbump_results']
    #mrbump_results.ResultsSummary().results_table()
    for r in mrb_results:
        name=r['ensemble_name']
        sec_id="sec_{0}".format(name)
        pyrvapi.rvapi_add_section(sec_id,"Results for: {0}".format(name),"results_tree",0,0,1,1,True )
        sec_table="sec_table_{0}".format(name)
        pyrvapi.rvapi_add_section(sec_table,"Results table: {0}".format(name),sec_id,0,0,1,1,True )
        tableId="table_{0}".format(name)
        pyrvapi.rvapi_add_table(tableId,"",sec_table,1,0,1,1,False)
        tdata=mrbump_results.ResultsSummary().results_table([r])
        fill_table(tableId,tdata)
        
        # PHASER
        if r['PHASER_logfile'] or (r['PHASER_pdbout'] and r['PHASER_mtzout']):
            sec_phaser="sec_phaser_{0}".format(name)
            pyrvapi.rvapi_add_section(sec_phaser,"PHASER Outputs",sec_id,0,0,1,1,True )
            if r['PHASER_pdbout'] and r['PHASER_mtzout']:
                data_phaser="data_phaser_out_{0}".format(name)
                pyrvapi.rvapi_add_data(data_phaser,
                                        "PHASER PDB",
                                        r['PHASER_pdbout'],
                                        "xyz:map",
                                        sec_phaser,
                                        2,0,1,1,True)
                pyrvapi.rvapi_append_to_data(data_phaser,r['PHASER_mtzout'],"xyz:map" )
            if r['PHASER_logfile']:
                pyrvapi.rvapi_add_data("data_phaser_logfile_{0}".format(name),
                                        "PHASER Logfile",
                                        r['PHASER_logfile'],
                                        "text",
                                        sec_phaser,
                                        2,0,1,1,True)
                
        # REFMAC
        if r['REFMAC_logfile'] or (r['REFMAC_pdbout'] and r['REFMAC_mtzout']):
            sec_refmac="sec_refmac_{0}".format(name)
            pyrvapi.rvapi_add_section(sec_refmac,"REFMAC Outputs",sec_id,0,0,1,1,True )
            if r['REFMAC_pdbout'] and r['REFMAC_mtzout']:
                data_refmac="data_refmac_out_{0}".format(name)
                pyrvapi.rvapi_add_data(data_refmac,
                                        "REFMAC PDB",
                                        r['REFMAC_pdbout'],
                                        "xyz:map",
                                        sec_refmac,
                                        2,0,1,1,True)
                pyrvapi.rvapi_append_to_data(data_refmac,r['REFMAC_mtzout'],"xyz:map" )
            if r['REFMAC_logfile']:
                pyrvapi.rvapi_add_data("data_refmac_logfile_{0}".format(name),
                                        "REFMAC Logfile",
                                        r['REFMAC_logfile'],
                                        "text",
                                        sec_refmac,
                                        2,0,1,1,True)
                
        # SHELXE
        if r['SHELXE_logfile'] or (r['SHELXE_pdbout'] and r['SHELXE_mtzout']):
            sec_shelxe="sec_shelxe_{0}".format(name)
            pyrvapi.rvapi_add_section(sec_shelxe,"SHELXE Outputs",sec_id,0,0,1,1,True )
            if r['SHELXE_pdbout'] and r['SHELXE_mtzout']:
                data_shelxe="data_shelxe_out_{0}".format(name)
                pyrvapi.rvapi_add_data(data_shelxe,
                                        "SHELXE PDB",
                                        r['SHELXE_pdbout'],
                                        "xyz:map",
                                        sec_shelxe,
                                        2,0,1,1,True)
                pyrvapi.rvapi_append_to_data(data_shelxe,r['SHELXE_mtzout'],"xyz:map" )
            if r['SHELXE_logfile']:
                pyrvapi.rvapi_add_data("data_shelxe_logfile_{0}".format(name),
                                        "SHELXE Logfile",
                                        r['SHELXE_logfile'],
                                        "text",
                                        sec_shelxe,
                                        2,0,1,1,True)
                
        
        pyrvapi.rvapi_set_tree_node("results_tree",sec_id,"{0}".format(name),"auto","" )
        pyrvapi.rvapi_flush()
    return
        
def log_tab(logfile):
    pyrvapi.rvapi_add_tab("log_tab","Log file",True) # Last arg is "open" - i.e. show or hide
    
    # Add watched (updatable) content to the log tab. Note that the
    # log file does not exist yet.
    pyrvapi.rvapi_append_content(logfile,True,"log_tab" )
    pyrvapi.rvapi_flush()
    return

def display_results(results_dict,logfile,run_dir=None):
    # Infrastructure to run
    ccp4 = os.environ["CCP4"]
    share_jsrview = os.path.join(ccp4, "share", "jsrview")
    jsrview = os.path.join(ccp4, "libexec", "jsrview")
    if not run_dir: run_dir=results_dict['work_dir']
    pyrvapi.rvapi_init_document ("TestRun",run_dir,"RVAPI Demo 1",1,7,share_jsrview,None,None,None)
    subprocess.Popen([jsrview, os.path.join("index.html")])
    pyrvapi.rvapi_add_header("AMPLE Results")
    summary_tab(results_dict)
    results_tab(results_dict)
    log_tab(logfile)
    return


pklfile="/opt/ample-dev1/tests/testfiles/resultsd1.pkl"
with open(pklfile) as f: results_dict=cPickle.load(f)
logfile="/opt/ample-dev1/python/pyrvapi_results.py"
display_results(results_dict,logfile,run_dir=os.getcwd())

