"""
    AMPLE.py: CCP4 GUI Project
    
    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.
    
    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    """

import cPickle
import os
from CCP4PluginScript import CPluginScript
from lxml import etree

# AMPLE imports
from ample.util import mrbump_util
from ample.util.ample_util import I2DIR

AMPLE_ROOT_NODE = 'AMPLE'
AMPLE_LOG_NODE = 'LogText'
LOGFILE_NAME = 'log.txt'
LOGFILE_NAME = os.path.join('AMPLE_0','AMPLE.log')

class AMPLE(CPluginScript):
    TASKNAME = 'AMPLE'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'jens.thomas@liv.ac.uk'
    ERROR_CODES = { 1 : {'description' : 'Something not very good has happened.' },
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                       ['log_mtzjoin.txt', 0]
                       ]
    TASKCOMMAND="ample"
    
    def __init__(self, *args, **kws):
        super(AMPLE, self).__init__(*args, **kws)

    def processInputFiles(self):
        #Preprocess reflections to generate an "HKLIN" file
        '''
        #makeHklin0 takes as arguments a list of sublists
        #Each sublist comprises 1) A reflection data object identifier (one of those specified in the inputData container 
        #                           the task in the corresponding .def.xml
        #                       2) The requested data representation type to be placed into the file that is generated
        #
        #makeHklin0 returns a tuple comprising:
        #                       1) the file path of the file that has been created
        #                       2) a list of strings, each of which contains a comma-separated list of column labels output from
        #                       the input data objects
        #                       3) A CCP4 Error object       
        ''' 
        import CCP4XtalData
        self.hklin, self.columns, error = self.makeHklin0([
            ['AMPLE_F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]
        ])
        self.columnsAsArray = self.columns.split(",")

        self.fasta = self.container.inputData.AMPLE_SEQIN
        self.mtz = self.container.inputData.AMPLE_F_SIGF
        
        import CCP4ErrorHandling
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        
        #Preprocess coordinates to extract a subset
        '''
        # The method "getSelectedAtomsPdbFile" applied to a coordinate data object
        # selects those atoms declared in the objects "selectionString" property and writes them into
        # a pruned down file, the name of which is provided in the argument
        self.selectedCoordinatesPath = os.path.join(self.getWorkDirectory(), "selected_xyzin.pdb")
        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.selectedCoordinatesPath)
        '''
        
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
        #self.appendCommandLine(self.getWorkDirectory())
        self.appendCommandLine('-fasta')
        self.appendCommandLine( self.fasta)
        self.appendCommandLine('-mtz')
        self.appendCommandLine( self.hklin)
        #self.appendCommandLine( self.mtz)
        self.appendCommandLine('-F')
        self.appendCommandLine( self.columnsAsArray[0])
        self.appendCommandLine('-SIGF')
        self.appendCommandLine( self.columnsAsArray[1])
        #self.appendCommandLine('-ideal_helices')
        #self.appendCommandLine('True')
        self.appendCommandLine('-models')
        self.appendCommandLine('/opt/ample.git/testfiles/models')
        #self.appendCommandLine('-nproc')
        #self.appendCommandLine(str(self.container.controlParameters.AMPLE_NPROC))
        #self.appendCommandLine(str(self.container.inputData.AMPLE_NPROC))
        self.appendCommandLine('-do_mr')
        self.appendCommandLine(False)

        self.xmlroot = etree.Element(AMPLE_ROOT_NODE)
        logFile = os.path.join(self.getWorkDirectory(),LOGFILE_NAME)
        #self.watchFile(logFile,self.handleLogChanged)
                
        return CPluginScript.SUCCEEDED

    def handleLogChanged(self, filename):
        with open(os.path.join(self.getWorkDirectory(),'foo.txt'),'a') as w:
            w.write('flushXML: {0}\n'.format(self.makeFileName('PROGRAMXML')))
        for ampleTxtNode in self.xmlroot.xpath(AMPLE_LOG_NODE):
            self.xmlroot.remove(ampleTxtNode)
        element = etree.SubElement(self.xmlroot,AMPLE_LOG_NODE)
        with open (filename,'r') as logFile:
            element.text = etree.CDATA(logFile.read())
        self.flushXML()

    def flushXML(self):
        tmpFilename = self.makeFileName('PROGRAMXML')+'_tmp'
        with open(tmpFilename,'w') as xmlFile:
            xmlFile.write(etree.tostring(self.xmlroot,pretty_print=True))
        if os.path.exists(self.makeFileName('PROGRAMXML')):
            os.remove(self.makeFileName('PROGRAMXML'))
        os.rename(tmpFilename, self.makeFileName('PROGRAMXML'))

    def processOutputFiles(self):
        #Associate the tasks output coordinate file with the output coordinate object XYZOUT:
        '''
        self.container.outputData.XYZOUT.setFullPath(os.path.join(self.getWorkDirectory(),"final.pdb"))
        '''
        
        #debug_console()
        
        # Split an MTZ file into minimtz data objects
        '''
        outputFilesToMake = ['FPHIOUT','DIFFPHIOUT']
        columnsToTake = ['FWT,PHWT','DELFWT,PHDELWT']
        infile = os.path.join(self.workDirectory,'final.mtz')
        error = self.splitHklout(outputFilesToMake, columnsToTake, infile=infile)
        import CCP4ErrorHandling
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        '''
        
        #Create (dummy) PROGRAMXML, which basically contains only the log text of the job
        #without this, a report will not be generated
#         logfilePath = os.path.join(self.getWorkDirectory(),LOGFILE_NAME)
#         with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
#             xmlStructure = etree.Element(AMPLE_ROOT_NODE)
#             logText = etree.SubElement(xmlStructure,AMPLE_LOG_NODE)
#             with open(logfilePath,"r") as logFile:
#                 logText.text = etree.CDATA(logFile.read())
#             programXMLFile.write(etree.tostring(xmlStructure))
        
        # results_summary.sortResults(mrb_results, prioritise="SHELXE_CC")[0:min(len(mrb_results),mrbump_util.TOP_KEEP)],
        top_files = mrbump_util.ResultsSummary(results_pkl=os.path.join(self.getWorkDirectory(), I2DIR, 'resultsd.pkl')).topFiles()
        if top_files:
            for d in top_files:
                self.container.outputData.XYZOUT.append(d['xyz'])
                self.container.outputData.XYZOUT[-1].annotation = 'PDB file of ' + d['info']
                self.container.outputData.HKLOUT.append(d['hkl'])
                self.container.outputData.HKLOUT[-1].annotation = 'MTZ file of ' + d['info']

        logPath = os.path.join(self.getWorkDirectory(),LOGFILE_NAME)
        if os.path.isfile(logPath):
            with open(logPath, 'r') as logFile:
                element = etree.SubElement(self.xmlroot,AMPLE_LOG_NODE)
                element.text = etree.CDATA(logFile.read())
        self.flushXML()

        return CPluginScript.SUCCEEDED
