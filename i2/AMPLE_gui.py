"""
    AMPLE_gui.py: CCP4 GUI Project
    
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

"""
Notes 

icons live in devel/qticons
for help.html file copy into devel/docs/tasks/AMPLE

"""

from CCP4TaskWidget import CTaskWidget
# David's handy debug console
from phil2etree import debug_console

from PyQt4 import QtCore
from multiprocessing import cpu_count

class AMPLE_gui(CTaskWidget):
    """
    
    FASTA
    MTZ
    PROTEIN TYPE
    HAVE MODELS? YES/NO
    
    ============================================
    YES MODELS - subFrame
    * MODEL TYPE
    ** AB INITIO
    ** HOMOLOGS
    ** SINGLE MODEL
    ** NMR ENSEMBLE
    
    IF SINGLE MODEL:
    * SCORE FILE
    ELIF NMR ENSEMBLE:
    * CONTACT FILE  (OPT)
    =============================================
    NO MODELS - subFrane
    * ROSETTA_DIR
    * FRAGS 3
    * FRAGS 9
    * CONTACT FILE (OPT)
    =============================================
    
    
    ADVANCED OPTIONS
    
    Need code to determine what type of models we've been given.
    Variables:
    AMPLE_PROTEIN_TYPE: Globular/Transmembrane
    AMPLE_EXISTING_MODELS: T/F
    AMPLE_MODEL_TYPE: abinitio, multiple_homologs, single_homolog, nmr_ensemble
    #GOT AMPLE_CONTACT_FILE: PATH
    AMPLE_SCORE_FILE: PATH
    
    """

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'AMPLE' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    TASKMODULE = [ 'molecular_replacement' ] #Section in the task list where this task will be listed e.g. 'refinement','model_building' for full list see MODULE_ORDER in core/CCP4TaskManager.py
    SHORTTASKTITLE='AMPLE Molecular Replacement Pipeline'
    TASKTITLE='Molecular Replacement with unconventional models - AMPLE'
    DESCRIPTION = '''This task is for running Molecular Replacement with unconventional models'''
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['coot_rebuild']
    
    def toggleIsFile(self):
        return self.container.inputData.AMPLE_RUN_MODE == 'nmr_ensemble' or \
            self.container.inputData.AMPLE_MODELS_SOURCE == 'file'
    def toggleIsDir(self):
        return self.container.inputData.AMPLE_MODELS_SOURCE == 'directory' and not \
            self.container.inputData.AMPLE_RUN_MODE == 'nmr_ensemble' 

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.openFolder(folderFunction='inputData',followFrom=False)

        self.createLine(['subtitle','Input sequence'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'tip','Input sequence','widget','AMPLE_SEQIN' ] )
        self.closeSubFrame()
        
        self.createLine(['subtitle','Input reflections'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'tip','Input reflections','widget','AMPLE_F_SIGF' ] )
        self.closeSubFrame()


        self.createLine(['subtitle','Class of protein:', 'widget', '-guiMode', 'radio', 'AMPLE_PROTEIN_CLASS'])


#         # AMPLE Run Mode
#         self.createLine(['subtitle','AMPLE Run Mode'])
#         self.openSubFrame(frame=True)
#         self.createLine( ['widget' ,'AMPLE_RUN_MODE' ] )
#         self.closeSubFrame()

#         # Model Selection
#         self.openSubFrame(frame=[True], toggle=['AMPLE_RUN_MODE', 'close', ['rosetta','ideal_helices'] ])
#         self.createLine(['subtitle','Model selection'])
#         self.createLine( ['widget', '-guiMode', 'radio', 'AMPLE_MODELS_SOURCE' ],
#                          toggle=['AMPLE_RUN_MODE', 'close', 'nmr_ensemble'] )
#         self.createLine( ['widget', 'AMPLE_MODELS_DIR'],
#                          toggleFunction=[self.toggleIsDir, ['AMPLE_RUN_MODE','AMPLE_MODELS_SOURCE']])
#         self.createLine( ['widget', 'AMPLE_MODELS_FILE'],
#                          toggleFunction=[self.toggleIsFile, ['AMPLE_RUN_MODE','AMPLE_MODELS_SOURCE']])
#         self.closeSubFrame()
#         
#         self.openSubFrame(frame=[True], toggle=['AMPLE_RUN_MODE', 'open', ['rosetta'] ])
#         self.createLine(['subtitle','Rosetta paths'])
#         self.createLine( ['widget', 'AMPLE_ROSETTA_DIR'])
#         self.createLine( ['widget', 'AMPLE_ROSETTA_FRAGS3'])
#         self.createLine( ['widget', 'AMPLE_ROSETTA_FRAGS9'])
#         self.createLine( ['widget', 'AMPLE_CONTACT_FILE'])
#         # Depending on the AMPLE_RUN_MODE we need to set certain files as required or not
#         self.connect(self.container.inputData.AMPLE_RUN_MODE,QtCore.SIGNAL('dataChanged'),self.AMPLE_RUN_MODEchanged)
#         self.AMPLE_RUN_MODEchanged()
#         self.closeSubFrame()
        
#         self.createLine(['subtitle','Use SHELXE'])
#         x = self.container.inputData.AMPLE_USE_SHELXE.qualifiers()['guiLabel']
#         self.createLine( ['label', x, 'widget', 'AMPLE_USE_SHELXE'])
        
        self.openSubFrame(frame=False)
        #self.createLine(['subtitle','Source of models:'])
        self.createLine(['subtitle','Model generation:', 'widget', 'AMPLE_EXISTING_MODELS'])
        #self.createLine(['widget', '-guiMode', 'radio', 'AMPLE_EXISTING_MODELS'])
        #self.createLine(['widget', '-guiLabel', 'LABEL', '-title', 'TITLE', 'AMPLE_EXISTING_MODELS'])
        #self.createLine(['widget', 'AMPLE_EXISTING_MODELS'])
        self.closeSubFrame() # Existing Models    

        # Model Selection
        #self.openSubFrame(frame=True, title='Models are from:', toggle=['AMPLE_EXISTING_MODELS', 'open', ['True'] ])
        self.openSubFrame(frame=True, toggle=['AMPLE_EXISTING_MODELS', 'open', ['True'] ])
        #self.createLine(['subtitle','Models are from:'])
        self.createLine( ['subtitle', 'Models are from:', 'widget', '-guiMode', 'radio', 'AMPLE_MODELS_SOURCE' ],
                         toggle=['AMPLE_RUN_MODE', 'close', 'nmr_ensemble'] )
        self.createLine( ['widget', 'AMPLE_MODELS_DIR'],
                         toggleFunction=[self.toggleIsDir, ['AMPLE_RUN_MODE','AMPLE_MODELS_SOURCE']])
        self.createLine( ['widget', 'AMPLE_MODELS_FILE'],
                         toggleFunction=[self.toggleIsFile, ['AMPLE_RUN_MODE','AMPLE_MODELS_SOURCE']])
        self.createLine( ['label' ,'Type of models:', 'widget', 'AMPLE_MODEL_TYPE' ] )
        self.closeSubFrame()

        self.openSubFrame(frame=True, title='Rosetta paths', toggle=['AMPLE_EXISTING_MODELS', 'close', ['True'] ])
        self.createLine( ['widget', 'AMPLE_ROSETTA_DIR'])
        self.createLine( ['widget', 'AMPLE_ROSETTA_FRAGS3'])
        self.createLine( ['widget', 'AMPLE_ROSETTA_FRAGS9'])
        self.createLine( ['widget', 'AMPLE_CONTACT_FILE'])
        # Depending on the AMPLE_RUN_MODE we need to set certain files as required or not
        #self.connect(self.container.inputData.AMPLE_RUN_MODE,QtCore.SIGNAL('dataChanged'),self.AMPLE_RUN_MODEchanged)
        #self.AMPLE_RUN_MODEchanged()
        self.closeSubFrame() # Model Selection

        self.drawOptions()
    
    def drawOptions(self):
        folder = self.openFolder(folderFunction='inputData',title='Options')
        self.createLine(['subtitle', 'Basic options' ])
        
        self.openSubFrame(frame=True)
        self.container.inputData.AMPLE_NPROC = cpu_count()
        x = self.container.inputData.AMPLE_NPROC.qualifiers()['guiLabel']
        self.createLine(['label', x, 'widget','AMPLE_NPROC'])
        self.closeSubFrame()

        self.openSubFrame(frame=True)
        x = self.container.inputData.AMPLE_ENSEMBLING_TYPE.qualifiers()['guiLabel']
        self.createLine(['label', x, 'widget', 'AMPLE_ENSEMBLING_TYPE'])
        self.closeSubFrame()

    def AMPLE_RUN_MODEchanged(self):
        if self.container.inputData.AMPLE_RUN_MODE == 'rosetta':
            self.container.inputData.AMPLE_ROSETTA_DIR.setQualifiers({'allowUndefined':False,'mustExist':True})
            self.container.inputData.AMPLE_ROSETTA_FRAGS3.setQualifiers({'allowUndefined':False,'mustExist':True})
            self.container.inputData.AMPLE_ROSETTA_FRAGS9.setQualifiers({'allowUndefined':False,'mustExist':True})
        else:
            self.container.inputData.AMPLE_ROSETTA_DIR.setQualifiers({'allowUndefined':True,'mustExist':False})
            self.container.inputData.AMPLE_ROSETTA_FRAGS3.setQualifiers({'allowUndefined':True,'mustExist':False})
            self.container.inputData.AMPLE_ROSETTA_FRAGS9.setQualifiers({'allowUndefined':True,'mustExist':False})
        self.getWidget('AMPLE_ROSETTA_DIR').validate()
        self.getWidget('AMPLE_ROSETTA_FRAGS3').validate()
        self.getWidget('AMPLE_ROSETTA_FRAGS9').validate()

