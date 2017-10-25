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
    * RESTRAINTS FILE  (OPT)
    =============================================
    NO MODELS - subFrane
    * ROSETTA_DIR
    * FRAGS 3
    * FRAGS 9
    * RESTRAINTS FILE (OPT)
    =============================================
    
    
    ADVANCED OPTIONS
    
    Need code to determine what type of models we've been given.
    Variables:
    AMPLE_PROTEIN_TYPE: Globular/Transmembrane
    AMPLE_EXISTING_MODELS: T/F
    AMPLE_MODEL_TYPE: abinitio, multiple_homologs, single_homolog, nmr_ensemble
    #GOT AMPLE_RESTRAINTS_FILE: PATH
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

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
#         self.createLine(['subtitle','Use SHELXE'])
#         x = self.container.inputData.AMPLE_USE_SHELXE.qualifiers()['guiLabel']
#         self.createLine( ['label', x, 'widget', 'AMPLE_USE_SHELXE'])

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

        # Existing Models
        self.openSubFrame(frame=False)
        self.createLine(['subtitle','Do you have existing models?', 'widget', 'AMPLE_EXISTING_MODELS'])
        self.closeSubFrame()     

        # Model Selection
        self.openSubFrame(frame=True, toggle=['AMPLE_EXISTING_MODELS', 'open', ['True'] ])
        self.createLine( ['subtitle', 'Models are from:', 'widget', '-guiMode', 'radio', 'AMPLE_MODELS_SOURCE' ])
        self.createLine( ['widget', 'AMPLE_MODELS_DIR'], toggle=['AMPLE_MODELS_SOURCE', 'open', ['directory']])
        self.createLine( ['widget', 'AMPLE_MODELS_FILE'], toggle=['AMPLE_MODELS_SOURCE', 'open', ['file']])
        self.createLine( ['subtitle' ,'What sort of models are these?', 'widget', 'AMPLE_MODEL_TYPE' ] )
        self.createLine(['subtitle','NMR remodelling:', 'widget', '-guiMode', 'radio', 'AMPLE_NMR_REMODEL'],
                        toggle=['AMPLE_MODEL_TYPE', 'open', ['nmr_ensemble'] ])
        self.closeSubFrame()

        self.openSubFrame(frame=True, toggle=['AMPLE_EXISTING_MODELS', 'close', ['True'] ])
        self.createLine(['subtitle','Model source:', 'widget', '-guiMode', 'radio', 'AMPLE_MODEL_GENERATION'])
        self.closeSubFrame() # Model Generation

        # Rosetta Paths
        self.openSubFrame(frame=True,
                          toggleFunction=[self.toggleRosettaFiles, ['AMPLE_EXISTING_MODELS','AMPLE_MODEL_GENERATION','AMPLE_NMR_REMODEL', 'AMPLE_MODEL_TYPE']])
        self.createLine(['subtitle', 'ROSETTA paths'])
        self.createLine(['advice', 'ROSETTA can be downloaded from: <a href="https://www.rosettacommons.org">https://www.rosettacommons.org</a>'])
        self.createLine( ['widget', 'AMPLE_ROSETTA_DIR'])
        self.createLine(['advice', 'Fragments should be created on the ROBETTA SERVER: <a href="http://robetta.bakerlab.org">http://robetta.bakerlab.org</a>'])
        self.createLine( ['widget', 'AMPLE_ROSETTA_FRAGS3'])
        self.createLine( ['widget', 'AMPLE_ROSETTA_FRAGS9'])
        self.createLine( ['widget', 'AMPLE_RESTRAINTS_FILE'])
        self.closeSubFrame()
        
        # Depending on the AMPLE_RUN_MODE we need to set certain files as required or not
        #self.connect(self.container.inputData.AMPLE_RUN_MODE,QtCore.SIGNAL('dataChanged'),self.AMPLE_RUN_MODEchanged)
        self.connect(self.container.inputData.AMPLE_EXISTING_MODELS,QtCore.SIGNAL('dataChanged'),self.setModelsRequired)
        self.connect(self.container.inputData.AMPLE_MODELS_SOURCE,QtCore.SIGNAL('dataChanged'),self.setModelsRequired)
        self.setModelsRequired()
        
        # Number of processors
        self.openSubFrame(frame=True)
        self.container.inputData.AMPLE_NPROC = cpu_count()
        x = self.container.inputData.AMPLE_NPROC.qualifiers()['guiLabel']
        self.createLine(['subtitle', x, 'widget','AMPLE_NPROC'])
        self.closeSubFrame()
        
        self.drawOptions()
    
    def drawOptions(self):
        folder = self.openFolder(folderFunction='inputData',title='Advanced Options')
        self.createLine(['subtitle', 'Rebuilding options' ])
        self.openSubFrame(frame=True)
        #x = self.container.inputData.AMPLE_ENSEMBLING_TYPE.qualifiers()['guiLabel']
        #self.createLine(['label', x, 'widget', 'AMPLE_ENSEMBLING_TYPE'])
        x = self.container.inputData.AMPLE_REFINE_REBUILD.qualifiers()['guiLabel']
        self.createLine(['label', x, 'widget', 'AMPLE_REFINE_REBUILD'])
        x = self.container.inputData.AMPLE_USE_SHELXE.qualifiers()['guiLabel']
        self.createLine(['label', x, 'widget', 'AMPLE_USE_SHELXE'])
        x = self.container.inputData.AMPLE_SHELXE_REBUILD.qualifiers()['guiLabel']
        self.createLine(['label', x, 'widget', 'AMPLE_SHELXE_REBUILD'])
        self.closeSubFrame()
        
        x = self.container.inputData.AMPLE_EXTRA_FLAGS.qualifiers()['guiLabel']
        self.createLine(['subtitle', x, 'widget', '-guiMode','multiLine', 'AMPLE_EXTRA_FLAGS' ])
        
            
    def toggleRosettaFiles(self):
        if self.container.inputData.AMPLE_EXISTING_MODELS == 'False'  and self.container.inputData.AMPLE_MODEL_GENERATION == 'rosetta' or \
            self.container.inputData.AMPLE_EXISTING_MODELS == 'True'  and self.container.inputData.AMPLE_MODEL_TYPE == 'nmr_ensemble' and self.container.inputData.AMPLE_NMR_REMODEL == 'nmr_remodel_true':
            self.container.inputData.AMPLE_ROSETTA_DIR.setQualifiers({'allowUndefined':False,'mustExist':True})
            self.container.inputData.AMPLE_ROSETTA_FRAGS3.setQualifiers({'allowUndefined':False,'mustExist':True})
            self.container.inputData.AMPLE_ROSETTA_FRAGS9.setQualifiers({'allowUndefined':False,'mustExist':True})
            self.getWidget('AMPLE_ROSETTA_DIR').validate()
            self.getWidget('AMPLE_ROSETTA_FRAGS3').validate()
            self.getWidget('AMPLE_ROSETTA_FRAGS9').validate()
            return True
        else:
            self.container.inputData.AMPLE_ROSETTA_DIR.setQualifiers({'allowUndefined':True,'mustExist':False})
            self.container.inputData.AMPLE_ROSETTA_FRAGS3.setQualifiers({'allowUndefined':True,'mustExist':False})
            self.container.inputData.AMPLE_ROSETTA_FRAGS9.setQualifiers({'allowUndefined':True,'mustExist':False})
            self.getWidget('AMPLE_ROSETTA_DIR').validate()
            self.getWidget('AMPLE_ROSETTA_FRAGS3').validate()
            self.getWidget('AMPLE_ROSETTA_FRAGS9').validate()
            return False
        
    def setModelsRequired(self):
        if self.container.inputData.AMPLE_EXISTING_MODELS == 'True':
            if self.container.inputData.AMPLE_MODELS_SOURCE == 'file':
                self.container.inputData.AMPLE_MODELS_FILE.setQualifiers({'allowUndefined':False,'mustExist':True})
                self.container.inputData.AMPLE_MODELS_DIR.setQualifiers({'allowUndefined':True,'mustExist':False})
            elif self.container.inputData.AMPLE_MODELS_SOURCE == 'directory':
                self.container.inputData.AMPLE_MODELS_DIR.setQualifiers({'allowUndefined':False,'mustExist':True})
                self.container.inputData.AMPLE_MODELS_FILE.setQualifiers({'allowUndefined':True,'mustExist':False})
            else: assert False,"Unrecognised Parameter: {0}".format(self.container.inputData.AMPLE_MODELS_SOURCE)
        else:
            self.container.inputData.AMPLE_MODELS_FILE.setQualifiers({'allowUndefined':True,'mustExist':False})
            self.container.inputData.AMPLE_MODELS_DIR.setQualifiers({'allowUndefined':True,'mustExist':False})
        self.getWidget('AMPLE_MODELS_DIR').validate()
        self.getWidget('AMPLE_MODELS_FILE').validate()
