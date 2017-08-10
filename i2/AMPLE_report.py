"""
    AMPLE_report.py: CCP4 GUI Project
    
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

from CCP4RvapiParser import RvapiReport

#from CCP4ReportParser import Report
# class AMPLE_report(Report):
#     # Specify which gui task and/or pluginscript this applies to
#     TASKNAME = 'AMPLE'
#     RUNNING = False
#     def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
#         Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, jobStatus=jobStatus, **kw)
#         clearingDiv = self.addDiv(style="clear:both;")
#         self.addDefaultReport(self)
#         
#     def addDefaultReport(self, parent=None):
#         if parent is None: parent=self
#         if len(self.xmlnode.xpath("LogText")) > 0:
#             newFold = parent.addFold(label="Log text", initiallyOpen=True)
#             newFold.addPre(text = self.xmlnode.xpath("LogText")[0].text)

class AMPLE_report(RvapiReport):
#class crank2_report(Report):
  TASKNAME="AMPLE"
  RUNNING = True
  SEPARATEDATA = True
  WATCHED_FILE = 'i2.xml'
  def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
      print "JMHT GOT RUNDIR ", jobInfo['fileroot']
      RvapiReport.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus,**kw)

