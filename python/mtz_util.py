'''
Created on 2 Dec 2014

@author: jmht
'''

import logging
import os
import shutil
import sys
import unittest

from iotbx import reflection_file_reader

import ample_util
import cif_parser

_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)

def delColumn(file_name, column, overwrite=True):
    """Delete a column from an mtz file and return a path to the file"""
    mtzDel = ample_util.filename_append(file_name, "d{0}".format(column) )
    cmd = [ "mtzutils", "hklin1", file_name, "hklout", mtzDel ]
    stdin = "EXCLUDE 1 {0}".format( column )
    logfile = os.path.join( os.getcwd(), "mtzutils.log" )
    retcode = ample_util.run_command(cmd, stdin=stdin, logfile=logfile)
    if retcode != 0:
        msg = "Error running mtzutils. Check the logfile: {0}".format(logfile)
        _logger.critical(msg)
        raise RuntimeError, msg
    
    if overwrite:
        shutil.move(mtzDel,file_name)
        return file_name
    else:
        return mtzDel

def addRfree(file_name,directory=None,overwrite=True):
    """Run uniqueify on mtz file to generate RFREE data column"""
    mtzUnique = ample_util.filename_append(file_name, "uniqueify", directory=directory)

    cmd = ['uniqueify', file_name, mtzUnique]
    #uniqueify {-p fraction} mydata.mtz.
    logfile = os.path.join( os.getcwd(), "uniqueify.log" )
    retcode = ample_util.run_command(cmd, logfile=logfile)
    if retcode != 0:
        msg = "Error running sfcif2mtz. Check the logfile: {0}".format(logfile)
        _logger.critical(msg)
        raise RuntimeError, msg

    if overwrite:
        shutil.move(mtzUnique,file_name)
        return file_name
    else:
        return mtzUnique

def getRfree(file_name):
    """Return the Rfree label"""

    reflection_file = reflection_file_reader.any_reflection_file(file_name=file_name)
    if not reflection_file.file_type()=="ccp4_mtz":
        msg="File is not of type ccp4_mtz: {0}".format(file_name)
        logging.critical(msg)
        raise RuntimeError,msg
    
    # Read the file
    content=reflection_file.file_content()
    return _getRfree(content)
    
def _getRfree(content):
    rfree_label=None
    #print "GOT ",content.column_labels()
    for label in content.column_labels():
        if 'free' in label.lower():
            column = content.get_column(label=label)
            #print "size ",column.array_size()
            #print "valid ",column.n_valid_values()
            selection_valid = column.selection_valid()
            flags = column.extract_values()
            sel_0 = (flags == 0)
            n0=( sel_0 & selection_valid).count(True)
            n1=(~sel_0 & selection_valid).count(True)
            #print "Number of 0 (work):",n0
            #print "Number of 1 (test):",n1
            #print float(n0)/float(n1)*100
            if n0>0 and n1>0:
                if rfree_label:
                    _logger.warning("FOUND >1 RFREE label in file!")
                rfree_label=label
                
    return rfree_label

def processReflectionFile( amoptd ):
    """Make sure we have a valid mtz file. If necessary convert a given cif file.
       Set the mtz variable in the given amoptd to the reflection file to use
       Return True if it all worked or raise an exception if it failed
    """

    # We've been given a sf_cif so convert to mtz
    if amoptd['sf_cif']:
        if not os.path.isfile( amoptd['sf_cif'] ):
            _logger.critical("Cannot find sf_cif file: {0}".format( amoptd['sf_cif'] ) )
            sys.exit(1)

        cp = cif_parser.CifParser()
        amoptd['mtz'] = cp.sfcif2mtz( amoptd['sf_cif'] )

    # Now have an mtz so check it's valid
    if not amoptd['mtz'] or not os.path.isfile( amoptd['mtz'] ):
        _logger.critical("Cannot find MTZ file: {0}".format( amoptd['mtz'] ) )
        sys.exit(1)

    # Get column label info
    reflection_file = reflection_file_reader.any_reflection_file(file_name=amoptd['mtz'])
    if not reflection_file.file_type()=="ccp4_mtz":
        _logger.critical("File is not of type ccp4_mtz: {0}".format( amoptd['mtz'] ) )
        sys.exit(1)
    
    # Read the file
    content=reflection_file.file_content()
    
    # Check any user-given flags
    for flag in ['F','SIGF','FREE']:
        if amoptd[flag] and amoptd[flag] not in content.column_labels():
            _logger.critical("Cannot find flag {0} label {1} in mtz file {2}".format( flag, amoptd[flag], amoptd['mtz'] ) )
            sys.exit(1)    
    
    # If any of the flags aren't given we set defaults based on what's in the file
    if not amoptd['F']:
        if 'F' not in content.column_types():
            _logger.critical("Cannot find column type F for flag F in mtz file: {0}".format( amoptd['mtz'] ) )
            sys.exit(1)
        amoptd['F']  = content.column_labels()[content.column_types().index('F')]
    if not amoptd['SIGF']:
        l='SIG'+amoptd['F']
        if not l in content.column_labels():
            _logger.critical("Cannot find column type {0} for flag SIGF in mtz file: {0}".format( l, amoptd['mtz'] ) )
            sys.exit(1)
        amoptd['SIGF']  = l
        
    if amoptd['FREE']:
        # Check is valid
        rfree=_getRfree(content)
        if not rfree or not rfree==amoptd['FREE']:
            _logger.critical("Given RFREE label {0} is not valid for mtz file: {0}".format( amoptd['FREE'], amoptd['mtz'] ) )
            sys.exit(1)
    else:
        # See if we can find a valid label in the file
        rfree=_getRfree(content)
        if not rfree:
            # Need to generate RFREE
            _logger.warning("Cannot find a valid FREE flag - running uniquefy to generate column with RFREE data." )
            amoptd['mtz'] = addRfree( amoptd['mtz'], directory=amoptd['work_dir'],overwrite=False)

            # Check file and get new FREE flag
            rfree=getRfree(amoptd['mtz'])
            if not rfree:
                _logger.critical("Cannot find valid rfree flag in mtz file {0} after running uniquiefy".format(amoptd['mtz']))
                sys.exit(1)
        amoptd['FREE']  = rfree

    return True


class Test(unittest.TestCase):
    """
    Unit test
    """

    def setUp(self):
        """
        Get paths need to think of a sensible way to do this
        """

        thisd =  os.path.abspath( os.path.dirname( __file__ ) )
        paths = thisd.split( os.sep )
        self.ampleDir = os.sep.join( paths[ : -1 ] )
        self.testfilesDir = os.sep.join( paths[ : -1 ] + [ 'tests', 'testfiles' ] )

        return

    def testProcessReflectionFile(self):
        """Get MTZ flags"""


        mtz = os.path.join( self.ampleDir, "examples", "toxd-example" , "1dtx.mtz" )


        d = { 'mtz'    : mtz,
              'sf_cif' : None,
              'F'      : None,
              'SIGF'   : None,
              'FREE'   : None,
               'work_dir': os.getcwd(),
             }

        processReflectionFile( d )

        self.assertEqual( 'FP', d['F'], "Correct F")
        self.assertEqual( 'SIGFP', d['SIGF'], "Correct SIGF")
        self.assertEqual( 'FreeR_flag', d['FREE'], "Correct FREE")

        return

    def testProcessReflectionFileNORFREE(self):
        """Get MTZ flags"""

        mtz = os.path.join( self.testfilesDir, "2uui_sigmaa.mtz" )

        d = { 'mtz'    : mtz,
              'sf_cif' : None,
              'F'      : None,
              'SIGF'   : None,
              'FREE'   : None,
              'work_dir': os.getcwd(),
             }

        processReflectionFile( d )

        self.assertEqual( 'F', d['F'], "Correct F")
        self.assertEqual( 'SIGF', d['SIGF'], "Correct SIGF")
        print "GOT ",d['FREE']
        self.assertEqual( 'FreeR_flag', d['FREE'], "Correct FREE")
        
        os.unlink('uniqueify.log')
        os.unlink('2uui_sigmaa_uniqueify.mtz')
        os.unlink('2uui_sigmaa_uniqueify.log')

        return

    def testProcessReflectionFileCIF(self):
        """Get MTZ flags"""

        cif = os.path.join( self.testfilesDir, "1x79-sf.cif" )

        d = { 'mtz'     : None,
              'sf_cif'  : cif,
              'F'       : None,
              'SIGF'    : None,
              'FREE'    : None,
              'work_dir': os.getcwd(),
             }

        processReflectionFile( d )

        self.assertEqual( 'FP', d['F'], "Correct F")
        self.assertEqual( 'SIGFP', d['SIGF'], "Correct SIGF")
        self.assertEqual( 'FreeR_flag', d['FREE'], "Correct FREE")

        os.unlink('cif2mtz.log')
        os.unlink('1x79-sf.mtz')
        os.unlink('mtzutils.log')
        #os.unlink('1x79-sf_dFREE.mtz')
        os.unlink('uniqueify.log')
        os.unlink('1x79-sf_uniqueify.mtz')
        os.unlink('1x79-sf_uniqueify.log')

        return
    
def testSuite():
    suite = unittest.TestSuite()
    suite.addTest(Test('testProcessReflectionFile'))
    suite.addTest(Test('testProcessReflectionFileNORFREE'))
    suite.addTest(Test('testProcessReflectionFileCIF'))
    return suite
    
#
# Run unit tests
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(testSuite())
