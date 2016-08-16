'''
Created on 2 Dec 2014

@author: jmht
'''

import logging
import os
import shutil
import sys

from iotbx import reflection_file_reader

import ample_util # Avoid circular dependencies
import exit_util
import cif_parser # Avoid circular dependencies

_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)

def del_column(file_name, column, overwrite=True):
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

def add_rfree(file_name,directory=None,overwrite=True):
    """Run uniqueify on mtz file to generate RFREE data column"""
    mtzUnique = ample_util.filename_append(file_name, "uniqueify", directory=directory)

    cmd = ['uniqueify', file_name, mtzUnique]
    logfile = os.path.join( os.getcwd(), "uniqueify.log" )
    retcode = ample_util.run_command(cmd, logfile=logfile)
    if retcode != 0:
        msg = "Error running command: {0}. Check the logfile: {1}".format(" ".join(cmd),logfile)
        _logger.critical(msg)
        raise RuntimeError, msg

    if overwrite:
        shutil.move(mtzUnique,file_name)
        return file_name
    else:
        return mtzUnique

def get_labels(file_name):
    """Return the F, FP and FREE column labels"""
    
    reflection_file = reflection_file_reader.any_reflection_file(file_name=file_name)
    if not reflection_file.file_type()=="ccp4_mtz":
        msg="File is not of type ccp4_mtz: {0}".format(file_name)
        logging.critical(msg)
        raise RuntimeError,msg
    
    content=reflection_file.file_content()
    ctypes=content.column_types()
    clabels=content.column_labels()
    ftype='F'
    if not ftype in ctypes:
        raise RuntimeError,"Cannot find any structure amplitudes in: {0}".format(file_name)
    F=clabels[ctypes.index(ftype)]
    
    # FP derived from F
    FP='SIG'+F
    if not FP in clabels:
        raise RuntimeError,"Cannot find label {0} in file: {1}".format(FP,file_name)
    
    FREE=_get_rfree(content)
    return F,FP,FREE

def get_rfree(file_name):
    """Return the Rfree label"""

    reflection_file = reflection_file_reader.any_reflection_file(file_name=file_name)
    if not reflection_file.file_type()=="ccp4_mtz":
        msg="File is not of type ccp4_mtz: {0}".format(file_name)
        logging.critical(msg)
        raise RuntimeError,msg
    
    return _get_rfree(reflection_file.file_content())

# def get_resolution(file_name):
#     if int(CCP4_VERSION[0]) >= 7:
#         hkl_info=clipper.HKL_info()
#         mtz_file=clipper.CCP4MTZfile()
#         mtz_file.open_read(file_name)
#         mtz_file.import_hkl_info(hkl_info)
# 
#         resolution =  "%.2lf" % hkl_info.resolution().limit()
#     else:
#         resolution = None
#     return resolution

def _get_rfree(content):
    rfree_label = None
    #print "GOT ",content.column_labels()
    for label in content.column_labels():
        if 'free' in label.lower():
            column = content.get_column(label=label)
            selection_valid = column.selection_valid()
            flags = column.extract_values()
            sel_0 = (flags == 0)
            # extract number of work/test reflections
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

def max_min_resolution(file_name):
    reflection_file = reflection_file_reader.any_reflection_file(file_name=file_name)
    if not reflection_file.file_type()=="ccp4_mtz":
        print("File is not of type ccp4_mtz: {0}".format( file_name ) )
        sys.exit(1)
    return reflection_file.file_content().max_min_resolution()

def to_hkl(mtz_file,hkl_file=None,directory=None,F=None,SIGF=None,FREE=None):
    
    if directory is None:
        directory=os.getcwd()
    
    if hkl_file is None:
        name=os.path.splitext(os.path.basename(mtz_file))[0]
        hkl_file=os.path.join(directory,name+".hkl")
        
    if F is None or SIGF is None or FREE is None:
        F,SIGF,FREE=get_labels(mtz_file)
        
    cmd=['mtz2various','HKLIN',mtz_file,'HKLOUT', hkl_file]
    logfile="mtz2various.log"
    stdin  = """LABIN FP={0} SIGFP={1} FREE={2}
OUTPUT SHELX
FSQUARED
END""".format(F,SIGF,FREE)
    
    ret = ample_util.run_command(cmd=cmd, logfile=logfile, directory=None, dolog=False, stdin=stdin)
    if not ret==0:
        raise RuntimeError,"Error converting {0} to HKL format - see log: {1}".format(mtz_file,logfile)
    else:
        os.unlink(logfile)
        
    return hkl_file

def processReflectionFile(amoptd):
    """Make sure we have a valid mtz file. If necessary convert a given cif file.
       Set the mtz variable in the given amoptd to the reflection file to use
       Return True if it all worked or raise an exception if it failed
    """

    # We've been given a sf_cif so convert to mtz
    if amoptd['sf_cif']:
        if not os.path.isfile(amoptd['sf_cif']):
            msg="Cannot find sf_cif file: {0}".format(amoptd['sf_cif'])
            exit_util.exit_error(msg)
        if not os.path.splitext(amoptd['sf_cif'])[1].lower() == ".cif":
            msg="Cif file extension is not .cif Please rename the file to give it a .cif extension."
            exit_util.exit_error(msg)

        cp = cif_parser.CifParser()
        amoptd['mtz'] = cp.sfcif2mtz(amoptd['sf_cif'])

    # Now have an mtz so check it's valid
    if not amoptd['mtz'] or not os.path.isfile( amoptd['mtz'] ):
        _logger.critical("Cannot find MTZ file: {0}".format( amoptd['mtz'] ) )
        sys.exit(1)

    # Get column label info
    reflection_file = reflection_file_reader.any_reflection_file(file_name=amoptd['mtz'])
    if not reflection_file.file_type() == "ccp4_mtz":
        _logger.critical("File is not of type ccp4_mtz: {0}".format( amoptd['mtz'] ) )
        sys.exit(1)
    
    # Read the file
    content = reflection_file.file_content()
    
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
        
    rfree=_get_rfree(content)
    if amoptd['FREE']:
        # Check is valid
        if not rfree or not rfree==amoptd['FREE']:
            _logger.critical("Given RFREE label {0} is not valid for mtz file: {0}".format( amoptd['FREE'], amoptd['mtz'] ) )
            sys.exit(1)
    else:
        # See if we can find a valid label in the file
        if not rfree:
            # Need to generate RFREE
            _logger.warning("Cannot find a valid FREE flag - running uniquefy to generate column with RFREE data." )
            amoptd['mtz'] = add_rfree( amoptd['mtz'], directory=amoptd['work_dir'],overwrite=False)

            # Check file and get new FREE flag
            rfree=get_rfree(amoptd['mtz'])
            if not rfree:
                _logger.critical("Cannot find valid rfree flag in mtz file {0} after running uniquiefy".format(amoptd['mtz']))
                sys.exit(1)
        amoptd['FREE']  = rfree
    
    # Output information to user and save to amoptd
    _logger.info("Using MTZ file: {0}".format(amoptd['mtz']))
    maxr, minr = content.max_min_resolution()
    amoptd['mtz_max_resolution'] = maxr
    amoptd['mtz_min_resolution'] = minr
    _logger.info("Resolution limits of MTZ file are: {0: > 6.3F} and {1: > 6.3F}".format(maxr, minr))

    return True
