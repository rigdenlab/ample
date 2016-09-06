#!/usr/bin/env ccp4-python

from collections import namedtuple
import datetime
import logging
import os

import phaser
import iotbx.pdb
import iotbx.mtz

from ample.util import pdb_edit
from ample.util import ample_util
from ample.util import mtz_util

_logger = logging.getLogger(__name__)

def _basis_str(origin_shift):
    """Return the string to get the sgtbx basis_op for the given origin shift"""
    bstr = ''
    for i, axis in enumerate(['x', 'y', 'z']):
        if origin_shift[i] == 0.0:
            s = axis
        else:
            s = '{0}{1:+F}'.format(axis,origin_shift[i])
        if i == 0:
            bstr = s
        else:
            bstr = bstr + ',' +s
    return bstr

def _check_mtz(mtz_file, labels=[]):
    """Return an iotbx.mtz object for the mtz_file making sure the file is valid and contains the given labels"""
    mtz = iotbx.mtz.object(file_name=mtz_file)
    if len(mtz.crystals()) > 2: raise RuntimeError("Cannot deal with multiple crystal in mtz")
    # # Assume the first crystal is always the base so we use the second
    if len(mtz.crystals()[1].datasets()) > 1: raise RuntimeError("Cannot deal with > 1 dataset in mtz")
    for label in labels:
        if not mtz.has_column(label):
            raise RuntimeError("Cannot find label :{0} in native mtz file: {1}".format(label, mtz_file))
    return mtz

def _get_miller_array_from_label(mtz_object, column_label):
    """Assuming only a single crystal/dataset, return the miller array for the given column label"""
    # Get the miller arrays for the columns we're interested in as dicts
    miller_dict = mtz_object.as_miller_arrays_dict()
    # Keys are 3-tuples where last item is array label
    try:
        col_key = [k for k in miller_dict.keys() if k[2] == column_label][0]
    except IndexError:
        raise RuntimeError("Cannot find column label: {0}".format(column_label))
    assert col_key
    return miller_dict[col_key]

def calc_phase_error_pdb(native_pdb, native_mtz, mr_mtz, f_label, sigf_label, fc_label='FC', cleanup=True, origin=None):
    """Phase error between native_pdb+native_mtz and mr_mtz
    """
    assert f_label and sigf_label,"Need f_label and sigf_label to be given!"
    native_mtz_phased = place_native_pdb(native_pdb, native_mtz, f_label, sigf_label)
    return calc_phase_error_mtz(native_mtz_phased, mr_mtz, f_label, sigf_label, fc_label, cleanup, origin)

def calc_phase_error_mtz(native_mtz_phased, mr_mtz, f_label=None, sigf_label=None, fc_label='FC', cleanup=True, origin=None):
    """Phase error between native_pdb+native_mtz and mr_mtz
    """
    if origin:
        # if we are given an origin shift we can calculate the phase difference directly with cctbx
        change_of_hand = False # hard-wired as we don't have this information
        origin_shift = origin

        native_object = _check_mtz(native_mtz_phased, labels = [fc_label])
        mr_object = _check_mtz(mr_mtz, labels = [fc_label])
        
        native_fc_miller_array = _get_miller_array_from_label(native_object, fc_label)
        mr_fc_miller_array = _get_miller_array_from_label(mr_object, fc_label)
        
        before_origin = native_fc_miller_array.mean_phase_error(mr_fc_miller_array)
        basis_str = _basis_str(origin)
        after_origin = native_fc_miller_array.mean_phase_error(mr_fc_miller_array.change_basis(basis_str))
    else:
        assert f_label and sigf_label,"Need f_label and sigf_label to be given!"
        # we use cphasematch to calculate the phase difference and origin shift
        r = make_merged_mtz(native_mtz_phased, mr_mtz, f_label)
        before_origin, after_origin, change_of_hand, origin_shift = run_cphasematch(r.merged_mtz,
                                                                                    f_label,
                                                                                    sigf_label,
                                                                                    fc_label=r.fc_label,
                                                                                    phifc_label=r.phifc_label,
                                                                                    fcalc_label=r.fcalc_label,
                                                                                    phifcalc_label=r.phifcalc_label)
        if cleanup:
            os.unlink(r.merged_mtz)
    return before_origin, after_origin, change_of_hand, origin_shift

def make_merged_mtz(native_mtz, mr_mtz, f_label, fc_label='FC', fcalc_label='FCALC'):
    """Create MTZ file with F from native_mtz and calculated phases from native_mtz and mr_mtz to enable phaser error calc by cphasematch"""
    # Check mtz files have required columns
    native_object = _check_mtz(native_mtz, labels = [f_label])
    mr_object = _check_mtz(mr_mtz, labels = [fc_label])

    merged_object = iotbx.mtz.object()
    merged_object.set_title("Calculated phases from {0} and {1}".format(native_mtz, mr_mtz))
    merged_object.add_history(line="Created by ample make_merged_mtz on: {0}".format(datetime.datetime.now()))
    merged_object.add_history(line="Created from: {0} and {1}".format(native_mtz, mr_mtz))
    
    f_miller_array = _get_miller_array_from_label(native_object, f_label)
    fc_miller_array = _get_miller_array_from_label(native_object, fc_label)
    
    unit_cell = fc_miller_array.unit_cell()
    merged_object.set_space_group_info(fc_miller_array.space_group_info())
    merged_object.set_hkl_base(unit_cell)
    crystal = merged_object.add_crystal(
                                        name="cphasematch_crystal",
                                        project_name="cphasematch_project",
                                        unit_cell = unit_cell
                                        )
    dataset = crystal.add_dataset(name = "test_dataset",
                                  wavelength=1) # FIX
    # The labels in add_miller_array automatically controlled by the label_decorator which adds 'PHI' 
    # as a prefix to the appropriate columns we therefore either need to make our own decorator or just stay with
    # using the PHI prefix
    dataset.add_miller_array(
                             miller_array=f_miller_array,
                             column_root_label=f_label,
                             )
    dataset.add_miller_array(
                             miller_array=fc_miller_array,
                             column_root_label=fc_label,
                             )
    
    # Now add the other calculated data
    mr_fc_miller_array = _get_miller_array_from_label(mr_object, fc_label)
    dataset.add_miller_array(
                             miller_array=mr_fc_miller_array,
                             column_root_label=fcalc_label,
                             )
    
    # as mtz dataset and then change labels?
    #dataset.add_column(label=fcalc_label, type="F").set_values(values=mr_object.get_column(fc_label).extract_values())
    #dataset.add_column(label=phicalc_label, type="P").set_values(values=mr_object.get_column(phic_label).extract_values())

    # Write out the file
    name1 = os.path.splitext(os.path.basename(native_mtz))[0]
    name2 = os.path.splitext(os.path.basename(mr_mtz))[0]
    merged_mtz = "{0}_{1}.mtz".format(name1, name2)
    merged_object.write(merged_mtz)
    
    phifc_label = 'PHI' + fc_label
    phifcalc_label = 'PHI' + fcalc_label
    results = namedtuple('results', ['merged_mtz','f_label', 'fc_label' ,'phifc_label','fcalc_label','phifcalc_label'])
    return results(merged_mtz, f_label, fc_label, phifc_label, fcalc_label, phifcalc_label)

# def make_merged_mtz(native_mtz, mr_mtz, fc_label = 'FC', phic_label = 'PHIC', fcalc_label='FCalc', phicalc_label='PHICalc'):
#     """Add fc_label and phic_label columns from mr_mtz to native_mtz and write out as a new mtz file.
#     """
#     # Add calculated phases from mr_mtz to those in native_mtz
#     native_iotbx = iotbx.mtz.object(file_name=native_mtz)
#     mr_iotbx = iotbx.mtz.object(file_name=mr_mtz)
#     
#     # Assume there is only one real crystal and one dataset
#     dataset = native_iotbx.crystals()[1].datasets()[0]
#     
#     dataset.add_column(label=fcalc_label, type="F").set_values(values=mr_iotbx.get_column(fc_label).extract_values())
#     dataset.add_column(label=phicalc_label, type="P").set_values(values=mr_iotbx.get_column(phic_label).extract_values())
#     
#     # Write out the file
#     name1 = os.path.splitext(os.path.basename(native_mtz))[0]
#     name2 = os.path.splitext(os.path.basename(mr_mtz))[0]
#     merged_mtz = "{0}_{1}.mtz".format(name1, name2)
#     native_iotbx.write(merged_mtz)
#     return merged_mtz

def place_native_pdb(native_pdb, native_mtz, f_label, sigf_label, cleanup=True):
    """Place native_pdb into data from native_mtz using phaser with f_label and sigf_label"""
    # get crystal info from pdb
    sym_obj = iotbx.pdb.pdb_input(file_name=native_pdb).crystal_symmetry()
    hall_symbol = sym_obj.space_group().type().hall_symbol()
    unit_cell = sym_obj.unit_cell().parameters()
    molecular_weight = pdb_edit.molecular_weight(native_pdb)
    
    # Get data from MTZ file
    mtz_object = _check_mtz(native_mtz, labels = [f_label, sigf_label])
    miller_indices = mtz_object.extract_miller_indices()
    fp_obs = mtz_object.get_column(f_label).extract_values().as_double()
    sigfp_obs = mtz_object.get_column(sigf_label).extract_values().as_double()
    
    # Create root_name for all files from native + mr_mtz
    name1 = os.path.splitext(os.path.basename(native_pdb))[0]
    name2 = os.path.splitext(os.path.basename(native_mtz))[0]
    root_name = "{0}_{1}".format(name1, name2)
    
    # Run phaser to position native
    phaser_input = phaser.InputMR_AUTO()
    phaser_input.setSPAC_HALL(hall_symbol)
    phaser_input.setCELL6(unit_cell)
    phaser_input.setREFL_F_SIGF(miller_indices, fp_obs, sigfp_obs)
    # Make sure labels same as native_mtz
    phaser_input.setLABI_F_SIGF(f_label,sigf_label)
    phaser_input.setROOT(root_name)
    phaser_input.addENSE_PDB_ID("native", native_pdb, 0.0)
    phaser_input.addCOMP_PROT_MW_NUM(molecular_weight, 1)
    phaser_input.addSEAR_ENSE_NUM("native", 1)
    phaser_input.setMUTE(True)
    phaser_run = phaser.runMR_AUTO(phaser_input)
    del phaser_input
    if not phaser_run.Success():
        raise RuntimeError("PHASER failed: {0} : {1}".format(phaser_run.ErrorName(), phaser_run.ErrorMessage()))
        
    if phaser_run.foundSolutions():
        #print "Phaser has found MR solutions"
        #print "Top LLG = %f" % phaser_run.getTopLLG()
        #print "Top PDB file = %s" % phaser_run.getTopPdbFile()
        #print "Top MTZ file = %s" % phaser_run.getTopMtzFile()
        native_placed_mtz = phaser_run.getTopMtzFile()
    else:
        raise RuntimeError("PHASER could not place native_pdb")
    
    if cleanup:
        os.unlink(phaser_run.getTopPdbFile())
    return native_placed_mtz

def parse_cphasematch_log(logfile):
    """Return phase error before_origin& after_origin from cphasematch logfile"""
    before_origin = None
    after_origin = None
    change_of_hand = False
    origin_shift = False
    with open(logfile, 'r') as f:
        for line in f:
            if line.startswith(' Mean phase error before origin fixing:'):
                before_origin = float(line.split()[-1])
            elif line.startswith(' Mean phase error after  origin fixing:'):
                after_origin = float(line.split()[-1])
            elif line.startswith(' Change of hand'):
                if line.strip().split()[4] == 'YES':
                    change_of_hand = True
            elif line.startswith(' Change of origin:'):
                # Change of origin: uvw = (         0,         0,         0)
                fields = line.strip().split()
                x = float(fields[6][:-1]) # strip trailing comma
                y = float(fields[7][:-1]) # strip trailing comma
                z = float(fields[8][:-1]) # strip trailing bracket
                origin_shift = [x,y,z]
                
    return before_origin, after_origin, change_of_hand, origin_shift

def run_cphasematch(merged_mtz,
                    f_label,
                    sigf_label,
                    fc_label='FC',
                    phifc_label='PHIFC',
                    fcalc_label='FCALC',
                    phifcalc_label='PHIFCALC',
                    resolution_bins=12,
                    cleanup=True):
    """run cphasematch to get phase error"""
    
    assert merged_mtz and f_label and sigf_label
    argd = { 'merged_mtz' : merged_mtz,
             'f_label' : f_label,
             'sigf_label' : sigf_label,
             'fc_label' : fc_label,
             'phic_label' : phifc_label,
             'fcalc_label' : fcalc_label,
             'phicalc_label' : phifcalc_label,
             'resolution_bins' : resolution_bins }
    
    stdin = """
mtzin {merged_mtz}
colin-fo /*/*/[{f_label},{sigf_label}]
colin-fc-1 /*/*/[{fc_label},{phic_label}]
colin-fc-2 /*/*/[{fcalc_label},{phicalc_label}]
resolution-bins {resolution_bins}
""".format(**argd)

    logfile = os.path.abspath("cphasematch.log")
    cmd= [ 'cphasematch',
          "-stdin" ]
    
    retcode = ample_util.run_command(cmd=cmd, stdin=stdin, logfile=logfile)
    if retcode != 0: raise RuntimeError, "Error running command: {0}\nCheck logfile: {1}".format( " ".join(cmd), logfile )
    
    before_origin, after_origin, change_of_hand, origin_shift = parse_cphasematch_log(logfile)
    
    if cleanup:os.unlink(logfile)
    
    return before_origin, after_origin, change_of_hand, origin_shift


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    import argparse
    parser = argparse.ArgumentParser(description='Run cphasematch on two mtz files', prefix_chars="-")

    #group = parser.add_argument_group()
    parser.add_argument('-native_mtz',
                       help="Input native MTZ file",
                       required=True)
    parser.add_argument('-mr_mtz',
                       help="Input MTZ file from Molecular Replacement",
                       required=True)
    parser.add_argument('-native_pdb',
                       help="Input native PDB file",
                       required=False)
    parser.add_argument('-f_label',
                       help="F label from native MTZ file",
                       required=False)
    parser.add_argument('-sigf_label',
                       help="SIGF label from native MTZ file",
                       required=False)

    args = parser.parse_args()
    
    # Extract F & SIGF from file if not give on command-line
    if args.f_label and args.sigf_label:
        f_label, sigf_label = args.f_label, args.sigf_label
    else:
        f_label, sigf_label, _ = mtz_util.get_labels(args.native_mtz)
    _logger.info("Using F, SIGF labels from mtz file: {0}, {1} ".format(f_label, sigf_label))
    
    if args.native_pdb:
        before_origin, after_origin, change_of_hand, origin_shift = calc_phase_error_pdb(args.native_pdb,
                                                                                         args.native_mtz,
                                                                                         args.mr_mtz,
                                                                                         f_label,
                                                                                         sigf_label)
    else:
        before_origin, after_origin, change_of_hand, origin_shift = calc_phase_error_mtz(args.native_mtz,
                                                                                         args.mr_mtz,
                                                                                         f_label,
                                                                                         sigf_label)
    
    print "Calculated phase error: {0}".format(after_origin)
    print "Calculated origin shift: {0}".format(origin_shift)
