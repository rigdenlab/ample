#!/usr/bin/env ccp4-python

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
            raise RuntimeError("Cannot find label: {0} in mtz file: {1}".format(label, mtz_file))
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

def calc_phase_error_pdb(native_pdb, native_mtz, mr_mtz, f_label, sigf_label, native_phase_labels=['FC', 'PHIC'], mr_phase_labels=['FC', 'PHIC'], cleanup=True, origin=None):
    """Phase error between native_pdb+native_mtz and mr_mtz
    """
    assert f_label and sigf_label,"Need f_label and sigf_label to be given!"
    native_mtz_phased = place_native_pdb(native_pdb, native_mtz, f_label, sigf_label)
    return calc_phase_error_mtz(native_mtz_phased,
                                mr_mtz,
                                f_label,
                                sigf_label,
                                native_phase_labels=native_phase_labels,
                                mr_phase_labels=mr_phase_labels,
                                cleanup=cleanup,
                                origin=origin)

def calc_phase_error_mtz(native_mtz_phased, mr_mtz, f_label=None, sigf_label=None, native_phase_labels=['FC', 'PHIC'], mr_phase_labels=['FC', 'PHIC'], cleanup=True, origin=None):
    """Phase error between native_pdb+native_mtz and mr_mtz
    """
    
    if origin:
        assert False
        # if we are given an origin shift we can calculate the phase difference directly with cctbx
        change_of_hand = False # hard-wired as we don't have this information
        origin_shift = origin

        native_object = _check_mtz(native_mtz_phased, labels=native_phase_labels)
        mr_object = _check_mtz(mr_mtz, labels = mr_phase_labels)
        
        native_fc_miller_array = _get_miller_array_from_label(native_object, native_phase_labels[0])
        mr_fc_miller_array = _get_miller_array_from_label(mr_object, mr_phase_labels[0])
        
        before_origin = native_fc_miller_array.mean_phase_error(mr_fc_miller_array)
        basis_str = _basis_str(origin)
        after_origin = native_fc_miller_array.mean_phase_error(mr_fc_miller_array.change_basis(basis_str))
    else:
        assert f_label and sigf_label,"Need f_label and sigf_label to be given!"
        # we use cphasematch to calculate the phase difference and origin shift
        merged_mtz, labels = merge_mtz(native_mtz_phased, [f_label, sigf_label] + native_phase_labels, mr_mtz, mr_phase_labels)
        assert len(labels)==6
        before_origin, after_origin, change_of_hand, origin_shift = run_cphasematch(merged_mtz,
                                                                                    labels[0:2],
                                                                                    native_phase_labels=labels[2:4],
                                                                                    mr_phase_labels=labels[4:])
        if cleanup: os.unlink(merged_mtz)
    return before_origin, after_origin, change_of_hand, origin_shift

def merge_mtz_cctbx(mtz1_path, mtz1_labels, mtz2_path, mtz2_labels):
    """TODO!
    Create MTZ file with F from native_mtz and calculated phases from native_mtz and mr_mtz to enable phaser error calc by cphasematch"""
    # Check mtz files have required columns
    mtz1 = _check_mtz(mtz1_path, labels = mtz1_labels)
    mtz2 = _check_mtz(mtz2_path, labels = mtz2_labels)

    merged_object = iotbx.mtz.object()
    merged_object.set_title("Calculated phases from {0} and {1}".format(mtz1_path, mtz2_path))
    merged_object.add_history(line="Created by ample merged_mtz on: {0}".format(datetime.datetime.now()))
    merged_object.add_history(line="Created from: {0} and {1}".format(mtz1_path, mtz2_path))
    
    miller_array1 = _get_miller_array_from_label(mtz1, mtz1_labels[0])
    
    unit_cell = miller_array1.unit_cell()
    merged_object.set_space_group_info(miller_array1.space_group_info())
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
    
    # Loop through the labels in each mtz
    labels = []
    for i, mtz in enumerate([mtz1, mtz2]):
        for label in [mtz1_labels, mtz2_labels][i]:
            miller_array = _get_miller_array_from_label(mtz, label)
            dataset.add_miller_array(
                                     miller_array=miller_array,
                                     column_root_label=label,
                                     )
            labels.append(label)
    
    # as mtz dataset and then change labels?
    #dataset.add_column(label=fcalc_label, type="F").set_values(values=mr_object.get_column(fc_label).extract_values())
    #dataset.add_column(label=phicalc_label, type="P").set_values(values=mr_object.get_column(phic_label).extract_values())

    # Write out the file
    name1 = os.path.splitext(os.path.basename(mtz1_path))[0]
    name2 = os.path.splitext(os.path.basename(mtz2_path))[0]
    merged_mtz = "{0}_{1}.mtz".format(name1, name2)
    merged_object.write(merged_mtz)
    
    return os.path.abspath(merged_mtz), labels

def merge_mtz(mtz1_path, mtz1_labels, mtz2_path, mtz2_labels):
    """Create MTZ file with columns from the given mtz files and mtz labels in each file"""
    
    # Can't have any duplicates in file labels
    assert len(mtz1_labels) == len(set(mtz1_labels)),"Duplicate labels in mtz1_labels"
    assert len(mtz2_labels) == len(set(mtz2_labels)),"Duplicate labels in mtz2_labels"

    name1 = os.path.splitext(os.path.basename(mtz1_path))[0]
    name2 = os.path.splitext(os.path.basename(mtz2_path))[0]
    merged_mtz = os.path.abspath("{0}_{1}.mtz".format(name1, name2))
    
    cmd = [ 'cad', 'hklin1', mtz1_path, 'hklin2', mtz2_path, 'hklout', merged_mtz ]

    # See if any labels are duplicate and need to be renamed
    rename = [] # List of (File_number, file_label_idx, orig_label, renamed_label)
    labels = []
    for i, mtz in enumerate([mtz1_path, mtz2_path]):
        for j, label in enumerate([mtz1_labels, mtz2_labels][i]):
            if label in labels:
                newlabel = label + str(i+1)
                rename.append((i+1,j+1, label, newlabel))
            else:
                newlabel = label
                rename.append((i+1,j+1, label, None))
            assert newlabel not in labels, "Too many duplicate label names: {0}".format(newlabel)
            labels.append(newlabel)

    # Build up the list of which labels to extract from which files
    stdin = ""
    last_fileno = None
    for fileno, labelno, orig_label, rename_label in rename:
        if fileno != last_fileno:
            if last_fileno is not None:
                stdin += '\n' # Need to terminate the line
            stdin += "LABIN FILE {0}".format(fileno)
            last_fileno = fileno
        stdin += " E{0}={1}".format(labelno, orig_label)
    stdin += '\n' # Need to terminate the line
    
    # Do any renaming for duplicate labels
    last_fileno = None
    for i, (fileno, label_idx, orig_label, rename_label) in enumerate(rename):
        if rename_label is not None:
            if last_fileno != fileno:
                stdin += 'LABOUT FILE_NUMBER {0}'.format(fileno)
                if last_fileno is not None:
                    # for anything other then then first, we need to terminate this block
                    stdin += '\n'
                last_fileno = fileno
            if fileno == last_fileno:
                stdin += ' E{0}={1}'.format(label_idx,rename_label)
    
    if fileno is not None:
        stdin += '\n' # Add last linebreak as we have added a rename clause
    
    logfile = os.path.abspath("cad.log")
    retcode = ample_util.run_command(cmd=cmd, stdin=stdin, logfile=logfile)
    if retcode != 0: raise RuntimeError, "Error running command: {0}\nCheck logfile: {1}".format( " ".join(cmd), logfile )
    
    return os.path.abspath(merged_mtz), labels

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
                    f_sigf_labels,
                    native_phase_labels,
                    mr_phase_labels,
                    resolution_bins=12,
                    cleanup=True):
    """run cphasematch to get phase error"""
    
    argd = { 'merged_mtz' : merged_mtz,
             'f_label' : f_sigf_labels[0],
             'sigf_label' : f_sigf_labels[1],
             'native_phase_label1' : native_phase_labels[0],
             'native_phase_label2' : native_phase_labels[1],
             'mr_phase_label1' : mr_phase_labels[0],
             'mr_phase_label2' : mr_phase_labels[1],
             'resolution_bins' : resolution_bins }
    
    stdin = """
mtzin {merged_mtz}
colin-fo /*/*/[{f_label},{sigf_label}]
colin-fc-1 /*/*/[{native_phase_label1},{native_phase_label2}]
colin-fc-2 /*/*/[{mr_phase_label1},{mr_phase_label2}]
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
