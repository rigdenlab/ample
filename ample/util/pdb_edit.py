"""Useful manipulations on PDB files"""

from __future__ import division, print_function

__author__ = "Adam Simpkin, Jens Thomas & Felix Simkovic"
__date__ = "21 Apr 2017"
__version__ = "2.0"

import glob
import logging
import numpy as np
import os
import string
import warnings

from cctbx.array_family import flex

import iotbx.pdb
import iotbx.pdb.amino_acid_codes

import ample_util
import chemistry
import pdb_model
import residue_map
import sequence_util

three2one = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter
one2three = iotbx.pdb.amino_acid_codes.three_letter_given_one_letter

logger = logging.getLogger()


def _cache(pdbin):
    """Cache the PDB input file"""
    pdb_input = iotbx.pdb.pdb_input(file_name=pdbin)
    crystal_symmetry = pdb_input.crystal_symmetry()
    hierarchy = pdb_input.construct_hierarchy()
    return pdb_input, hierarchy, crystal_symmetry


def _first_chain_only(h):
    """Remove everything from hierarchy but the first chain"""
    for i, m in enumerate(h.models()):
        if i != 0:
            h.remove_model(m)
    m = h.models()[0]
    for i, c in enumerate(m.chains()):
        if i != 0:
            m.remove_chain(c)


def _most_prob(hierarchy, always_keep_one=True, reset_i_seq=True):
    """Remove alternate conforms from a hierarchy"""
    for m in hierarchy.models():
        for c in m.chains():
            for rg in c.residue_groups():
                if rg.have_conformers():
                    atom_groups = rg.atom_groups()
                    if always_keep_one:
                        if len(atom_groups) == 1 and atom_groups[0].altloc == '':
                            continue
                        atom_groups_and_occupancies = []
                        for ag in atom_groups:
                            if ag.altloc == '':
                                continue
                            mean_occ = flex.mean(ag.atoms().extract_occ())
                            atom_groups_and_occupancies.append((ag, mean_occ))
                        atom_groups_and_occupancies.sort(lambda a, b: cmp(b[1], a[1]))
                        for atom_group, occ in atom_groups_and_occupancies[1:]:
                            rg.remove_atom_group(atom_group=atom_group)
                        single_conf, occ = atom_groups_and_occupancies[0]
                        single_conf.altloc = ''
                    elif rg.have_conformers():
                        for ag in atom_groups:
                            if ag.altloc not in ["", "A"]:
                                rg.remove_atom_group(atom_group=ag)
                            else:
                                atom_group.altloc = ""
                    # Essential so we don't split residue groups by atom groups
                    rg.merge_atom_groups(atom_groups[0], atom_groups[1])
                if rg.atom_groups_size() == 0:
                    c.remove_residue_group(rg)
            if c.residue_groups_size() == 0:
                m.remove_chain(c)
        if m.chains_size() == 0:
            hierarchy.remove_model(m)
    # Assign occupancies of 1.0 to all remaining atoms
    new_occ = flex.double(hierarchy.atoms().size(), 1.0)
    hierarchy.atoms().set_occ(new_occ)
    # Reset the atom i_seq entry
    if reset_i_seq:
        hierarchy.atoms().reset_i_seq()


def _natm_nres_mw(hierarchy, first=False):
    """Function to extract the number of atoms, number of residues and molecular weight
    from a PDB structure
    """
    # Define storage variables
    natm, nres, mw = 0, 0, 0
    hydrogen_atoms, other_atoms = 0, 0
    water_atoms, water_hydrogen_atoms = 0, 0
    mw = 0

    # Pick chain
    if first:
        _first_chain_only(hierarchy)

    # Collect all the data using the hierarchy
    for m in hierarchy.models():
        for c in m.chains():
            for rg in c.residue_groups():
                resseq = None
                for ag in rg.atom_groups():
                    if ag.resname in three2one and resseq != rg.resseq:
                        nres += 1
                        resseq = rg.resseq
                        hydrogen_atoms += chemistry.atomic_composition[ag.resname].H
                    for atom in ag.atoms():
                        if ag.resname.strip() == 'HOH' or ag.resname.strip() == 'WAT':
                            water_hydrogen_atoms += (2.0 * atom.occ)
                            water_atoms += (1.0 * atom.occ)
                        else:
                            other_atoms += atom.occ
                            # Be carful, models might not have the last element column
                            if atom.element.strip():
                                aname = atom.element.strip()
                            else:
                                aname = atom.name.strip()
                                aname = aname.translate(None, string.digits)[0]
                            mw += chemistry.periodic_table(aname).weight() * atom.occ

    mw += hydrogen_atoms * chemistry.periodic_table('H').weight()

    # Compute the number of atoms and number of residues
    natm = int(other_atoms + hydrogen_atoms + water_atoms + water_hydrogen_atoms - 0.5)

    return natm, nres, mw


def _parse_modres(modres_text):
    """
COLUMNS        DATA TYPE     FIELD       DEFINITION
--------------------------------------------------------------------------------
 1 -  6        Record name   "MODRES"
 8 - 11        IDcode        idCode      ID code of this entry.
13 - 15        Residue name  resName     Residue name used in this entry.
17             Character     chainID     Chain identifier.
19 - 22        Integer       seqNum      Sequence number.
23             AChar         iCode       Insertion code.
25 - 27        Residue name  stdRes      Standard residue name.
30 - 70        String        comment     Description of the residue modification.
    """
    modres = []
    for line in modres_text:
        assert line[0:6] == "MODRES", "Line did not begin with an MODRES record!: {0}".format(line)
        id_code = line[7:11]
        resname = line[12:15].strip()
        # Use for all so None means an empty field
        chain_id = line[16] if line[16].strip() else ""
        seq_num = int(line[18:22])
        i_code = line[22] if line[22].strip() else ""
        std_res = line[24:27].strip()
        comment = line[29:70].strip() if line[29:70].strip() else ""
        modres.append([id_code, resname, chain_id, seq_num, i_code, std_res, comment])
    return modres


def _rename_chains(hierarchy, table):
    """Rename all chains in a hierarchy using the provided conversion table"""
    for chain in hierarchy.chains():
        if chain.id in table:
            chain.id = table[chain.id]


def _renumber(hierarchy, start):
    """Renumber the residue sequence"""
    for model in hierarchy.models():
        for chain in model.chains():
            for idx, residue_group in enumerate(chain.residue_groups()):
                residue_group.resseq = idx + start


def _renumber_residues_gaps(hierarchy, gaps, start):
    """Renumber the residue sequence respecting gaps"""
    for model in hierarchy.models():
        for chain in model.chains():
            resseq = 0
            for idx, is_gap in enumerate(gaps):
                if is_gap:
                    continue
                residue_group = chain.residue_groups()[resseq]
                residue_group.resseq = idx + start
                resseq += 1


def _resseq(hierarchy):
    """Extract the sequence of residues from a pdb file"""
    chain2data = _sequence_data(hierarchy)
    return dict((k, chain2data[k][1]) for k in chain2data.keys())


def _rog_side_chain_treatment(hierarchy, scores, del_orange):
    """Remove side chains based on the Red/Orange/Green system"""
    def _remove(rg):
        for ag in rg.atom_groups():
            for atom in ag.atoms():
                if atom.name.strip() not in ['N', 'CA', 'C', 'O', 'CB']:
                    ag.remove_atom(atom=atom)

    for model in hierarchy.models():
        for chain in model.chains():
            for rg in chain.residue_groups():
                if not del_orange:
                    # Remove just the red scoring side chains
                    res = rg.resseq_as_int()
                    for i, j in scores:
                        if i == res and j == 'red':
                            _remove(rg)
                else:
                    # Only keep the green scoring side chains
                    res = rg.resseq_as_int()
                    for i, j in scores:
                        if i == res and j != 'green':
                            _remove(rg)


def _save(pdbout, hierarchy, crystal_symmetry=None, remarks=[]):
    """Save the CCTBX hierarchy to a file"""
    with open(pdbout, 'w') as f_out:
        for remark in remarks:
            f_out.write("REMARK %s" % remark + os.linesep)
        f_out.write(hierarchy.as_pdb_string(
            anisou=False, write_scale_records=True, crystal_symmetry=crystal_symmetry
        ))


def _select(hierarchy, construct):
    """Select atoms from hierarchy using a selection string"""
    selection = hierarchy.atom_selection_cache().selection(construct)
    return hierarchy.select(selection)


def _sequence(hierarchy):
    """Extract the sequence of residues from a pdb file."""
    chain2data = _sequence_data(hierarchy)
    return dict((k, chain2data[k][0]) for k in chain2data.keys())


def _sequence1(hierarchy):
    """Return sequence of the first chain"""
    d = _sequence(hierarchy)
    return d[sorted(d.keys())[0]]


def _sequence_data(hierarchy):
    """Extract the sequence of residues and resseqs from a pdb file."""
    chain2data = {}
    for chain in set(hierarchy.models()[0].chains()):  # only the first model
        if not chain.is_protein():
            continue
        got = False
        seq = ""
        resseq = []
        for residue in chain.conformers()[0].residues():  # Just look at the first conformer
            # See if any of the atoms are non-hetero - if so we add this residue
            if any([not atom.hetero for atom in residue.atoms()]):
                got = True
                seq += three2one[residue.resname]
                # resseq.append(int(residue.resseq.strip()))
                resseq.append(residue.resseq_as_int())
        if got:
            chain2data[chain.id] = (seq, resseq)
    return chain2data


def _strip(hierarchy, reset_i_seq=True, hetatm=False, hydrogen=False, atom_types=[]):
    """Remove all hetatoms from pdbfile"""
    for m in hierarchy.models():
        for c in m.chains():
            for rg in c.residue_groups():
                for ag in rg.atom_groups():
                    for a in ag.atoms():
                        if hydrogen and a.element.strip().upper() in ["H"]:
                            ag.remove_atom(a)
                        elif len(atom_types) > 0 and a.name.strip().upper() in atom_types:
                            ag.remove_atom(a)
                        elif hetatm and a.hetero:
                            ag.remove_atom(a)
                    if ag.atoms_size() == 0:
                        rg.remove_atom_group(ag)
                if rg.atom_groups_size() == 0:
                    c.remove_residue_group(rg)
            if c.residue_groups_size() == 0:
                m.remove_chain(c)
        if m.chains_size() == 0:
            hierarchy.remove_model(m)
    if reset_i_seq:
      hierarchy.atoms().reset_i_seq()


def _translate(hierarchy, vector):
    """Translate all atoms in the hierarchy by the provided vector"""
    for atom in hierarchy.atoms():
        atom.set_xyz(np.asarray(atom.xyz) + vector)
        

def backbone(pdbin, pdbout):
    """Only output backbone atoms
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    
    """
    _, hierarchy, symmetry = _cache(pdbin)
    hierarchy = _select(hierarchy, "name n or name ca or name c or name o or name cb")
    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=['Original file: {0}'.format(pdbin)])


def calpha_only(pdbin, pdbout):
    """Strip PDB to c-alphas only
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    
    """
    _, hierarchy, symmetry = _cache(pdbin)
    hierarchy = _select(hierarchy, "name ca")
    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=['Original file: {0}'.format(pdbin)])


def check_pdb_directory(directory, single=True, allsame=True, sequence=None):
    """Check a directory of structure files

    Parameters
    ----------
    directory : str
       The path to the structure file directory
    single : bool, optional
       Each file contains a single model [default: True]
    allsame : bool, optional
       All structures contain the same sequence [default: True]
    sequence : str, optional
       The sequence of all models

    Returns
    -------
    bool
       A status describing if all structure files are okay

    """
    logger.info("Checking pdbs in directory: %s", directory)
    if os.path.isdir(directory):
        models = glob.glob(os.path.join(directory, "*.pdb"))
        if len(models) > 0 and not (single or sequence or allsame):
            return True
        elif len(models) > 0:
            return check_pdbs(models, sequence=sequence, single=single, allsame=allsame)
        else:
            logger.critical("Cannot find any pdb files in directory: %s", directory)
    else:
        logger.critical("Cannot find directory: %s", directory)
    return False


def check_pdbs(models, single=True, allsame=True, sequence=None):
    """Check a set of structure files

    Parameters
    ----------
    models : list
       A list of paths to structure files
    single : bool, optional
       Each file contains a single model [default: True]
    allsame : bool, optional
       All structures contain the same sequence [default: True]
    sequence : str, optional
       The sequence of all models

    Returns
    -------
    bool
       A status describing if all structure files are okay

    """
    # Get sequence from first model
    if allsame and not sequence:
        try:
            _, h, _ = _cache(models[0])
        except Exception as e:
            s = "*** ERROR reading sequence from first pdb: {0}\n{1}".format(models[0], e)
            logger.critical(s)
            return False
        sequence = _sequence1(h)  # only one model/chain

    # Store info as simple array - errors, multi, no_protein, sequence_err
    summary = np.zeros((0, 5), dtype=np.uint8)
    for idx, pdb in enumerate(models):
        entry = np.zeros((1, 5), dtype=np.uint8)
        entry[0][0] = idx
        try:
            _, h, _ = _cache(models[0])
        except Exception:
            entry[0][1] = 1
            continue
        if single and h.models_size() != 1 and h.models()[0].chains_size() != 1:
            entry[0][2] = 1
        elif single and not h.models()[0].chains()[0].is_protein():
            entry[0][3] = 1
        elif sequence and sequence != _sequence(h).values()[0]:
            entry[0][4] = 1
        summary = np.concatenate((summary, entry), axis=0)

    # The summary table has no error messages (indicated by 0s)
    if np.count_nonzero(summary[:, 1:]) == 0:
        logger.info("check_pdb_directory - pdb files all seem valid")
        return True
    # The summary table has error messages (indicated by non-0s)
    else:
        s = "\n*** ERROR ***\n"
        if np.count_nonzero(summary[:, 1]) != 0:
            s += "The following pdb files have errors:\n\n"
            for idx in np.nonzero(summary[:, 1])[0]:
                s += "\t{0}\n".format(models[idx])
        elif np.count_nonzero(summary[:, 2]) != 0:
            s += "The following pdb files have more than one chain:\n\n"
            for idx in np.nonzero(summary[:, 2])[0]:
                s += "\t{0}\n".format(models[idx])
        elif np.count_nonzero(summary[:, 3]) != 0:
            s += "The following pdb files do not appear to contain any protein:\n\n"
            for idx in np.nonzero(summary[:, 3])[0]:
                s += "\t{0}\n".format(models[idx])
        elif np.count_nonzero(summary[:, 4]) != 0:
            s += "The following pdb files have diff sequences from the ref sequence: {0}\n\n".format(sequence)
            for idx in np.nonzero(summary[:, 4])[0]:
                s += "\t{0}\n".format(models[idx])
        logger.critical(s)
        return False


def extract_chain(pdbin, pdbout, chain_id, new_chain_id=None, c_alpha=False, renumber=False):
    """Extract a chain from the input PDB
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    chain_id : str
       The chain to extract
    new_chain_id : str, optional
       The chain ID to rename ``chain_id`` to
    c_alpha : bool, optional
       Strip chain residues back to c-alpha atoms [defaut: False]
    renumber : bool, optional
       Renumber the chain [default: False]
    
    """
    _, hierarchy, symmetry = _cache(pdbin)

    sel_string = "chain %s and not hetero" % chain_id
    if c_alpha:
        sel_string += " and name ca" 
    hierarchy = _select(hierarchy, sel_string)

    if new_chain_id:
        for model in hierarchy.models():
            for chain in model.chains():
                chain.id = new_chain_id

    if renumber:
        _renumber(hierarchy, 1)
        hierarchy.atoms().reset_serial()

    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=["Original file: %s" % pdbin])


def extract_model(pdbin, pdbout, model_id):
    """Extract a model from the input PDB

    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    model_id : str
       The model to extract

    """
    _, hierarchy, symmetry = _cache(pdbin)
    hierarchy = _select(hierarchy, "model {0}".format(model_id))
    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=['Original file: {0}'.format(pdbin)])


def extract_resSeq(pdbin, chain_id=None):
    """Extract a residue numbers of the input PDB

    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    chain_id : str, optional
       The chain to extract

    Returns
    -------
    list
       A list of the residue numbers

    """
    _, hierarchy, _ = _cache(pdbin)
    pdb_input = iotbx.pdb.pdb_input(file_name=pdbin)
    hierarchy = pdb_input.construct_hierarchy()

    # We only consider chains from the first model in the hierarchy
    chains = {}
    for chain in hierarchy.models()[0].chains():
        # Check required to avoid overwriting ATOM chain with HETATM one
        if chain.id not in chains:
            chains[chain.id] = chain

    if chain_id is None:
        chain_id = hierarchy.models()[0].chains()[0].id

    return [rg.resseq_as_int() for rg in chains[chain_id].residue_groups()]


def keep_residues(pdbin, pdbout, residue_range, chain_id):
    """Given a range relative to the first residue for a specific chain ID,
    keeps only the residues in that given range

    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    residue_range : list, tuple
       The range of residues to keep
    chain_id : str
       The chain to extract
    
    """
    _, hierarchy, symmetry = _cache(pdbin)

    # Only keep the specified chain ID
    hierarchy = _select(hierarchy, "chain %s" % chain_id)
    
    # Renumber the chain
    _renumber(hierarchy, start=1)

    # Select only residues in the desired range
    residues_to_keep_string = " or ".join(["resseq %d" % r for r in range(residue_range[0], residue_range[1] + 1)])
    hierarchy = _select(hierarchy, residues_to_keep_string)

    # remove hetatms
    _strip(hierarchy, hetatm=True)

    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=['Original file: {0}'.format(pdbin)])


def keep_matching(refpdb=None, targetpdb=None, outpdb=None, resSeqMap=None):
    """Only keep those atoms in targetpdb that are in refpdb and write the result to outpdb.
    We also take care of renaming any chains.
    """

    assert refpdb and targetpdb and outpdb and resSeqMap

    # Paranoid check
    if False:
        refinfo = get_info(refpdb)
        targetinfo = get_info(targetpdb)
        if len(refinfo.models) > 1 or len(targetinfo.models) > 1:
            raise RuntimeError("PDBS contain more than 1 model!")

        if refinfo.models[0].chains != targetinfo.models[0].chains:
            raise RuntimeError("Different numbers/names of chains {0}->{1} between {2} and {3}!".format(
                refinfo.models[0].chains,
                targetinfo.models[0].chains,
                refpdb,
                targetpdb
            ))
            # Now we do our keep matching
    tmp1 = ample_util.tmp_file_name(suffix=".pdb")  # pdbcur insists names have a .pdb suffix

    _keep_matching(refpdb, targetpdb, tmp1, resSeqMap=resSeqMap)

    # now renumber with pdbcur
    logfile = tmp1 + ".log"
    cmd = "pdbcur xyzin {0} xyzout {1}".format(tmp1, outpdb).split()
    stdint = """sernum
"""
    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdint)

    if retcode == 0:
        # remove temporary files
        os.unlink(tmp1)
        os.unlink(logfile)

    return retcode


def _keep_matching(refpdb=None, targetpdb=None, outpdb=None, resSeqMap=None):
    """Create a new pdb file that only contains that atoms in targetpdb that are
    also in refpdb. It only considers ATOM lines and discards HETATM lines in the target.

    Args:
    refpdb: path to pdb that contains the minimal set of atoms we want to keep
    targetpdb: path to the pdb that will be stripped of non-matching atoms
    outpdb: output path for the stripped pdb
    """

    assert refpdb and targetpdb and outpdb and resSeqMap

    def _output_residue(refResidues, targetAtomList, resSeqMap, outfh):
        """Output a single residue only outputting matching atoms, shuffling the atom order and changing the resSeq num"""

        # Get the matching list of atoms
        targetResSeq = targetAtomList[0].resSeq

        refResSeq = resSeqMap.ref2target(targetResSeq)

        # Get the atomlist for the reference
        for (rid, alist) in refResidues:
            if rid == refResSeq:
                refAtomList = alist
                break

        # Get ordered list of the ref atom names for this residue
        rnames = [x.name for x in refAtomList]

        if len(refAtomList) > len(targetAtomList):
            s = "Cannot keep matching as refAtomList is > targetAtomList for residue {0}\nRef: {1}\nTrg: {2}".format(
                targetResSeq,
                rnames,
                [x.name for x in targetAtomList]
            )
            raise RuntimeError(s)

        # Remove any not matching in the target
        alist = []
        for atom in targetAtomList:
            if atom.name in rnames:
                alist.append(atom)

        # List now only contains matching atoms
        targetAtomList = alist

        # Now just have matching so output in the correct order
        for refname in rnames:
            for i, atom in enumerate(targetAtomList):
                if atom.name == refname:
                    # Found the matching atom

                    # Change resSeq and write out
                    atom.resSeq = refResSeq
                    outfh.write(atom.toLine() + "\n")

                    # now delete both this atom and the line
                    targetAtomList.pop(i)

                    # jump out of inner loop
                    break
        return

    # Go through refpdb and find which refResidues are present
    refResidues = []
    targetResSeq = []  # ordered list of tuples - ( resSeq, [ list_of_atoms_for_that_residue ] )

    last = None
    chain = -1
    for line in open(refpdb, 'r'):

        if line.startswith("MODEL"):
            raise RuntimeError("Multi-model file!")

        if line.startswith("TER"):
            break

        if line.startswith("ATOM"):
            a = pdb_model.PdbAtom(line)

            # First atom/chain
            if chain == -1:
                chain = a.chainID

            if a.chainID != chain:
                raise RuntimeError("ENCOUNTERED ANOTHER CHAIN! {0}".format(line))

            if a.resSeq != last:
                last = a.resSeq

                # Add the corresponding resSeq in the target
                targetResSeq.append(resSeqMap.target2ref(a.resSeq))
                refResidues.append((a.resSeq, [a]))
            else:
                refResidues[-1][1].append(a)

    # Now read in target pdb and output everything bar the atoms in this file that
    # don't match those in the refpdb
    t = open(targetpdb, 'r')
    out = open(outpdb, 'w')

    chain = None  # The chain we're reading
    residue = None  # the residue we're reading
    targetAtomList = []

    for line in t:

        if line.startswith("MODEL"):
            raise RuntimeError("Multi-model file!")

        if line.startswith("ANISOU"):
            raise RuntimeError("I cannot cope with ANISOU! {0}".format(line))

        # Stop at TER
        if line.startswith("TER"):
            _output_residue(refResidues, targetAtomList, resSeqMap, out)
            # we write out our own TER
            out.write("TER\n")
            continue

        if line.startswith("ATOM"):

            atom = pdb_model.PdbAtom(line)

            # First atom/chain
            if chain == None:
                chain = atom.chainID

            if atom.chainID != chain:
                raise RuntimeError("ENCOUNTERED ANOTHER CHAIN! {0}".format(line))

            if atom.resSeq in targetResSeq:

                # If this is the first one add the empty tuple and reset residue
                if atom.resSeq != residue:
                    if residue != None:  # Dont' write out owt for first atom
                        _output_residue(refResidues, targetAtomList, resSeqMap, out)
                    targetAtomList = []
                    residue = atom.resSeq

                # If not first keep adding
                targetAtomList.append(atom)

                # We don't write these out as we write them with _output_residue
                continue

            else:
                # discard this line as not a matching atom
                continue

        # For time being exclude all HETATM lines
        elif line.startswith("HETATM"):
            continue
        # Endif line.startswith("ATOM")

        # Output everything else
        out.write(line)

    # End reading loop

    t.close()
    out.close()

    return


def get_info(pdbin):
    """Read a PDB and extract as much information as possible into a PdbInfo object
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    
    Returns
    -------
    :obj:`PdbInfo <ample.util.pdb_model.PdbInfo>`
       A :obj:`PdbInfo <ample.util.pdb_model.PdbInfo>` object
    
    """
    pdb_input = iotbx.pdb.pdb_input(file_name=pdbin)

    # Create a new PdbInfo object
    info = pdb_model.PdbInfo()
    info.pdb = pdbin

    # HEADER & TITLE
    for title in pdb_input.title_section():
        if title.startswith("HEADER"):
            info.pdbCode = title[62:66].strip()
        elif title.startswith('TITLE') and not info.title:
            info.title = title[10:-1].strip()

    # CRYST1
    info.crystalInfo = pdb_model.CrystalInfo()
    info.crystalInfo.z = pdb_input.extract_cryst1_z_columns()
    cryst1 = pdb_input.crystal_symmetry_from_cryst1()
    if cryst1:
        info.crystalInfo.spaceGroup = cryst1.space_group_info()
        info.crystalInfo.unit_cell = cryst1.unit_cell().parameters()

    # REMARK   2
    info.resolution = -1
    for remark in pdb_input.extract_remark_iii_records(2):
        if "RESOLUTION" in remark:
            info.resolution = float(remark.split()[3])

    # REMARK 2809
    info.solventContent = pdb_input.get_solvent_content()
    info.matthewsCoefficient = pdb_input.get_matthews_coeff()

    # MODEL & ATOM
    hierarchy = pdb_input.construct_hierarchy()
    _most_prob(hierarchy, always_keep_one=True)
    _strip(hierarchy, hetatm=True, hydrogen=True)

    for m in hierarchy.models():
        model = pdb_model.PdbModel()
        model.serial = m.id

        for c in m.chains():
            chain_atoms = [pdb_model.PdbAtom(a.format_atom_record()) for a in c.atoms() if not a.hetero]

            if len(chain_atoms) > 0:
                chain_resseq = [rg.resseq_as_int() for rg in c.residue_groups()]
                chain_sequence = "".join([three2one[ag.resname.strip()] for ag in c.atom_groups()])

                model.chains.append(c.id)
                model.atoms.append(chain_atoms)
                model.resSeqs.append(chain_resseq)
                model.sequences.append(chain_sequence)

                model.caMask.append([])
                model.bbMask.append([])

                for ag in c.atom_groups():
                    atoms = [a.name.strip() for a in ag.atoms()]

                    # If    all backbone atoms present
                    # Elif  some backbone atoms missing but C-alpha is there
                    # Else  some backbone atoms and C-alpha missing
                    if all(a in atoms for a in ['N', 'CA', 'C', 'O', 'CB']):
                        foo, bar = False, False
                    elif "CA" in atoms:
                        foo, bar = True, False
                    else:
                        foo, bar = True, True

                    model.bbMask[-1].append(foo)
                    model.caMask[-1].append(bar)

        info.models.append(model)

    return info


def match_resseq(targetPdb=None, outPdb=None, resMap=None, sourcePdb=None):
    """

    """

    assert sourcePdb or resMap
    assert not (sourcePdb and resMap)

    if not resMap:
        resMap = residue_map.residueSequenceMap(targetPdb, sourcePdb)

    target = open(targetPdb, 'r')
    out = open(outPdb, 'w')

    chain = None  # The chain we're reading
    residue = None  # the residue we're reading

    for line in target:

        if line.startswith("MODEL"):
            raise RuntimeError("Multi-model file!")

        if line.startswith("ANISOU"):
            raise RuntimeError("I cannot cope with ANISOU! {0}".format(line))

        # Stop at TER
        if line.startswith("TER"):
            # we write out our own TER
            # out.write("TER\n")
            # break
            pass

        if line.startswith("ATOM"):

            atom = pdb_model.PdbAtom(line)

            # First atom/chain
            if chain == None:
                chain = atom.chainID

            if atom.chainID != chain:
                pass
                # raise RuntimeError, "ENCOUNTERED ANOTHER CHAIN! {0}".format( line )

            # Get the matching resSeq for the model
            modelResSeq = resMap.ref2target(atom.resSeq)
            if modelResSeq == atom.resSeq:
                out.write(line)
            else:
                atom.resSeq = modelResSeq
                out.write(atom.toLine() + "\n")
            continue
        # Endif line.startswith("ATOM")

        # Output everything else
        out.write(line)

    # End reading loop

    target.close()
    out.close()

    return


def merge(pdbin1, pdbin2, pdbout):
    """Merge two pdb files into one
    
    The chains of ``pdbin2`` are appended to ``pdbin1`` corresponding to the models
    
    Parameters
    ----------
    pdbin1 : str
       The path to the first input PDB    
    pdbin2 : str
       The path to the second input PDB
    pdbout : str
       The path to the output PDB
    
    Raises
    ------
    ValueError
       Cannot handle multiple models yet
    
    """
    # Don't save symmetry or anything otherwise header only applicable to one
    _, h1, _ = _cache(pdbin1)
    _, h2, _ = _cache(pdbin2)

    if h1.models_size() > 1 or h2.models_size() > 1:
        msg = "Cannot handle multiple models yet"
        raise ValueError(msg)

    alphabet = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    chains_in_first = {c.id for c in h1.only_model().chains()}
    chain_mapping = {}
    for c, l in zip(h2.only_model().chains(), alphabet[len(chains_in_first):h2.only_model().chains_size()]):
        if all(a.hetero for a in c.atoms()):
            continue
        else:
            chain_mapping[c.id] = l

    m = h1.only_model()
    for c in h2.only_model().chains():
        addable = c.detached_copy()
        addable.id = chain_mapping[addable.id]
        m.append_chain(addable)

    # Reset the atom sequence
    h1.atoms().reset_i_seq()
    _save(pdbout, h1, remarks=['Merge between {0} and {1}'.format(pdbin1, pdbin2),
                               'Chains renamed for {0}'.format(pdbin2)])


def most_prob(pdbin, pdbout, always_keep_one_conformer=True):
    """Remove alternate conforms from the structure file
    
    Depending on the value of always_keep_one_conformer, this will either 
    remove any atom_group with altloc other than blank or 'A', or it will 
    remove any atom_group beyond the first conformer found.

    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    always_keep_one_conformer : bool, optional
       Keep at least a single conformer [default: True]

    """
    _, hierarchy, symmetry = _cache(pdbin)
    _most_prob(hierarchy, always_keep_one=always_keep_one_conformer)
    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=['Original file: {0}'.format(pdbin)])


def molecular_weight(pdbin, first=False):
    """Determine the molecular weight of a pdb
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    first : bool
       Consider the first chain in the first model only [default: False]

    Returns
    -------
    float
       The molecular weight of the extracted atoms

    Notes
    -----
    This function ignores water molecules.
    
    """
    _, hierarchy, _ = _cache(pdbin)
    _, _, mw = _natm_nres_mw(hierarchy, first)
    return mw


def num_atoms_and_residues(pdbin, first=False):
    """Determine the number of atoms and residues in a pdb file.
    
    If all is True, return all atoms and residues, else just for the first chain in the first model

    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    first : bool
       Consider the first chain in the first model only [default: False]

    Returns
    -------
    int
       The number of atoms
    int
       The number of residues

    """
    _, hierarchy, _ = _cache(pdbin)
    natoms, nresidues, _ = _natm_nres_mw(hierarchy, first)
    return natoms, nresidues


def prepare_nmr_model(nmr_model_in, models_dir):
    """Split an nmr pdb into its constituent parts and standardise the lengths
    
    Parameters
    ----------
    nmr_model_in : str
       The path to the input NMR ensemble
    models_dir : str
       The path to the directory for individual NMR model storage
       
    Returns
    -------
    list
       A list of paths to each individual model
    
    """
    if not os.path.isdir(models_dir):
        os.mkdir(models_dir)
    split_pdbs = split_pdb(nmr_model_in, models_dir)

    # We can only work with equally sized PDBS so we pick the most numerous if there are different sizes
    lengths = {}
    lmax = 0
    for pdb in split_pdbs:
        h = iotbx.pdb.pdb_input(pdb).construct_hierarchy()
        l = h.models()[0].chains()[0].residue_groups_size()
        if l not in lengths:
            lengths[l] = [pdb]
        else:
            lengths[l].append(pdb)
        lmax = max(lmax, l)

    if len(lengths) > 1:
        # The pdbs were of different lengths
        to_keep = lengths[lmax]
        logger.info('All NMR models were not of the same length, only %d will be kept.', len(to_keep))
        # Delete any that are not of most numerous length
        for p in [p for p in split_pdbs if p not in to_keep]:
            os.unlink(p)
        split_pdbs = to_keep

    return split_pdbs


def reliable_sidechains(pdbin, pdbout):
    """Only output non-backbone atoms for certain residues

    This function strips side chain atoms of residues not defined in the
    following list:
       ['MET', 'ASP', 'PRO', 'GLN', 'LYS', 'ARG', 'GLU', 'SER']

    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB

    """
    _, hierarchy, symmetry = _cache(pdbin)

    # Remove sidechains that are in res_names where the atom name is not in atom_names
    res_names = ['MET', 'ASP', 'PRO', 'GLN', 'LYS', 'ARG', 'GLU', 'SER']
    atom_names = ['N', 'CA', 'C', 'O', 'CB']
    select_string = "({residues}) or not ({residues}) and ({atoms})".format(
        atoms=" or ".join(['name %s' % atm.lower() for atm in atom_names]),
        residues=" or ".join(['resname %s' % res.upper() for res in res_names]),
    )
    hierarchy = _select(hierarchy, select_string)
    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=['Original file: {0}'.format(pdbin)])


def rename_chains(pdbin, pdbout, fromChain, toChain):
    """Rename chains
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    from_chain : list, tuple
       A list of reference chain names
    to_chain : list, tuple
       A list of target chain names
    
    Raises
    ------
    ValueError
       Renaming lists need to be of equal shape
    
    """
    if len(fromChain) != len(toChain):
        raise ValueError("Renaming lists need to be of equal shape")
    _, hierarchy, symmetry = _cache(pdbin)
    table = dict(zip(fromChain, toChain))
    _rename_chains(hierarchy, table)
    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=['Original file: {0}'.format(pdbin)])


def resseq(pdbin):
    """Get the residue sequence
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    
    Returns
    -------
    dict
       A dictionary of chains and the corresponding sequences
    
    """
    _, hierarchy, _ = _cache(pdbin)
    return _resseq(hierarchy)


def renumber_residues(pdbin, pdbout, start=1):
    """Renumber the residues in a structure file
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    start : int, optional
       The starting number [default: 1]
    
    """
    _, hierarchy, symmetry = _cache(pdbin)
    _renumber(hierarchy, start)
    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=['Original file: {0}'.format(pdbin)])


def renumber_residues_gaps(pdbin, pdbout, gaps, start=1):
    """
    Renumber the residues in the chain based on specified gaps

    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    gaps : list
        List containing True/False for gaps
    start : int, optional
       The starting number [default: 1]

    """
    _, hierarchy, symmetry = _cache(pdbin)
    _renumber_residues_gaps(hierarchy,  gaps, start)
    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=['Original file: {0}'.format(pdbin)])


def rog_side_chain_treatment(pdbin, pdbout, rog_data, del_orange=False):
    """Remove side chains using the ROG scores
        
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    rog_data : list, tuple
       The per-residue ROG data
    del_orange : bool, optional
       Delete residue groups with an orange label [default: False]
       
    """
    resSeq_data = extract_resSeq(pdbin)
    scores = zip(resSeq_data, rog_data)
    _, hierarchy, symmetry = _cache(pdbin)
    _rog_side_chain_treatment(hierarchy, scores, del_orange)
    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=['Original file: {0}'.format(pdbin)])


def select_residues(pdbin, pdbout, delete=None, tokeep=None, delete_idx=None, tokeep_idx=None):
    """Select specified residues in a given PDB structure

    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    delete : list, tuple, optional
       A list of residues to delete
    tokeep : list, tuple, optional
       A list of residues to keep

    """
    _, hierarchy, symmetry = _cache(pdbin)

    if len(hierarchy.models()) > 1 or len(hierarchy.models()[0].chains()) > 1:
        print("pdb {0} has > 1 model or chain - only first model/chain will be kept".format(pdbin))
        _first_chain_only(hierarchy)
    
    chain = hierarchy.models()[0].chains()[0]

    idx = -1
    for residue_group in chain.residue_groups():
        # We ignore hetatms when indexing as we are concerned with residue indexes
        if delete_idx or tokeep_idx:
            if any([atom.hetero for atom in residue_group.atoms()]):
                continue
        idx += 1

        if delete and residue_group.resseq_as_int() not in delete:
            continue
        elif delete_idx and idx not in delete:
            continue
        elif tokeep and residue_group.resseq_as_int() in tokeep: 
            continue
        elif tokeep_idx and idx in tokeep_idx:
            continue
        
        chain.remove_residue_group(residue_group)

    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=['Original file: {0}'.format(pdbin)])


def sequence(pdbin):
    """Get the residue sequence
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    
    Returns
    -------
    dict
       A dictionary of chains and the corresponding sequences

    """
    return _sequence(iotbx.pdb.pdb_input(pdbin).construct_hierarchy())


def sequence_data(pdbin):
    """Get the residue sequence data

    Parameters
    ----------
    pdbin : str
       The path to the input PDB

    Returns
    -------
    dict
       A dictionary of chains and the corresponding sequences

    """
    return _sequence_data(iotbx.pdb.pdb_input(pdbin).construct_hierarchy())


def split_pdb(pdbin, directory=None):
    """Split a pdb file into its separate models
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    directory : str, optional
       A path to a directory to store models in

    Returns
    -------
    list
       The list of split pdb models
    
    Raises
    ------
    RuntimeError
       split_into_chains only works with single-model pdbs

    '"""

    if directory is None: 
        directory = os.path.dirname(pdbin)

    _, hierarchy, symmetry = _cache(pdbin)

    # Nothing to do
    n_models = hierarchy.models_size()
    if n_models == 1:
        raise RuntimeError("split_pdb {0} only contained 1 model!".format(pdbin))

    output_files = []
    for k, model in enumerate(hierarchy.models()):
        k += 1
        hierarchy = iotbx.pdb.hierarchy.root()
        hierarchy.append_model(model.detached_copy())
        if model.id == "":
            model_id = str(k)
        else:
            model_id = model.id.strip()

        output_file = ample_util.filename_append(pdbin, model_id, directory)
        _save(output_file, hierarchy, crystal_symmetry=symmetry,
              remarks=['Model %d of %d' % (k, n_models), 'Original file: %s' % pdbin])
        output_files.append(output_file)

    return output_files


def split_into_chains(pdbin, chain=None, chain_id=None, directory=None):
    """Split a pdb file into its separate chains
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    chain_id : str, optional
       Specify a single chain to extract 
    directory : str, optional
       A path to a directory to store models in

    Returns
    -------
    list
       The list of split pdb models
    
    Raises
    ------
    RuntimeError
       split_into_chains only works with single-model pdbs

    '"""

    if directory is None: 
        directory = os.path.dirname(pdbin)

    if chain:
        warnings.warn("Keyword deprecated - please use chain_id instead")
        chain_id = chain

    _, hierarchy, symmetry = _cache(pdbin)
    
    # Nothing to do
    n_models = hierarchy.models_size()
    if n_models != 1: 
        raise RuntimeError("split_into_chains only works with single-model pdbs!")

    output_files = []
    for i, hchain in enumerate(hierarchy.models()[0].chains()):
        if not hchain.is_protein():
            continue
        if chain_id and chain_id != hchain.id:
            continue
        hierarchy = iotbx.pdb.hierarchy.root()
        model = iotbx.pdb.hierarchy.model()
        hierarchy.append_model(model)
        model.append_chain(hchain.detached_copy())

        output_file = ample_util.filename_append(pdbin, hchain.id, directory)
        _save(output_file, hierarchy, crystal_symmetry=symmetry,
              remarks=['Original file: %s' % pdbin, 'Chain' % hchain.id])

    if not len(output_files): 
        raise RuntimeError("split_into_chains could not find any chains to split")

    return output_files


def standardise(pdbin, pdbout, chain=None, chain_id=None, del_hetatm=False):
    """Standardize a PDB input structure
    
    Standarization includes:
        - renaming of any non-standard AA
        - removal of solvent
        - deletion of less probable rotamer conformations
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    chain_id : str, optional
       The chain to extract
    del_hetatm : bool, optional
       Remove HETATM entries [default: False]
    
    """
    if chain:
        warnings.warn("Keyword deprecated - please use chain_id instead")
        chain_id = chain

    _, hierarchy, symmetry = _cache(pdbin)

    # Remove solvents defined below
    sol_select = " or ".join(
        ["resname {0}".format(sol) for sol in
        {'ADE', 'CYT', 'GUA', 'INO', 'THY', 'URA', 'WAT', 'HOH', 'TIP', 'H2O', 'DOD', 'MOH'}
    ])
    sol_exclude = "not ({0})".format(sol_select)
    hierarchy = _select(hierarchy, sol_exclude)

    # Keep the most probably conformer
    _most_prob(hierarchy, always_keep_one=True)

    # Extract one of the chains
    if chain_id:
        hierarchy = _select(hierarchy, "chain {0}".format(chain))

    tmpf = ample_util.tmp_file_name(delete=False, suffix=".pdb")
    _save(tmpf, hierarchy)

    # Standardise AA names and then remove any remaining HETATMS
    std_residues(tmpf, pdbout, del_hetatm=del_hetatm)
    os.unlink(tmpf)


def std_residues(pdbin, pdbout, del_hetatm=False):
    """Map all residues in MODRES section to their standard counterparts
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    del_hetatm : bool, optional
       Remove HETATM entries [default: False]
    
    """
    pdb_input, hierarchy, symmetry = _cache(pdbin)

    # Get MODRES Section & build up dict mapping the changes
    modres_text = [l.strip() for l in pdb_input.primary_structure_section()
                   if l.startswith("MODRES")]
    modres = {}
    for id, resname, chain, resseq, icode, stdres, comment in _parse_modres(modres_text):
        if not chain in modres:
            modres[chain] = {}
            modres[chain][int(resseq)] = (resname, stdres)

    # TODO: move this to its own function
    for model in hierarchy.models():
        for chain in model.chains():
            for residue_group in chain.residue_groups():
                resseq = residue_group.resseq_as_int()
                for atom_group in residue_group.atom_groups():
                    resname = atom_group.resname
                    if chain.id in modres and resseq in modres[chain.id] and modres[chain.id][resseq][0] == resname:
                        # Change modified name to std name
                        atom_group.resname = modres[chain.id][resseq][1]
                        # If any of the atoms are hetatms, set them to be atoms
                        for atom in atom_group.atoms():
                            if atom.hetero:
                                atom.hetero = False

    if del_hetatm:
        _strip(hierarchy, hetatm=True)

    _save(pdbout, hierarchy, crystal_symmetry=symmetry,
          remarks=['Original file: %s' % pdbin])


def strip(pdbin, pdbout, hetatm=False, hydrogen=False, atom_types=[]):
    """Remove atom types from a structure file
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    hetatm : bool, optional
       Strip the PDB structure of any HETATM entries
    hydrogen : bool, optional
       Strip the PDB structure of all hydrogen elements
    atom_types : list, tuple, optional
       Strip the PDB structure of any specified elements
       
    Raises
    ------
    ValueError
       Define which atoms to strip
    
    """
    if not (hetatm or hydrogen or atom_types):
        msg = "Define which atoms to strip"
        raise ValueError(msg)

    hierarchy, symmetry = _cache(pdbin)
    _strip(hierarchy, hetatm=hetatm, hydrogen=hydrogen, atom_types=atom_types)
    _save(pdbout, hierarchy, crystal_symmetry=symmetry, remarks=['Original file: %s' % pdbin])


def to_single_chain(pdbin, pdbout):
    """Condense a single-model multi-chain pdb to a single-chain pdb
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    
    """
    _, hierarchy, symmetry = _cache(pdbin)
    _first_chain_only(hierarchy)
    _save(pdbout, hierarchy, crystal_symmetry=symmetry)


def translate(pdbin, pdbout, ftranslate):
    """Translate all atoms in a structure file by the provided vector
    
    Parameters
    ----------
    pdbin : str
       The path to the input PDB
    pdbout : str
       The path to the output PDB
    ftranslate : list, tuple
       The vector of fractional coordinates to shift by
    
    """
    _, hierarchy, symmetry = _cache(pdbin)

    # Obtain information about the fractional coordinates
    crystal_info = get_info(pdbin).crystalInfo
    ftranslate = np.asarray([crystal_info.a, crystal_info.b, crystal_info.c]) * np.asarray(ftranslate)
    _translate(hierarchy, ftranslate)
    _save(pdbout, hierarchy, crystal_symmetry=symmetry,
          remarks=['Original file: %s' % pdbin, 'Translated using: [%s]' % ", ".join(map(str, ftranslate))])


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Manipulate PDB files', prefix_chars="-")

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-ren', action='store_true',
                       help="Renumber the PDB")
    group.add_argument('-std', action='store_true',
                       help='Standardise the PDB')
    group.add_argument('-seq', action='store_true',
                       help='Write a fasta of the found AA to stdout')
    group.add_argument('-split_models', action='store_true',
                       help='Split a pdb into constituent models')
    group.add_argument('-split_chains', action='store_true',
                       help='Split a pdb into constituent chains')
    parser.add_argument('input_file',  # nargs='?',
                        help='The input file - will not be altered')
    parser.add_argument('-o', dest='output_file',
                        help='The output file - will be created')
    parser.add_argument('-chain', help='The chain to use')

    args = parser.parse_args()

    # Get full paths to all files
    args.input_file = os.path.abspath(args.input_file)
    if not os.path.isfile(args.input_file):
        raise RuntimeError("Cannot find input file: {0}".format(args.input_file))

    if args.output_file:
        args.output_file = os.path.abspath(args.output_file)
    else:
        n = os.path.splitext(os.path.basename(args.input_file))[0]
        args.output_file = n + "_std.pdb"

    if args.ren:
        renumber_residues(args.input_file, args.output_file, start=1)
    elif args.std:
        standardise(args.input_file, args.output_file, del_hetatm=True, chain=args.chain)
    elif args.seq:
        print(sequence_util.Sequence(pdb=args.input_file).fasta_str())
    elif args.split_models:
        print(split_pdb(args.input_file))
    elif args.split_chains:
        print(split_into_chains(args.input_file, chain=args.chain))

