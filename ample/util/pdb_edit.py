#!/usr/bin/env ccp4-python
'''
Useful manipulations on PDB files
'''

# Python imports
from collections import defaultdict
import copy
import glob
import logging
import os
import re
import sys
import unittest

import iotbx.file_reader
import iotbx.pdb
#iotbx.pdb.amino_acid_codes.one_letter_given_three_letter

import ample_util
import pdb_model
import residue_map
import sequence_util

three2one = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
    'UNK': 'X'
}

# http://stackoverflow.com/questions/3318625/efficient-bidirectional-hash-table-in-python
#aaDict.update( dict((v, k) for (k, v) in aaDict.items()) )
one2three = dict((v, k) for (k, v) in three2one.items())

logger = logging.getLogger(__name__)


def backbone(inpath=None, outpath=None):
    """Only output backbone atoms.
    """

    # pdbcur segfaults with long pathnames
    inpath = os.path.relpath(inpath)
    outpath = os.path.relpath(outpath)

    logfile = outpath + ".log"
    cmd = "pdbcur xyzin {0} xyzout {1}".format(inpath, outpath).split()

    # Build up stdin
    stdin = 'lvatom "N,CA,C,O,CB[N,C,O]"'
    #stdin='lvatom "N,CA,C,O[N,C,O]"'
    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)

    if retcode == 0:
        # remove temporary files
        os.unlink(logfile)
    else:
        raise RuntimeError, "Error stripping PDB to backbone atoms. See log:{0}".format(logfile)

    return


def calpha_only(inpdb, outpdb):
    """Strip PDB to c-alphas only"""

    logfile = outpdb + ".log"
    cmd = "pdbcur xyzin {0} xyzout {1}".format(inpdb, outpdb).split()

    # Build up stdin
    stdin = 'lvatom "CA[C]:*"'
    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)

    if retcode == 0:
        # remove temporary files
        os.unlink(logfile)
    else:
        raise RuntimeError("Error stripping PDB to c-alpha atoms")
    return


def check_pdb_directory(directory, single=True, allsame=True, sequence=None):
    """Check a directory of pdb files to ensure they are valid

    Parameters
    ----------
    direcotry : str
       Path to a directory containing PDB files
    sequence : str
       single-letter protein sequence - if given a check will be made that all
       models are of this sequence
    single : bool
       if True check each pdb only contains a single model
    allsame : bool
       only extract a file if the suffix is in the list

    Returns
    -------
    bool
        True if models all validated
    """
    logger.info("Checking pdbs in directory: {0}".format(directory))
    if not os.path.isdir(directory):
        logger.critical("Cannot find directory: {0}".format(directory))
        return False
    models = glob.glob(os.path.join(directory, "*.pdb"))
    if not len(models):
        logger.critical("Cannot find any pdb files in directory: {0}".format(directory))
        return False
    if not (single or sequence or allsame):
        return True
    return check_pdbs(models, sequence=sequence, single=single, allsame=allsame)


def check_pdbs(models, single=True, allsame=True, sequence=None):
    """Check a list of PDB files to ensure they are valid

    Parameters
    ----------
    models : list
       A list of PDB files
    sequence : str
       single-letter protein sequence - if given a check will be made that all
       models are of this sequence
    single : bool
       if True check each pdb only contains a single model
    allsame : bool
       only extract a file if the suffix is in the list

    Returns
    -------
    bool
        True if models all validated
    """
    if allsame and not sequence:
        # Get sequence from first model
        try:
            h = iotbx.pdb.pdb_input(models[0]).construct_hierarchy()
        except Exception as e:
            s = "*** ERROR reading sequence from first pdb: {0}\n{1}".format(models[0], e)
            logger.critical(s)
            return False
        sequence = _sequence1(h)  # only one model/chain
    errors = []
    multi = []
    no_protein = []
    sequence_err = []
    unnamed_chain = []
    for pdb in models:
        try:
            h = iotbx.pdb.pdb_input(pdb).construct_hierarchy()
        except Exception as e:
            errors.append((pdb, e))
            continue
        if not single:
            continue
        if not (h.models_size() == 1 and h.models()[0].chains_size() == 1):
            multi.append(pdb)
            continue
        # single chain from one model so check is protein
        if not h.models()[0].chains()[0].is_protein():
            no_protein.append(pdb)
            continue
        if not h.models()[0].chains()[0].id.strip():  # Check we have a named chain
            unnamed_chain.append(pdb)
            continue
        if sequence:
            s = _sequence1(h)  # only one chain/model
            if not s == sequence:
                sequence_err.append((pdb, s))

    if not (len(errors) or len(multi) or len(sequence_err) or len(no_protein) or len(unnamed_chain)):
        logger.info("check_pdb_directory - pdb files all seem valid")
        return True

    s = "\n"
    if len(errors):
        s = "*** ERROR ***\n"
        s += "The following pdb files have errors:\n\n"
        for pdb, e in errors:
            s += "{0}: {1}\n".format(pdb, e)

    if len(multi):
        s += "\n"
        s += "The following pdb files have more than one chain:\n\n"
        for pdb in multi:
            s += "{0}\n".format(pdb)

    if len(no_protein):
        s += "\n"
        s += "The following pdb files do not appear to contain any protein:\n\n"
        for pdb in no_protein:
            s += "{0}\n".format(pdb)

    if len(unnamed_chain):
        s += "\n"
        s += "The following pdb files do not have named chains:\n\n"
        for pdb in unnamed_chain:
            s += "{0}\n".format(pdb)

    if len(sequence_err):
        s += "\n"
        s += "The following pdb files have differing sequences from the reference sequence:\n\n{0}\n\n".format(sequence)
        for pdb, seq in sequence_err:
            s += "PDB: {0}\n{1}\n".format(pdb, seq)

    logger.critical(s)
    return False


def extract_chain(inpdb, outpdb, chainID=None, newChainID=None, cAlphaOnly=False, renumber=True):
    """Extract chainID from inpdb and renumner.
    If cAlphaOnly is set, strip down to c-alpha atoms
    """
    logfile = outpdb + ".log"
    cmd = "pdbcur xyzin {0} xyzout {1}".format(inpdb, outpdb).split()

    # Build up stdin
    stdin = "lvchain {0}\n".format(chainID)
    if newChainID:
        stdin += "renchain {0} {1}\n".format(chainID, newChainID)
    if cAlphaOnly:
        stdin += 'lvatom "CA[C]:*"\n'
    if renumber:
        stdin += "sernum\n"

    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)

    if retcode == 0:
        # remove temporary files
        os.unlink(logfile)
    else:
        raise RuntimeError, "Error extracting chain {0}".format(chainID)

    return


def extract_model(inpdb, outpdb, modelID):
    """Extract modelID from inpdb into outpdb"""

    assert modelID > 0

    logfile = outpdb + ".log"
    cmd = "pdbcur xyzin {0} xyzout {1}".format(inpdb, outpdb).split()

    # Build up stdin
    stdin = "lvmodel /{0}\n".format(modelID)
    #stdin += "sernum\n"

    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)

    if retcode != 0:
        raise RuntimeError, "Problem extracting model with cmd: {0}".format

    # remove temporary files
    os.unlink(logfile)

    return


def extract_header_pdb_code(pdb_input):
    for line in pdb_input.title_section():
        if line.startswith("HEADER ") and len(line) >= 65: return line[62:66]
    return None


def extract_header_title(pdb_input):
    for line in pdb_input.title_section():
        if line.startswith('TITLE'): return line[10:-1].strip()
    return None


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
            raise RuntimeError, "PDBS contain more than 1 model!"

        if refinfo.models[0].chains != targetinfo.models[0].chains:
            raise RuntimeError, "Different numbers/names of chains {0}->{1} between {2} and {3}!".format(
                refinfo.models[0].chains, targetinfo.models[0].chains, refpdb, targetpdb)
        # Now we do our keep matching
    tmp1 = ample_util.tmp_file_name() + ".pdb"  # pdbcur insists names have a .pdb suffix

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
                targetResSeq, rnames, [x.name for x in targetAtomList])
            raise RuntimeError, s

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
            raise RuntimeError, "Multi-model file!"

        if line.startswith("TER"):
            break

        if line.startswith("ATOM"):
            a = pdb_model.PdbAtom(line)

            # First atom/chain
            if chain == -1:
                chain = a.chainID

            if a.chainID != chain:
                raise RuntimeError, "ENCOUNTERED ANOTHER CHAIN! {0}".format(line)

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
            raise RuntimeError, "Multi-model file!"

        if line.startswith("ANISOU"):
            raise RuntimeError, "I cannot cope with ANISOU! {0}".format(line)

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
                raise RuntimeError, "ENCOUNTERED ANOTHER CHAIN! {0}".format(line)

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
        #Endif line.startswith("ATOM")

        # Output everything else
        out.write(line)

    # End reading loop

    t.close()
    out.close()

    return


def get_info(inpath):
    """Read a PDB and extract as much information as possible into a PdbInfo object
    """

    info = pdb_model.PdbInfo()
    info.pdb = inpath

    currentModel = None
    currentChain = -1

    modelAtoms = []  # list of models, each of which is a list of chains with the list of atoms

    # Go through refpdb and find which ref_residues are present
    f = open(inpath, 'r')
    line = f.readline()
    while line:

        # First line of title
        if line.startswith('HEADER'):
            info.pdbCode = line[62:66].strip()

        # First line of title
        if line.startswith('TITLE') and not info.title:
            info.title = line[10:-1].strip()

        if line.startswith("REMARK"):

            try:
                numRemark = int(line[7:10])
            except ValueError:
                line = f.readline()
                continue

            # Resolution
            if numRemark == 2:
                line = f.readline()
                if line.find("RESOLUTION") != -1:
                    try:
                        info.resolution = float(line[25:30])
                    except ValueError:
                        # RESOLUTION. NOT APPLICABLE.
                        info.resolution = -1

            # Get solvent content
            if numRemark == 280:

                maxread = 5
                # Clunky - read up to maxread lines to see if we can get the information we're after
                # We assume the floats are at the end of the lines
                for _ in range(maxread):
                    line = f.readline()
                    if line.find("SOLVENT CONTENT") != -1:
                        try:
                            info.solventContent = float(line.split()[-1])
                        except ValueError:
                            # Leave as None
                            pass
                    if line.find("MATTHEWS COEFFICIENT") != -1:
                        try:
                            info.matthewsCoefficient = float(line.split()[-1])
                        except ValueError:
                            # Leave as None
                            pass
        #End REMARK

        if line.startswith("CRYST1"):
            try:
                info.crystalInfo = pdb_model.CrystalInfo(line)
            except ValueError as e:
                logger.critical("ERROR READING CRYST1 LINE in file %s\":%s\"\n%s", inpath, line.rstrip(), e)
                info.crystalInfo = None

        if line.startswith("MODEL"):
            if currentModel:
                # Need to make sure that we have an id if only 1 chain and none given
                if len(currentModel.chains) <= 1:
                    if currentModel.chains[0] == None:
                        currentModel.chains[0] = 'A'

                info.models.append(currentModel)

            # New/first model
            currentModel = pdb_model.PdbModel()
            # Get serial
            currentModel.serial = int(line.split()[1])

            currentChain = None
            modelAtoms.append([])

        # Count chains (could also check against the COMPND line if present?)
        if line.startswith('ATOM'):

            # Create atom object
            atom = pdb_model.PdbAtom(line)

            # Check for the first model
            if not currentModel:
                # This must be the first model and there should only be one
                currentModel = pdb_model.PdbModel()
                modelAtoms.append([])

            if atom.chainID != currentChain:
                currentChain = atom.chainID
                currentModel.chains.append(currentChain)
                modelAtoms[-1].append([])

            modelAtoms[-1][-1].append(atom)

        # Can ignore TER and ENDMDL for time being as we'll pick up changing chains anyway,
        # and new models get picked up by the models line

        line = f.readline()
        # End while loop

    # End of reading loop so add the last model to the list
    info.models.append(currentModel)

    f.close()

    bbatoms = ['N', 'CA', 'C', 'O', 'CB']

    # Now process the atoms
    for modelIdx, model in enumerate(info.models):

        chainList = modelAtoms[modelIdx]

        for chainIdx, atomList in enumerate(chainList):

            # Paranoid check
            assert model.chains[chainIdx] == atomList[0].chainID

            # Add list of atoms to model
            model.atoms.append(atomList)

            # Initialise new chain
            currentResSeq = atomList[0].resSeq
            currentResName = atomList[0].resName
            model.resSeqs.append([])
            model.sequences.append("")
            model.caMask.append([])
            model.bbMask.append([])

            atomTypes = []
            for i, atom in enumerate(atomList):

                aname = atom.name.strip()
                if atom.resSeq != currentResSeq and i == len(atomList) - 1:
                    # Edge case - last residue containing one atom
                    atomTypes = [aname]
                else:
                    if aname not in atomTypes:
                        atomTypes.append(aname)

                if atom.resSeq != currentResSeq or i == len(atomList) - 1:
                    # End of reading the atoms for a residue
                    model.resSeqs[chainIdx].append(currentResSeq)
                    model.sequences[chainIdx] += three2one[currentResName]

                    if 'CA' not in atomTypes:
                        model.caMask[chainIdx].append(True)
                    else:
                        model.caMask[chainIdx].append(False)

                    missing = False
                    for bb in bbatoms:
                        if bb not in atomTypes:
                            missing = True
                            break

                    if missing:
                        model.bbMask[chainIdx].append(True)
                    else:
                        model.bbMask[chainIdx].append(False)

                    currentResSeq = atom.resSeq
                    currentResName = atom.resName
                    atomTypes = []

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
            raise RuntimeError, "Multi-model file!"

        if line.startswith("ANISOU"):
            raise RuntimeError, "I cannot cope with ANISOU! {0}".format(line)

        # Stop at TER
        if line.startswith("TER"):
            # we write out our own TER
            #out.write("TER\n")
            #break
            pass

        if line.startswith("ATOM"):

            atom = pdb_model.PdbAtom(line)

            # First atom/chain
            if chain == None:
                chain = atom.chainID

            if atom.chainID != chain:
                pass
                #raise RuntimeError, "ENCOUNTERED ANOTHER CHAIN! {0}".format( line )

            # Get the matching resSeq for the model
            modelResSeq = resMap.ref2target(atom.resSeq)
            if modelResSeq == atom.resSeq:
                out.write(line)
            else:
                atom.resSeq = modelResSeq
                out.write(atom.toLine() + "\n")
            continue
        #Endif line.startswith("ATOM")

        # Output everything else
        out.write(line)

    # End reading loop

    target.close()
    out.close()

    return


#     def match_resseq(self, targetPdb, sourcePdb, keepAtoms="all", workdir=None, resSeqMap=None ):
#         """Given a native pdb file and a model pdb file, create a copy of the native that can be directly compared with the model
#
#         args:
#         targetPdb:
#         sourcePdb:
#         keepAtoms: all, backbone or calpha - the atoms which are to be kept for the comparision
#
#         """
#
#         if not workdir:
#             workdir = os.curdir()
#
#         if not resSeqMap:
#             # Calculate the RefSeqMap - need to do this before we reduce to c-alphas
#             resSeqMap = residue_map.residueSequenceMap( targetPdb, sourcePdb )
#
#         # Find out if there are atoms in the model that we need to remove
#         modelIncomparable = resSeqMap.modelIncomparable()
#         if len( modelIncomparable ):
#
#             n = os.path.splitext( os.path.basename( targetPdb ) )[0]
#             nativePdbCut = os.path.join( workdir, n+"_cut.pdb" )
#
#             logfile = "{0}.log".format( nativePdbCut )
#             cmd="pdbcur xyzin {0} xyzout {1}".format( targetPdb, nativePdbCut ).split()
#
#             # Build up stdin - I'm too thick to work out the selection syntax for a discrete list
#             stdin = ""
#             for e in modelIncomparable:
#                 stdin += "delresidue {0}\n".format( e )
#
#             retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=workdir, dolog=False, stdin=stdin)
#
#             if retcode == 0:
#                 # remove temporary files
#                 os.unlink(logfile)
#             else:
#                 raise RuntimeError,"Error deleting residues {0}".format( modelIncomparable )
#
#             targetPdb = nativePdbCut
#
#
#         if keepAtoms == "calpha":
#             # If only alpha atoms are required, we create a copy of the model with only alpha atoms
#             n = os.path.splitext( os.path.basename( targetPdb ) )[0]
#             tmp = os.path.join( workdir, n+"_cAlphaOnly.pdb" )
#             calpha_only( targetPdb, tmp )
#             targetPdb = tmp
#         elif keepAtoms == "backbone":
#             # Strip down to backbone atoms
#             n = os.path.splitext( os.path.basename( targetPdb ) )[0]
#             tmp = os.path.join( workdir, n+"_backbone.pdb" )
#             backbone( targetPdb, tmp  )
#             targetPdb = tmp
#         elif keepAtoms == "all":
#             pass
#         else:
#             raise RuntimeError,"Unrecognised keepAtoms: {0}".format( keepAtoms )
#
#         # Now create a PDB with the matching atoms from native that are in refined
#         n = os.path.splitext( os.path.basename( targetPdb ) )[0]
#         nativePdbMatch = os.path.join( workdir, n+"_matched.pdb" )
#         keep_matching( refpdb=refinedPdb, targetpdb=targetPdb, outpdb=nativePdbMatch, resSeqMap=resSeqMap )
#
#         return


def merge(pdb1=None, pdb2=None, pdbout=None):
    """Merge two pdb files into one"""

    logfile = pdbout + ".log"
    cmd = ['pdb_merge', 'xyzin1', pdb1, 'xyzin2', pdb2, 'xyzout', pdbout]

    # Build up stdin
    stdin = 'nomerge'
    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)

    if retcode == 0:
        # remove temporary files
        os.unlink(logfile)
    else:
        raise RuntimeError, "Error merging pdbs: {0} {1}".format(pdb1, pdb2)

    return


def molecular_weight(pdbin):
    logfile = "rwcontents.log"
    _run_rwcontents(pdbin, logfile)
    _, _, mw = _parse_rwcontents(logfile)
    os.unlink(logfile)
    return mw


def num_atoms_and_residues(pdbin, first=False):
    """"Return number of atoms and residues in a pdb file.
    If all is True, return all atoms and residues, else just for the first chain in the first model'
    """
    #pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdbin)
    #model = pdb_obj.hierarchy.models()[0]
    #return sum(  [ len( chain.residues() ) for chain in model.chains() ]  )

    if not first:
        logfile = "rwcontents.log"
        _run_rwcontents(pdbin, logfile)
        natoms, nresidues, _ = _parse_rwcontents(logfile)
        os.unlink(logfile)
    else:
        pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdbin)
        model = pdb_obj.hierarchy.models()[0]
        nresidues = len(model.chains()[0].residues())
        natoms = len(model.chains()[0].atoms())

    assert natoms > 0 and nresidues > 0

    return (natoms, nresidues)


def _only_equal_sizes(hierarchy):
    """If a hiearchy contains different size models, only keep models of the most numerous size"""
    lengths = defaultdict(list)
    lmax = 0
    for i, model in enumerate(hierarchy.models()):
        l = model.chains()[0].residue_groups_size()
        lengths[l].append(i)
        lmax = max(lmax, l)

    if len(lengths) > 1:
        # The pdbs were of different lengths
        to_keep = lengths[lmax]
        logger.debug('All models were not of the same length, only {0} will be kept.'.format(len(to_keep)))
        # Delete any that are not of most numerous length
        for i, model in enumerate(hierarchy.models()):
            if i not in to_keep: hierarchy.remove_model(model)
    return hierarchy


def _parse_rwcontents(logfile):
    natoms = 0
    nresidues = 0
    molecular_weight = 0
    with open(logfile) as f:
        for line in f:
            if line.startswith(" Number of amino-acids residues"):
                nresidues = int(line.strip().split()[5])
            #Total number of protein atoms (including hydrogens)
            if line.startswith(" Total number of         atoms (including hydrogens)"):
                natoms = int(float(line.strip().split()[6]))
            if line.startswith(" Molecular Weight of protein:"):
                molecular_weight = float(line.strip().split()[4])
    return natoms, nresidues, molecular_weight


def _run_rwcontents(pdbin, logfile):
    logfile = os.path.abspath(logfile)
    cmd = ['rwcontents', 'xyzin', pdbin]
    stdin = ''  # blank to trigger EOF
    retcode = ample_util.run_command(cmd=cmd, directory=os.getcwd(), logfile=logfile, stdin=stdin)
    if retcode != 0:
        raise RuntimeError, "Error running cmd {0}\nSee logfile: {1}".format(cmd, logfile)
    return


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

        idCode = line[7:11]
        resName = line[12:15].strip()
        # Use for all so None means an empty field
        if line[16].strip(): chainID = line[16]
        seqNum = int(line[18:22])
        iCode = ""
        if line[22].strip(): iCode = line[22]
        stdRes = line[24:27].strip()
        comment = ""
        if line[29:70].strip(): comment = line[29:70].strip()

        modres.append([idCode, resName, chainID, seqNum, iCode, stdRes, comment])

    return modres


def reliable_sidechains(inpath=None, outpath=None):
    """Only output non-backbone atoms for residues in the res_names list.
    """

    # Remove sidechains that are in res_names where the atom name is not in atom_names
    res_names = ['MET', 'ASP', 'PRO', 'GLN', 'LYS', 'ARG', 'GLU', 'SER']
    atom_names = ['N', 'CA', 'C', 'O', 'CB']

    pdb_in = open(inpath, "r")
    pdb_out = open(outpath, "w")

    for pdbline in pdb_in:
        pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
        pdb_result = pdb_pattern.match(pdbline)

        # Check ATOM line and for residues in res_name, skip any that are not in atom names
        if pdb_result:
            pdb_result2 = re.split(pdb_pattern, pdbline)
            if pdb_result2[3] in res_names and not pdb_result2[2] in atom_names:
                continue

        # Write out everything else
        pdb_out.write(pdbline)

    #End for
    pdb_out.close()
    pdb_in.close()

    return


def reliable_sidechains_cctbx(pdbin=None, pdbout=None):
    """Only output non-backbone atoms for residues in the res_names list.
    """

    # Remove sidechains that are in res_names where the atom name is not in atom_names
    res_names = ['MET', 'ASP', 'PRO', 'GLN', 'LYS', 'ARG', 'GLU', 'SER']
    atom_names = ['N', 'CA', 'C', 'O', 'CB']

    pdb_input = iotbx.pdb.pdb_input(pdbin)
    hierachy = pdb_input.construct_hierarchy()
    # Remove HETATMS
    for model in hierachy.models():
        for chain in model.chains():
            for residue_group in chain.residue_groups():
                assert not residue_group.have_conformers(), "Fix for conformers"
                if residue_group.unique_resnames()[0] not in res_names:
                    # removing whilst looping through?!? - maybe...
                    chain.remove_residue_group(residue_group)
                    continue
                for atom_group in residue_group.atom_groups():
                    # Can't use below as it uses indexes which change as we remove atoms
                    # ag.atoms().extract_hetero()]
                    todel = [a for a in atom_group.atoms() if a.name.strip() in atom_names]
                    for a in todel:
                        atom_group.remove_atom(a)

    # Need to get crystal info and include
    hierachy.write_pdb_file(pdbout, anisou=False)
    return


def rename_chains(inpdb=None, outpdb=None, fromChain=None, toChain=None):
    """Rename Chains
    """

    assert len(fromChain) == len(toChain)

    logfile = outpdb + ".log"
    cmd = "pdbcur xyzin {0} xyzout {1}".format(inpdb, outpdb).split()

    # Build up stdin
    stdin = ""
    for i in range(len(fromChain)):
        stdin += "renchain {0} {1}\n".format(fromChain[i], toChain[i])

    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)

    if retcode == 0:
        # remove temporary files
        os.unlink(logfile)
    else:
        raise RuntimeError, "Error renaming chains {0}".format(fromChain)

    return


def resseq(pdbin):
    return _resseq(iotbx.pdb.pdb_input(pdbin).construct_hierarchy())


def _resseq(hierarchy):
    """Extract the sequence of residues from a pdb file."""
    chain2data = _sequence_data(hierarchy)
    return dict((k, chain2data[k][1]) for k in chain2data.keys())


def renumber_residues(pdbin, pdbout, start=1):
    """ Renumber the residues in the chain """
    pdb_input = iotbx.pdb.pdb_input(file_name=pdbin)
    hierarchy = pdb_input.construct_hierarchy()

    _renumber(hierarchy, start)

    with open(pdbout, 'w') as f:
        f.write("REMARK Original file:\n")
        f.write("REMARK   {0}\n".format(pdbin))
        f.write(hierarchy.as_pdb_string(anisou=False))
    return


def _renumber(hierarchy, start):
    for model in hierarchy.models():
        for chain in model.chains():
            for idx, residue_group in enumerate(chain.residue_groups()):
                residue_group.resseq = idx + start
    return


def renumber_residues_gaps(pdbin, pdbout, gaps, start=1):
    """
    Renumber the residues in the chain based on specified gaps

    Parameters
    ----------
    pdbin : str
    pdbout : str
    gaps : list
        List containing True/False for gaps
    """
    pdb_input = iotbx.pdb.pdb_input(file_name=pdbin)
    hierarchy = pdb_input.construct_hierarchy()

    for model in hierarchy.models():
        for chain in model.chains():
            resseq = 0
            for idx, is_gap in enumerate(gaps):
                if is_gap:
                    continue
                try:
                    residue_group = chain.residue_groups()[resseq]
                except:
                    pass
                else:
                    residue_group.resseq = idx + start
                finally:
                    resseq += 1

    with open(pdbout, 'w') as f:
        f.write("REMARK Original file:\n")
        f.write("REMARK   {0}\n".format(pdbin))
        f.write(hierarchy.as_pdb_string(anisou=False))
    return


def select_residues(pdbin, pdbout, delete=None, tokeep=None, delete_idx=None, tokeep_idx=None):

    pdbf = iotbx.file_reader.any_file(pdbin, force_type="pdb")
    pdbf.check_file_type("pdb")
    hierarchy = pdbf.file_object.construct_hierarchy()

    crystal_symmetry = pdbf.file_object.crystal_symmetry()

    if len(hierarchy.models()) > 1 or len(hierarchy.models()[0].chains()) > 1:
        logger.debug("pdb %s has > 1 model or chain - only first model/chain will be kept", pdbin)

    if len(hierarchy.models()) > 1:
        for i, m in enumerate(hierarchy.models()):
            if i != 0: hierarchy.remove_model(m)
    model = hierarchy.models()[0]

    if len(model.chains()) > 1:
        for i, c in enumerate(model.chains()):
            if i != 0: model.remove_chain(c)
    chain = model.chains()[0]

    idx = -1
    for residue_group in chain.residue_groups():
        # We ignore hetatms when indexing as we are concerned with residue indexes
        if delete_idx or tokeep_idx:
            if any([atom.hetero for atom in residue_group.atoms()]): continue
        idx += 1

        remove = False
        if delete:
            if residue_group.resseq_as_int() in delete: remove = True
        elif delete_idx:
            if idx in delete: remove = True
        elif tokeep:
            if residue_group.resseq_as_int() not in tokeep: remove = True
        elif tokeep_idx:
            if idx not in tokeep_idx: remove = True

        if remove:
            chain.remove_residue_group(residue_group)

    #hierarchy.write_pdb_file(pdbout,anisou=False)
    with open(pdbout, 'w') as f:
        f.write("REMARK Original file:\n")
        f.write("REMARK   {0}\n".format(pdbin))
        if (crystal_symmetry is not None):
            f.write(
                iotbx.pdb.format_cryst1_and_scale_records(crystal_symmetry=crystal_symmetry, write_scale_records=True) +
                "\n")
        f.write(hierarchy.as_pdb_string(anisou=False))
    return


def sequence(pdbin):
    return _sequence(iotbx.pdb.pdb_input(pdbin).construct_hierarchy())


def _sequence(hierarchy):
    """Extract the sequence of residues from a pdb file."""
    chain2data = _sequence_data(hierarchy)
    return dict((k, chain2data[k][0]) for k in chain2data.keys())


def _sequence1(hierarchy):
    """Return sequence of the first chain"""
    d = _sequence(hierarchy)
    return d[sorted(d.keys())[0]]


def sequence_data(pdbin):
    return _sequence_data(iotbx.pdb.pdb_input(pdbin).construct_hierarchy())


def _sequence_data(hierarchy):
    """Extract the sequence of residues and resseqs from a pdb file."""
    chain2data = {}
    for chain in set(hierarchy.models()[0].chains()):  # only the first model
        if not chain.is_protein(): continue
        got = False
        seq = ""
        resseq = []
        for residue in chain.conformers()[0].residues():  # Just look at the first conformer
            # See if any of the atoms are non-hetero - if so we add this residue
            if any([not atom.hetero for atom in residue.atoms()]):
                got = True
                seq += three2one[residue.resname]
                #resseq.append(int(residue.resseq.strip()))
                resseq.append(residue.resseq_as_int())
        if got: chain2data[chain.id] = (seq, resseq)
    return chain2data


def split_pdb(pdbin, directory=None, strip_hetatm=False, same_size=False):
    """Split a pdb file into its separate models

    Parameters
    ----------
    pdbin : str
      path to input pdbf file
    directory : str
      path to directory where pdb files will be created
    strip_hetatm : bool
        remove HETATMS if true
    same_size : bool
      Only output models of equal length (the most numerous length is selected)
    """

    if directory is None: directory = os.path.dirname(pdbin)
    if not os.path.isdir(directory): os.mkdir(directory)

    # Largely stolen from pdb_split_models.py in phenix
    #http://cci.lbl.gov/cctbx_sources/iotbx/command_line/pdb_split_models.py

    pdbf = iotbx.file_reader.any_file(pdbin, force_type="pdb")
    pdbf.check_file_type("pdb")
    hierarchy = pdbf.file_object.construct_hierarchy()

    # Nothing to do
    n_models = hierarchy.models_size()
    if n_models == 1: raise RuntimeError, "split_pdb {0} only contained 1 model!".format(pdbin)

    if same_size: _only_equal_sizes(hierarchy)
    crystal_symmetry = pdbf.file_object.crystal_symmetry()
    output_files = []
    for k, model in enumerate(hierarchy.models()):
        k += 1
        new_hierarchy = iotbx.pdb.hierarchy.root()
        new_hierarchy.append_model(model.detached_copy())
        if strip_hetatm: _strip(
                new_hierarchy,
                hetatm=True,
        )
        if (model.id == ""):
            model_id = str(k)
        else:
            model_id = model.id.strip()
        output_file = ample_util.filename_append(pdbin, model_id, directory)
        with open(output_file, "w") as f:
            if (crystal_symmetry is not None):
                print >> f, iotbx.pdb.format_cryst1_and_scale_records(
                    crystal_symmetry=crystal_symmetry, write_scale_records=True)
            print >> f, "REMARK Model %d of %d" % (k, n_models)
            if (pdbin is not None):
                print >> f, "REMARK Original file:"
                print >> f, "REMARK   %s" % pdbin
            f.write(new_hierarchy.as_pdb_string())
        output_files.append(output_file)
    return output_files


def split_into_chains(pdbin, chain=None, directory=None):
    """Split a pdb file into its separate chains"""

    if directory is None: directory = os.path.dirname(pdbin)

    # Largely stolen from pdb_split_models.py in phenix
    #http://cci.lbl.gov/cctbx_sources/iotbx/command_line/pdb_split_models.py
    pdbf = iotbx.file_reader.any_file(pdbin, force_type="pdb")
    pdbf.check_file_type("pdb")
    hierarchy = pdbf.file_object.construct_hierarchy()

    # Nothing to do
    n_models = hierarchy.models_size()
    if n_models != 1: raise RuntimeError, "split_into_chains only works with single-mdoel pdbs!"

    crystal_symmetry = pdbf.file_object.crystal_symmetry()

    output_files = []
    n_chains = len(hierarchy.models()[0].chains())
    for i, hchain in enumerate(hierarchy.models()[0].chains()):
        if not hchain.is_protein(): continue
        if chain and not hchain.id == chain: continue
        new_hierarchy = iotbx.pdb.hierarchy.root()
        new_model = iotbx.pdb.hierarchy.model()
        new_hierarchy.append_model((new_model))
        new_model.append_chain(hchain.detached_copy())
        output_file = ample_util.filename_append(pdbin, hchain.id, directory)
        with open(output_file, "w") as f:
            if (crystal_symmetry is not None):
                print >> f, iotbx.pdb.format_cryst1_and_scale_records(
                    crystal_symmetry=crystal_symmetry, write_scale_records=True)
            print >> f, "REMARK Chain %d of %d" % (i, n_chains)
            if (pdbin is not None):
                print >> f, "REMARK Original file:"
                print >> f, "REMARK   %s" % pdbin
            f.write(new_hierarchy.as_pdb_string())

        output_files.append(output_file)

    if not len(output_files): raise RuntimeError, "split_into_chains could not find any chains to split"

    return output_files


def standardise(pdbin, pdbout, chain=None, del_hetatm=False):
    """Rename any non-standard AA, remove solvent and only keep most probably conformation.
    """

    tmp1 = ample_util.tmp_file_name() + ".pdb"  # pdbcur insists names have a .pdb suffix

    # Now clean up with pdbcur
    logfile = tmp1 + ".log"
    cmd = "pdbcur xyzin {0} xyzout {1}".format(pdbin, tmp1).split()
    stdin = """delsolvent
noanisou
mostprob
"""
    # We are extracting one  of the chains
    if chain: stdin += "lvchain {0}\n".format(chain)

    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)
    if retcode == 0: os.unlink(logfile)  # remove temporary files
    else: raise RuntimeError, "Error standardising pdb!"

    # Standardise AA names and then remove any remaining HETATMs
    std_residues_cctbx(tmp1, pdbout, del_hetatm=del_hetatm)
    os.unlink(tmp1)

    return retcode


def std_residues_cctbx(pdbin, pdbout, del_hetatm=False):
    """Map all residues in MODRES section to their standard counterparts
    optionally delete all other HETATMS"""

    pdb_input = iotbx.pdb.pdb_input(pdbin)
    crystal_symmetry = pdb_input.crystal_symmetry()

    # Get MODRES Section & build up dict mapping the changes
    modres_text = [ l.strip() for l in pdb_input.primary_structure_section() \
                    if l.startswith("MODRES")]
    modres = {}
    for id, resname, chain, resseq, icode, stdres, comment in _parse_modres(modres_text):
        if not chain in modres:
            modres[chain] = {}
        modres[chain][int(resseq)] = (resname, stdres)

    hierachy = pdb_input.construct_hierarchy()
    for model in hierachy.models():
        for chain in model.chains():
            for residue_group in chain.residue_groups():
                resseq = residue_group.resseq_as_int()
                for atom_group in residue_group.atom_groups():
                    resname = atom_group.resname
                    if chain.id in modres and resseq in modres[chain.id] and modres[chain.id][resseq][0] == resname:
                        # Change modified name to std name
                        #assert modres[chain.id][resseq][0]==resname,\
                        #"Unmatched names: {0} : {1}".format(modres[chain.id][resseq][0],resname)
                        atom_group.resname = modres[chain.id][resseq][1]
                        # If any of the atoms are hetatms, set them to be atoms
                        for atom in atom_group.atoms():
                            if atom.hetero: atom.hetero = False

    if del_hetatm: _strip(hierachy, hetatm=True)

    with open(pdbout, 'w') as f:
        f.write("REMARK Original file:\n")
        f.write("REMARK   {0}\n".format(pdbin))
        if (crystal_symmetry is not None):
            f.write(
                iotbx.pdb.format_cryst1_and_scale_records(crystal_symmetry=crystal_symmetry, write_scale_records=True) +
                "\n")
        f.write(hierachy.as_pdb_string(anisou=False))
    return


def strip(pdbin, pdbout, hetatm=False, hydrogen=False, atom_types=[]):
    assert hetatm or hydrogen or atom_types, "Need to set what to strip!"

    pdb_input = iotbx.pdb.pdb_input(pdbin)
    crystal_symmetry = pdb_input.crystal_symmetry()

    hierachy = pdb_input.construct_hierarchy()
    _strip(hierachy, hetatm=hetatm, hydrogen=hydrogen, atom_types=atom_types)

    with open(pdbout, 'w') as f:
        f.write("REMARK Original file:\n")
        f.write("REMARK   {0}\n".format(pdbin))
        if (crystal_symmetry is not None):
            f.write(
                iotbx.pdb.format_cryst1_and_scale_records(crystal_symmetry=crystal_symmetry, write_scale_records=True) +
                "\n")
        f.write(hierachy.as_pdb_string(anisou=False))
    return


def _strip(hierachy, hetatm=False, hydrogen=False, atom_types=[]):
    """Remove all hetatoms from pdbfile"""

    def remove_atom(atom, hetatm=False, hydrogen=False, atom_types=[]):
        return (hetatm and atom.hetero) or (hydrogen and atom.element_is_hydrogen()) or atom.name.strip() in atom_types

    for model in hierachy.models():
        for chain in model.chains():
            for residue_group in chain.residue_groups():
                for atom_group in residue_group.atom_groups():
                    to_del = [
                        a for a in atom_group.atoms()
                        if remove_atom(a, hetatm=hetatm, hydrogen=hydrogen, atom_types=atom_types)
                    ]
                    for atom in to_del:
                        atom_group.remove_atom(atom)
    return


def to_single_chain(inpath, outpath):
    """Condense a single-model multi-chain pdb to a single-chain pdb"""

    o = open(outpath, 'w')

    firstChainID = None
    currentResSeq = None  # current residue we are reading - assume it always starts from 1
    globalResSeq = None
    globalSerial = -1
    for line in open(inpath):

        # Remove any HETATOM lines and following ANISOU lines
        if line.startswith("HETATM") or line.startswith("MODEL") or line.startswith("ANISOU"):
            raise RuntimeError, "Cant cope with the line: {0}".format(line)

        # Skip any TER lines
        if line.startswith("TER"):
            continue

        if line.startswith("ATOM"):
            changed = False

            atom = pdb_model.PdbAtom(line)

            # First atom/residue
            if globalSerial == -1:
                globalSerial = atom.serial
                firstChainID = atom.chainID
                globalResSeq = atom.resSeq
                currentResSeq = atom.resSeq
            else:
                # Change residue numbering and chainID
                if atom.chainID != firstChainID:
                    atom.chainID = firstChainID
                    changed = True

                # Catch each change in residue
                if atom.resSeq != currentResSeq:
                    # Change of residue
                    currentResSeq = atom.resSeq
                    globalResSeq += 1

                # Only change if don't match global
                if atom.resSeq != globalResSeq:
                    atom.resSeq = globalResSeq
                    changed = True

                # Catch each change in numbering
                if atom.serial != globalSerial + 1:
                    atom.serial = globalSerial + 1
                    changed = True

                if changed:
                    line = atom.toLine() + "\n"

                # Increment counter for all atoms
                globalSerial += 1

        o.write(line)

    o.close()

    return


def translate(inpdb=None, outpdb=None, ftranslate=None):
    """translate pdb
    args:
    ftranslate -- vector of fractional coordinates to shift by
    """

    logfile = outpdb + ".log"
    cmd = "pdbcur xyzin {0} xyzout {1}".format(inpdb, outpdb).split()

    # Build up stdin
    stdin = 'translate * frac {0:F} {1:F} {2:F}'.format(ftranslate[0], ftranslate[1], ftranslate[2])
    retcode = ample_util.run_command(cmd=cmd, logfile=logfile, directory=os.getcwd(), dolog=False, stdin=stdin)

    if retcode == 0:
        # remove temporary files
        os.unlink(logfile)
    else:
        raise RuntimeError, "Error translating PDB"

    return


def xyz_coordinates(pdbin):
    ''' Extract xyz for all atoms '''
    pdb_input = iotbx.pdb.pdb_input(file_name=pdbin)
    hierarchy = pdb_input.construct_hierarchy()
    return _xyz_coordinates(hierarchy)


def _xyz_coordinates(hierarchy):
    res_lst, tmp = [], []

    for residue_group in hierarchy.models()[0].chains()[0].residue_groups():
        for atom_group in residue_group.atom_groups():
            for atom in atom_group.atoms():
                tmp.append(atom.xyz)
            res_lst.append([residue_group.resseq_as_int(), tmp])
            tmp = []

    return res_lst


def xyz_cb_coordinates(pdbin):
    ''' Extract xyz for CA/CB atoms '''
    pdb_input = iotbx.pdb.pdb_input(file_name=pdbin)
    hierarchy = pdb_input.construct_hierarchy()

    res_dict = _xyz_cb_coordinates(hierarchy)

    cb_lst = []
    for i in xrange(len(res_dict)):
        if len(res_dict[i]) > 1:
            cb_lst.append(res_dict[i][1])
        elif len(res_dict[i]) == 1:
            cb_lst.append(res_dict[i][0])

    return cb_lst


def _xyz_cb_coordinates(hierarchy):
    res_lst = []

    for residue_group in hierarchy.models()[0].chains()[0].residue_groups():
        for atom_group in residue_group.atom_groups():
            xyz_lst = _xyz_atom_coords(atom_group)
            res_lst.append([residue_group.resseq_as_int(), xyz_lst])

    return res_lst


def _xyz_atom_coords(atom_group):
    ''' Use this method if you need to identify if CB is present
        in atom_group and if not return CA
    '''

    tmp_dict = {}
    for atom in atom_group.atoms():
        if atom.name.strip() in set(["CA", "CB"]):
            tmp_dict[atom.name.strip()] = atom.xyz

    if 'CB' in tmp_dict: return tmp_dict['CB']
    elif 'CA' in tmp_dict: return tmp_dict['CA']
    else: return (float('inf'), float('inf'), float('inf'))


class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Set up paths. Need to do this with setUpClass, as otherwise the __file__
        variable is updated whenever the cwd is changed in a test and the next test
        gets the wrong paths.
        """
        cls.thisd = os.path.abspath(os.path.dirname(__file__))
        paths = cls.thisd.split(os.sep)
        cls.ample_dir = os.sep.join(paths[:-1])
        cls.tests_dir = os.path.join(cls.ample_dir, "tests")
        cls.testfiles_dir = os.path.join(cls.tests_dir, 'testfiles')
        return

    def testGetInfo1(self):
        """"""

        pdbfile = os.path.join(self.testfiles_dir, "1GU8.pdb")

        info = get_info(pdbfile)

        self.assertEqual(info.pdbCode, "1GU8")
        self.assertEqual(len(info.models), 2)

        m1 = info.models[0]
        self.assertEqual(m1.chains[0], 'A')
        self.assertEqual(m1.resSeqs[0], [
            2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
            31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
            58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
            85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108,
            109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129,
            130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150,
            151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171,
            172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192,
            193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213,
            214, 215, 216, 217, 218, 219
        ])
        self.assertEqual(
            m1.sequences[0],
            'VGLTTLFWLGAIGMLVGTLAFAWAGRDAGSGERRYYVTLVGISGIAAVAYVVMALGVGWVPVAERTVFAPRYIDWILTTPLIVYFLGLLAGLDSREFGIVITLNTVVMLAGFAGAMVPGIERYALFGMGAVAFLGLVYYLVGPMTESASQRSSGIKSLYVRLRNLTVILWAIYPFIWLLGPPGVALLTPTVDVALIVYLDLVTKVGFGFIALDAAATL'
        )

        self.assertEqual(m1.caMask[0], [False] * 218)
        self.assertEqual(m1.bbMask[0], [
            False, True, False, False, False, False, False, False, False, True, False, False, True, False, False, False,
            True, False, False, False, False, False, False, False, True, False, False, False, True, False, True, False,
            False, False, False, False, False, False, False, False, True, False, False, True, False, False, False,
            False, False, False, False, False, False, False, False, True, False, True, False, False, False, False,
            False, False, False, False, False, False, False, False, False, False, False, False, False, False, False,
            False, False, False, False, False, False, False, False, False, True, False, False, False, True, False,
            False, False, False, False, False, True, False, False, False, False, False, False, False, False, False,
            False, False, False, True, False, False, True, False, False, False, False, True, False, False, False, False,
            False, False, False, True, False, True, False, False, False, False, False, True, False, False, False, False,
            False, False, True, False, False, False, False, False, False, False, False, False, False, False, True,
            False, False, False, False, False, False, False, False, False, False, False, False, False, False, False,
            False, False, False, False, False, False, False, False, False, False, True, False, False, True, False,
            False, False, False, False, False, False, False, False, False, False, False, False, False, False, False,
            False, False, False, False, False, False, True, False, True, False, False, False, False, False, False,
            False, False, False, True
        ])

        self.assertEqual(info.numAtoms(modelIdx=0), 1621)
        self.assertEqual(info.numCalpha(modelIdx=0), 218)

        m2 = info.models[1]
        self.assertEqual(m2.chains[0], 'A')
        self.assertEqual(m2.resSeqs[0], [
            2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
            31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
            58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
            85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108,
            109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129,
            130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150,
            151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171,
            172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192,
            193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213,
            214, 215, 216, 217, 218, 219
        ])
        self.assertEqual(
            m2.sequences[0],
            'VGLTTLFWLGAIGMLVGTLAFAWAGRDAGSGERRYYVTLVGISGIAAVAYVVMALGVGWVPVAERTVFAPRYIDWILTTPLIVYFLGLLAGLDSREFGIVITLNTVVMLAGFAGAMVPGIERYALFGMGAVAFLGLVYYLVGPMTESASQRSSGIKSLYVRLRNLTVILWAIYPFIWLLGPPGVALLTPTVDVALIVYLDLVTKVGFGFIALDAAATL'
        )

        self.assertEqual(info.numAtoms(modelIdx=1), 1621)
        self.assertEqual(info.numCalpha(modelIdx=1), 218)

        return

    def testGetInfo2(self):
        """"""

        pdbfile = os.path.join(self.testfiles_dir, "2UUI.pdb")

        info = get_info(pdbfile)

        self.assertEqual(len(info.models), 1)

        m1 = info.models[0]
        self.assertEqual(m1.chains[0], 'A')
        self.assertEqual(m1.resSeqs[0], [i for i in range(-5, 150)])
        self.assertEqual(
            m1.sequences[0],
            'MHHHHHHKDEVALLAAVTLLGVLLQAYFSLQVISARRAFRVSPPLTTGPPEFERVYRAQVNCSEYFPLFLATLWVAGIFFHEGAAALCGLVYLFARLRYFQGYARSAQLRLAPLYASARALWLLVALAALGLLAHFLPAALRAALLGRLRTLLPW'
        )
        self.assertEqual(m1.caMask[0], [False] * 154 + [True])
        self.assertEqual(m1.bbMask[0], [
            False, False, False, False, False, False, False, False, False, False, False, False, False, False, False,
            False, False, False, False, False, True, False, False, False, False, False, False, False, False, False,
            False, False, False, False, False, False, False, False, False, False, False, False, False, False, False,
            False, False, True, False, False, False, False, False, False, False, False, False, False, False, False,
            False, False, False, False, False, False, False, False, False, False, False, False, False, False, False,
            False, True, False, False, False, False, False, True, False, False, False, False, False, True, False, False,
            False, False, False, False, False, False, False, False, False, False, True, False, False, False, False,
            False, False, False, False, False, False, False, False, False, False, False, False, False, False, False,
            False, False, False, False, False, False, False, False, False, True, False, False, False, False, False,
            False, False, False, False, False, False, False, False, False, False, True, False, False, False, False,
            True, True, True, True
        ])

        self.assertEqual(info.numAtoms(modelIdx=0), 1263)

        return

    def testCheckPdbs(self):
        logging.basicConfig()
        logging.getLogger().setLevel(logging.DEBUG)

        pdbs = glob.glob(os.path.join(self.testfiles_dir, "models", "*.pdb"))
        self.assertTrue(check_pdbs(pdbs))

        self.assertFalse(check_pdbs(pdbs, single=True, sequence="AABBCC"))

        pdbs += [os.path.join(self.testfiles_dir, "1GU8.pdb")]
        self.assertFalse(check_pdbs(pdbs, single=True, sequence="AABBCC"))

        return

    def testSelectResidues(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "testSelectResidues1.pdb"
        to_delete = [5, 10, 15, 20]

        b4 = set(resseq(pdbin)['A'])

        select_residues(pdbin=pdbin, pdbout=pdbout, delete=to_delete)

        after = set(resseq(pdbout)['A'])
        self.assertEqual(after, b4.difference(set(to_delete)))

        os.unlink(pdbout)

        return

    def testSelectResiduesKeepIdxs(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "testSelectResidues2.pdb"
        tokeep_idx = [0, 5, 10, 15, 20]

        b4 = [r for i, r in enumerate(resseq(pdbin)['A']) if i in tokeep_idx]
        select_residues(pdbin=pdbin, pdbout=pdbout, tokeep_idx=tokeep_idx)

        #hierachy=iotbx.pdb.pdb_input(pdbout).construct_hierarchy()
        #print [a.name for a in hierachy.models()[0].chains()[0].atoms() ]
        #print "GOT ",resseq(pdbout)
        after = resseq(pdbout)['A']
        self.assertEqual(after, b4)

        os.unlink(pdbout)

        return

    def testSequence1(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        ref = {
            'A': 'GEIAALKQEIAALKKEIAALKEIAALKQGYY',
            'B': 'GEIAALKQEIAALKKEIAALKEIAALKQGYY',
            'C': 'GEIAALKQEIAALKKEIAALKEIAALKQGYY'
        }
        s = sequence(pdbin)
        self.assertEqual(ref, s, "Bad _sequecne: {0}".format(s))
        return

    def XtestSplit(self):
        pdbin = os.path.join(self.testfiles_dir, "1GU8.pdb")
        Xsplit(pdbin)
        #os.unlink(pdbout)
        return

    def testStdResidues(self):

        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "std.pdb"

        std_residues_cctbx(pdbin, pdbout)

        # Check it's valid
        pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdbout)

        #Get list of all the residue names in chain 1
        resnames = [g.unique_resnames()[0] for g in pdb_obj.hierarchy.models()[0].chains()[0].residue_groups()]
        ref = [
            'ACE', 'GLY', 'GLU', 'ILE', 'ALA', 'ALA', 'LEU', 'LYS', 'GLN', 'GLU', 'ILE', 'ALA', 'ALA', 'LEU', 'LYS',
            'LYS', 'GLU', 'ILE', 'ALA', 'ALA', 'LEU', 'LYS', 'PHE', 'GLU', 'ILE', 'ALA', 'ALA', 'LEU', 'LYS', 'GLN',
            'GLY', 'TYR', 'TYR'
        ]
        self.assertEqual(resnames, ref)

        os.unlink(pdbout)

        return

    def testStdResiduesCctbx(self):

        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        pdbout = "std.pdb"

        std_residues_cctbx(pdbin, pdbout)

        # Check it's valid
        pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdbout)

        #Get list of all the residue names in chain 1
        resnames = [g.unique_resnames()[0] for g in pdb_obj.hierarchy.models()[0].chains()[0].residue_groups()]
        ref = [
            'ACE', 'GLY', 'GLU', 'ILE', 'ALA', 'ALA', 'LEU', 'LYS', 'GLN', 'GLU', 'ILE', 'ALA', 'ALA', 'LEU', 'LYS',
            'LYS', 'GLU', 'ILE', 'ALA', 'ALA', 'LEU', 'LYS', 'PHE', 'GLU', 'ILE', 'ALA', 'ALA', 'LEU', 'LYS', 'GLN',
            'GLY', 'TYR', 'TYR'
        ]
        self.assertEqual(resnames, ref)

        os.unlink(pdbout)

        return

    def testStripHetatm(self):
        pdbin = os.path.join(self.testfiles_dir, "1BYZ.pdb")
        pdbout = 'strip_het.pdb'
        hierachy = iotbx.pdb.pdb_input(pdbin).construct_hierarchy()
        _strip(hierachy, hetatm=True, hydrogen=False)
        hierachy.write_pdb_file(pdbout, anisou=False)
        with open(pdbout) as f:
            got = any([True for l in f.readlines() if l.startswith('HETATM')])
        self.assertFalse(got, "Found HETATMS")
        os.unlink(pdbout)
        return

    def testStripHydrogen(self):
        pdbin = os.path.join(self.testfiles_dir, "1BYZ.pdb")
        pdbout = 'strip_H.pdb'
        hierachy = iotbx.pdb.pdb_input(pdbin).construct_hierarchy()
        _strip(hierachy, hetatm=False, hydrogen=True)
        hierachy.write_pdb_file(pdbout, anisou=False)
        with open(pdbout) as f:
            got = any([True for l in f.readlines() if l.startswith('ATOM') and l[13] == 'H'])
        self.assertFalse(got, "Found Hydrogens")
        os.unlink(pdbout)
        return

    def testStripAtomTypes(self):
        pdbin = os.path.join(self.testfiles_dir, "1BYZ.pdb")
        pdbout = 'strip_types.pdb'
        hierachy = iotbx.pdb.pdb_input(pdbin).construct_hierarchy()
        _strip(hierachy, hetatm=False, hydrogen=False, atom_types=['CB'])
        hierachy.write_pdb_file(pdbout, anisou=False)
        with open(pdbout) as f:
            got = any([True for l in f.readlines() if l.startswith('ATOM') and l[12:15].strip() == 'CB'])
        self.assertFalse(got, "Found Atom Types")
        os.unlink(pdbout)
        return

    def testReliableSidechains(self):

        pdbin = os.path.join(self.testfiles_dir, "1GU8.pdb")
        pdbout = "std.pdb"

        reliable_sidechains(pdbin, pdbout)

        # Check it's valid
        pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdbout)

        #Get list of all the residue names in chain 1
        resnames = [g.unique_resnames()[0] for g in pdb_obj.hierarchy.models()[0].chains()[0].residue_groups()]
        ref = [
            'VAL', 'GLY', 'LEU', 'THR', 'THR', 'LEU', 'PHE', 'TRP', 'LEU', 'GLY', 'ALA', 'ILE', 'GLY', 'MET', 'LEU',
            'VAL', 'GLY', 'THR', 'LEU', 'ALA', 'PHE', 'ALA', 'TRP', 'ALA', 'GLY', 'ARG', 'ASP', 'ALA', 'GLY', 'SER',
            'GLY', 'GLU', 'ARG', 'ARG', 'TYR', 'TYR', 'VAL', 'THR', 'LEU', 'VAL', 'GLY', 'ILE', 'SER', 'GLY', 'ILE',
            'ALA', 'ALA', 'VAL', 'ALA', 'TYR', 'VAL', 'VAL', 'MET', 'ALA', 'LEU', 'GLY', 'VAL', 'GLY', 'TRP', 'VAL',
            'PRO', 'VAL', 'ALA', 'GLU', 'ARG', 'THR', 'VAL', 'PHE', 'ALA', 'PRO', 'ARG', 'TYR', 'ILE', 'ASP', 'TRP',
            'ILE', 'LEU', 'THR', 'THR', 'PRO', 'LEU', 'ILE', 'VAL', 'TYR', 'PHE', 'LEU', 'GLY', 'LEU', 'LEU', 'ALA',
            'GLY', 'LEU', 'ASP', 'SER', 'ARG', 'GLU', 'PHE', 'GLY', 'ILE', 'VAL', 'ILE', 'THR', 'LEU', 'ASN', 'THR',
            'VAL', 'VAL', 'MET', 'LEU', 'ALA', 'GLY', 'PHE', 'ALA', 'GLY', 'ALA', 'MET', 'VAL', 'PRO', 'GLY', 'ILE',
            'GLU', 'ARG', 'TYR', 'ALA', 'LEU', 'PHE', 'GLY', 'MET', 'GLY', 'ALA', 'VAL', 'ALA', 'PHE', 'LEU', 'GLY',
            'LEU', 'VAL', 'TYR', 'TYR', 'LEU', 'VAL', 'GLY', 'PRO', 'MET', 'THR', 'GLU', 'SER', 'ALA', 'SER', 'GLN',
            'ARG', 'SER', 'SER', 'GLY', 'ILE', 'LYS', 'SER', 'LEU', 'TYR', 'VAL', 'ARG', 'LEU', 'ARG', 'ASN', 'LEU',
            'THR', 'VAL', 'ILE', 'LEU', 'TRP', 'ALA', 'ILE', 'TYR', 'PRO', 'PHE', 'ILE', 'TRP', 'LEU', 'LEU', 'GLY',
            'PRO', 'PRO', 'GLY', 'VAL', 'ALA', 'LEU', 'LEU', 'THR', 'PRO', 'THR', 'VAL', 'ASP', 'VAL', 'ALA', 'LEU',
            'ILE', 'VAL', 'TYR', 'LEU', 'ASP', 'LEU', 'VAL', 'THR', 'LYS', 'VAL', 'GLY', 'PHE', 'GLY', 'PHE', 'ILE',
            'ALA', 'LEU', 'ASP', 'ALA', 'ALA', 'ALA', 'THR', 'LEU'
        ]

        self.assertEqual(resnames, ref)

        reliable_sidechains_cctbx(pdbin, pdbout)
        pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdbout)
        self.assertEqual(resnames, ref)
        os.unlink(pdbout)

        return

    def testXyzCoordinates(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        test_hierarchy = iotbx.pdb.pdb_input(file_name=pdbin).construct_hierarchy()
        xyz_lst = _xyz_coordinates(test_hierarchy)

        ref_data_start = [(0, [(25.199, 11.913, -9.25), (25.201, 10.666, -9.372),
                               (26.454, 12.702, -9.001)]), (1, [(24.076, 12.643, -9.179), (22.806, 12.124, -9.698),
                                                                (22.170, 11.067, -8.799), (22.404, 11.024, -7.580)]),
                          (2, [(21.377, 10.190, -9.397), (20.675, 9.156, -8.637), (21.614, 8.106, -7.996),
                               (21.337, 7.619, -6.898), (19.625, 8.485, -9.531), (18.637, 7.595, -8.790),
                               (17.652, 8.361, -7.951), (17.724, 9.603, -7.887), (16.786, 7.706, -7.365)])]

        for idx in xrange(len(ref_data_start)):  # Stuff that needs to be true
            self.assertEqual(ref_data_start[idx][0], xyz_lst[idx][0])
            self.assertSequenceEqual(ref_data_start[idx][1], xyz_lst[idx][1])
        nr_atoms = sum(len(i[1]) for i in xyz_lst)
        self.assertEqual(252, nr_atoms)
        self.assertEqual(35, len(xyz_lst))

    def testXyzCbCoordinates(self):
        pdbin = os.path.join(self.testfiles_dir, "4DZN.pdb")
        test_hierarchy = iotbx.pdb.pdb_input(file_name=pdbin).construct_hierarchy()
        xyz_cb_lst = _xyz_cb_coordinates(test_hierarchy)

        ref_data_start = [(0, (float('inf'), float('inf'), float('inf'))), (1, (22.806, 12.124,
                                                                                -9.698)), (2, (19.625, 8.485, -9.531)),
                          (3, (24.783, 6.398, -9.051)), (4, (25.599, 10.846, -6.036)), (5, (20.430, 10.143, -4.644))]

        self.assertSequenceEqual(ref_data_start[1], xyz_cb_lst[1][:6])
        self.assertEqual(35, len(xyz_cb_lst))


if __name__ == "__main__":
    #unittest.TextTestRunner(verbosity=2).run(testSuite())
    #
    # Command-line handling
    #
    #unittest.TextTestRunner(verbosity=2).run(testSuite())
    import argparse
    parser = argparse.ArgumentParser(description='Manipulate PDB files', prefix_chars="-")

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-ren', action='store_true', help="Renumber the PDB")
    group.add_argument('-std', action='store_true', help='Standardise the PDB')
    group.add_argument('-seq', action='store_true', help='Write a fasta of the found AA to stdout')
    group.add_argument('-split_models', action='store_true', help='Split a pdb into constituent models')
    group.add_argument('-split_chains', action='store_true', help='Split a pdb into constituent chains')
    parser.add_argument(
        'input_file',  #nargs='?',
        help='The input file - will not be altered')
    parser.add_argument('-o', dest='output_file', help='The output file - will be created')
    parser.add_argument('-chain', help='The chain to use')
    parser.add_argument('-test', action='store_true', help='Run unittests')

    args = parser.parse_args()

    if args.test:
        print unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])
        sys.exit(unittest.TextTestRunner().run(unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])))

    # Get full paths to all files
    args.input_file = os.path.abspath(args.input_file)
    if not os.path.isfile(args.input_file):
        raise RuntimeError, "Cannot find input file: {0}".format(args.input_file)

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
        print sequence_util.Sequence(pdb=args.input_file).fasta_str()
    elif args.split_models:
        print split_pdb(args.input_file)
    elif args.split_chains:
        print split_into_chains(args.input_file, chain=args.chain)
