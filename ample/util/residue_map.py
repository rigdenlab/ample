"""Useful manipulations on PDB files"""

import types

from ample.util import ample_util, pdb_model


class residueSequenceMap(object):
    """Class for handling mapping between model and native residue indices."""

    def __init__(self, refPdb=None, targetPdb=None):

        self.refResSeq = []
        self.refSequence = None
        self.refCAlphaMask = []
        self.refBbMask = []
        self.refOffset = None
        self._refIncomparable = None  # List of atoms in the model that cannot be compared to the model

        self.targetResSeq = []
        self.targetSequence = None
        self.targetCAlphaMask = []
        self.targetBbMask = []
        self.targetOffset = None  # Where the matched part of the sequences starts in the model
        self._targetIncomparable = None  # List of atoms in the model that cannot be compared to the native

        self.lenMatch = None

        # The window of AA we used to check for a match
        self.probeLen = 10

        # maxInset is the max number of AA into the sequence that we will go searching for a match - i.e. if more
        # then maxInset AA at the start are non-matching, we won't find the match
        self.maxInset = 50

        # Like this just for testing
        if refPdb and targetPdb:
            self.calc_map(refPdb, targetPdb)

        return

    def ref2target(self, refResSeq):
        """Return the target resSeq for the given reference resSeq.
        This will calculate a resSeq in the target if there isn't one.
        """

        # Work out how many residues from the start of the matching region this residue is in the target
        indent = refResSeq - self.refResSeq[self.refOffset]

        # calculate the corresponding index in the reference
        targetResSeq = self.targetResSeq[self.targetOffset] + indent

        ## paranoid check
        # if 0 < self.targetOffset + indent < len( self.targetResSeq ):
        #    assert targetResSeq == self.targetResSeq[ self.targetOffset + indent ]

        return targetResSeq

    def target2ref(self, targetResSeq):
        """Return the referece resSeq for the given target resSeq.
        This will calculate a resSeq in the reference if there isn't one.
        """

        # Work out how many residues from the start of the matching region this residue is in the target
        indent = targetResSeq - self.targetResSeq[self.targetOffset]

        refResSeq = self.refResSeq[self.refOffset] + indent

        ## paranoid check
        # if 0 < self.refOffset + indent < len( self.refResSeq ):
        #    assert refResSeq == self.refResSeq[ self.refOffset + indent ]

        return refResSeq

    def targetIncomparable(self, cAlphaMask=True, bbMask=False):
        """Return a list of the resSeq in the target that cannot be compared to the reference.
        This includes any where there isn't a corresponding residue in the reference, or there isn't a c-alpha
        or backbone atom in either (if cAlphaMask or bbMask is set)
        """

        if self._targetIncomparable == None:

            self._targetIncomparable = []
            for i, resSeq in enumerate(self.targetResSeq):

                # Before the start of the matching region
                if i < self.targetOffset:
                    self._targetIncomparable.append(resSeq)
                    continue

                # After end of matching region
                if i > self.lenMatch + 1:
                    self._targetIncomparable.append(resSeq)
                    continue

                # In matching region but no C-alpha
                if cAlphaMask:
                    if self.targetCAlphaMask[i]:
                        self._targetIncomparable.append(resSeq)
                        continue

                # In matching region but no complete bbatoms
                if bbMask:
                    if self.targetBbMask[i]:
                        self._targetIncomparable.append(resSeq)
                        continue

                # Matching residues in reference
                refResSeq = self.target2ref(resSeq)
                try:
                    j = self.refResSeq.index(refResSeq)
                except ValueError:
                    # A residue that isn't actually in the reference
                    self._targetIncomparable.append(resSeq)
                    continue

                # No C-Alpha
                if cAlphaMask:
                    if self.refCAlphaMask[j]:
                        self._targetIncomparable.append(resSeq)
                        continue

                # No bbMask
                if bbMask:
                    if self.refBbMask[j]:
                        self._targetIncomparable.append(resSeq)
                        continue

        return self._targetIncomparable

    def refIncomparable(self, cAlphaMask=True, bbMask=False):
        """Return a list of the resSeq in the reference that cannot be compared to the target.
        This includes any where there isn't a corresponding residue in the target, or there isn't a c-alpha
        or backbone atom in either (if cAlphaMask or bbMask is set)
        """

        if self._refIncomparable == None:
            self._refIncomparable = []
            for i, resSeq in enumerate(self.refResSeq):

                # Before the start of the matching region
                if i < self.refOffset:
                    self._refIncomparable.append(resSeq)
                    continue

                # After end of matching region
                if i > self.lenMatch + 1:
                    self._refIncomparable.append(resSeq)
                    continue

                # In matching region but no C-alpha
                if cAlphaMask:
                    if self.refCAlphaMask[i]:
                        self._refIncomparable.append(resSeq)
                        continue

                # In matching region but no complete bbatoms
                if bbMask:
                    if self.refBbMask[i]:
                        self._refIncomparable.append(resSeq)
                        continue

                # Matching residues in reference
                targetResSeq = self.ref2target(resSeq)
                try:
                    j = self.targetResSeq.index(targetResSeq)
                except ValueError:
                    # A residue that isn't actually in the reference
                    self._refIncomparable.append(resSeq)
                    continue

                # No C-Alpha
                if cAlphaMask:
                    if self.targetCAlphaMask[j]:
                        self._refIncomparable.append(resSeq)
                        continue

                # No bbMask
                if bbMask:
                    if self.targetBbMask[j]:
                        self._refIncomparable.append(resSeq)
                        continue

        return self._refIncomparable

    def __str__(self):
        me = {}
        for slot in dir(self):
            attr = getattr(self, slot)
            if not slot.startswith("__") and not (
                isinstance(attr, types.MethodType) or isinstance(attr, types.FunctionType)
            ):
                me[slot] = attr

        s = self.__repr__() + "\n"
        for k in sorted(me.keys()):
            s += "{0}: {1}\n".format(k, me[k])
        return s

    def calc_map(self, nativePdb, modelPdb):

        self.refSequence, self.refResSeq, self.refCAlphaMask = self.read_pdb(nativePdb)
        self.targetSequence, self.targetResSeq, self.targetCAlphaMask = self.read_pdb(modelPdb)

        self._calc_map()

        return

    def _calc_map(self):
        """Return a ResSeqMap mapping the index of a residue in the model to the corresponding residue in the native.
        Only works if 1 chain in either file and with standard residues
        """
        self.refOffset, self.targetOffset = self._calcOffset(self.refSequence, self.targetSequence)
        self.lenMatch = self._lenMatch()
        self._checkContiguous()

    def _checkContiguous(self):
        """
        Check that there is congiuous resSeq numbering as otherwise the map in current form won't work
        Can only check from the start while residues are matching - won't work after the first gap
        """
        p_refResSeq = p_targetResSeq = None
        for count in range(self.lenMatch):

            refIdx = self.refOffset + count
            refAA = self.refSequence[refIdx]
            refResSeq = self.refResSeq[refIdx]

            targetIdx = self.targetOffset + count
            targetAA = self.targetSequence[targetIdx]
            targetResSeq = self.targetResSeq[targetIdx]

            if count != 0:

                # Need to stop after the first gap
                if refAA != targetAA:
                    break

                # Complain if the residues match but the numbering doesn't
                if refResSeq != p_refResSeq + 1 or targetResSeq != p_targetResSeq + 1:
                    raise RuntimeError(
                        "Non-contiguous residue numbering: {}->{} and {}->{}".format(
                            p_refResSeq, refResSeq, p_targetResSeq, targetResSeq
                        )
                    )
            p_refResSeq = refResSeq
            p_targetResSeq = targetResSeq

    def _calcOffset(self, refSequence, targetSequence, reverse=False):

        # Probe length dependant on protein size - is length of shortest protein or the default
        shortest = min(len(targetSequence), len(refSequence))
        probeLen = min(self.probeLen, shortest)

        # Work out how far we can probe into each sequence
        if len(targetSequence) - probeLen >= self.maxInset:
            targetMaxInset = self.maxInset
        else:
            targetMaxInset = len(targetSequence) - probeLen

        if len(refSequence) - probeLen >= self.maxInset:
            refMaxInset = self.maxInset
        else:
            refMaxInset = len(refSequence) - probeLen

        # If checking from the end, reverse the strings
        if reverse:
            refSequence = refSequence[::-1]
            targetSequence = targetSequence[::-1]

        got = False
        for targetOffset in range(targetMaxInset + 1):
            probe = targetSequence[targetOffset : targetOffset + probeLen]
            for refOffset in range(refMaxInset + 1):
                if refSequence[refOffset : refOffset + probeLen] == probe:
                    got = True
                    break

            if got:
                break

        if not got:
            raise RuntimeError("Could not calculate map for:\n{}\n{}".format(refSequence, targetSequence))

        return refOffset, targetOffset

    def _lenMatch(self):
        refBackOffset, targetBackOffset = self._calcOffset(self.refSequence, self.targetSequence, reverse=True)
        # Calculate match from the residue numbers - use reference for now
        length = self.refResSeq[len(self.refSequence) - 1 - refBackOffset] - self.refResSeq[self.refOffset] + 1

        if length > len(self.refResSeq) and length > len(self.targetResSeq):
            raise RuntimeError("Match of {} is longer than both of the sequences!".format(length))

        return length

    def fromInfo(self, refInfo=None, refChainID=None, targetInfo=None, targetChainID=None, modelIdx=0):
        """Create a map from 2 info objects"""

        # Determine index of chain so we know where to get the data from
        nativeIdx = refInfo.models[modelIdx].chains.index(refChainID)

        self.refResSeq = refInfo.models[modelIdx].resSeqs[nativeIdx]
        self.refSequence = refInfo.models[modelIdx].sequences[nativeIdx]
        self.refCAlphaMask = refInfo.models[modelIdx].caMask[nativeIdx]
        self.refBbMask = refInfo.models[modelIdx].bbMask[nativeIdx]
        self.refOffset = None
        self._refIncomparable = None

        targetIdx = targetInfo.models[modelIdx].chains.index(targetChainID)
        self.targetResSeq = targetInfo.models[modelIdx].resSeqs[targetIdx]
        self.targetSequence = targetInfo.models[modelIdx].sequences[targetIdx]
        self.targetCAlphaMask = targetInfo.models[modelIdx].caMask[targetIdx]
        self.targetBbMask = targetInfo.models[modelIdx].bbMask[targetIdx]
        self.targetOffset = None
        self._targetIncomparable = None

        self._calc_map()

    def read_pdb(self, pdb):
        """Get sequence as string of 1AA
        get list of matching resSeq
        """

        atomTypes = []

        resSeq = []
        resName = []
        _atomTypes = []
        atomTypesList = []

        chain = None
        readingResSeq = None
        readingResName = None
        for line in open(pdb):

            if line.startswith("MODEL"):
                raise RuntimeError("FOUND MULTI_MODEL FILE!")

            if line.startswith("TER"):
                break

            if line.startswith("ATOM"):

                atom = pdb_model.PdbAtom(line)

                if not chain:
                    chain = atom.chainID

                if atom.chainID != chain:
                    raise RuntimeError("FOUND ADDITIONAL CHAIN")

                # First atom in first residue
                if readingResSeq == None:
                    readingResSeq = atom.resSeq
                    readingResName = atom.resName
                    _atomTypes.append(atom.name.strip())
                    continue

                if readingResSeq != atom.resSeq:
                    resName.append(readingResName)
                    resSeq.append(readingResSeq)
                    atomTypesList.append(_atomTypes)

                    readingResSeq = atom.resSeq
                    readingResName = atom.resName
                    _atomTypes = [atom.name.strip()]
                else:
                    if atom.name not in _atomTypes:
                        _atomTypes.append(atom.name.strip())

        resName.append(readingResName)
        resSeq.append(readingResSeq)
        atomTypesList.append(_atomTypes)

        sequence = ""
        for n in resName:
            sequence += ample_util.three2one[n]

        cAlphaMask = ['CA' not in atomTypes for atomTypes in atomTypesList]
        return sequence, resSeq, cAlphaMask

    def resSeqMatch(self):
        """Return true if the residue numbering between the model and native over the aligned region is the same"""
        return (
            self.targetResSeq[self.targetOffset : self.targetOffset + self.lenMatch]
            == self.refResSeq[self.refOffset : self.refOffset + self.lenMatch]
        )
