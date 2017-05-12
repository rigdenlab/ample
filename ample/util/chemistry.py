"""Storage for chemical information"""

__author__ = "Felix Simkovic"
__date__ = "26 Apr 2017"
__version__ = "0.1"

from cctbx.eltbx import tiny_pse

__all__ = ['atomic_composition', 'periodic_table']


class AtomicComposition(object):

    def __init__(self):
        acids = [
            ('ALA', 'A', {'H':5,  'C':3,  'N':1, 'O':1}),
            ('CYS', 'C', {'H':5,  'C':3,  'N':1, 'O':1, 'S':1}),
            ('ASP', 'D', {'H':4,  'C':4,  'N':1, 'O':3}),
            ('GLU', 'E', {'H':6,  'C':5,  'N':1, 'O':3}),
            ('PHE', 'F', {'H':9,  'C':9,  'N':1, 'O':1}),
            ('GLY', 'G', {'H':3,  'C':2,  'N':1, 'O':1}),
            ('HIS', 'H', {'H':8,  'C':6,  'N':1, 'O':1}),
            ('ILE', 'I', {'H':11, 'C':6,  'N':1, 'O':1}),
            ('LYS', 'K', {'H':13, 'C':6,  'N':2, 'O':1}),
            ('LEU', 'L', {'H':11, 'C':6,  'N':1, 'O':1}),
            ('MET', 'M', {'H':9,  'C':5,  'N':1, 'O':1, 'S':1}),
            ('ASN', 'N', {'H':6,  'C':4,  'N':2, 'O':2}),
            ('PRO', 'P', {'H':7,  'C':5,  'N':1, 'O':1}),
            ('GLN', 'Q', {'H':8,  'C':5,  'N':2, 'O':2}),
            ('ARG', 'R', {'H':13, 'C':6,  'N':4, 'O':1}),
            ('SER', 'S', {'H':5,  'C':3,  'N':1, 'O':2}),
            ('THR', 'T', {'H':7,  'C':4,  'N':1, 'O':2}),
            ('VAL', 'V', {'H':9,  'C':5,  'N':1, 'O':1}),
            ('TRP', 'W', {'H':10, 'C':11, 'N':2, 'O':1}),
            ('TYR', 'Y', {'H':9,  'C':9,  'N':1, 'O':2}),
        ]
        self._aadict = {}
        for three, one, prop in acids:
            self._aadict[three] = self._aadict[one] = _AminoAcidComposition(**prop)

    def __getitem__(self, k):
        if k.upper() in self._aadict:
            return self._aadict[k.upper()]
        return None


class _AminoAcidComposition(object):

    __slots__ = ['H', 'C', 'N', 'O', 'S']

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __repr__(self):
        return "{0}({1})".format(
            self.__class__.__name__, ", ".join(["{0}={1}".format(k, v) for k, v in self.__dict__.items()])
        )


# Instantiate some stuff here so we can call it immediately
atomic_composition = AtomicComposition()
periodic_table = tiny_pse.table
