from __future__ import print_function

__author__ = 'Felix Simkovic'
__date__ = '02 Dec 2016'
__version__ = '0.1'

import logging

logger = logging.getLogger(__name__)


class DynamicDistances(object):
    """Dynamic distance calculator

    Description
    -----------
    The data stored in this class corresponds to Supplementary Table 3
    in [1]_.

    Examples
    --------
    >>> DynamicDistances.cutoff('A', 'Y')
    9.121

    >>> DynamicDistances.percentile('A', 'Y')
    0.443

    References
    ----------
    .. [1] Kamisetty et al. (2013). Assessing the utility of coevolution based residue-residue
       contact predictions in a sequence and structure rich era. PNAS 110(39), 15674-9.

    """

    __slots__ = ()

    _CB_CB_CUTOFF = {
        'A': {
            'A': 5.381,
            'C': 6.057,
            'E': 7.124,
            'D': 6.388,
            'G': 5.201,
            'F': 8.162,
            'I': 6.587,
            'H': 7.591,
            'K': 8.327,
            'M': 7.605,
            'L': 6.707,
            'N': 6.766,
            'Q': 7.583,
            'P': 6.412,
            'S': 5.829,
            'R': 9.365,
            'T': 5.982,
            'W': 9.252,
            'V': 5.854,
            'Y': 9.121,
        },
        'C': {
            'A': 6.057,
            'C': 6.426,
            'E': 7.449,
            'D': 6.985,
            'G': 5.777,
            'F': 9.026,
            'I': 7.476,
            'H': 8.422,
            'K': 8.494,
            'M': 8.265,
            'L': 7.685,
            'N': 7.205,
            'Q': 7.962,
            'P': 7.157,
            'S': 6.59,
            'R': 9.46,
            'T': 6.801,
            'W': 9.752,
            'V': 6.941,
            'Y': 9.362,
        },
        'D': {
            'A': 6.388,
            'C': 6.985,
            'E': 8.945,
            'D': 8.001,
            'G': 6.135,
            'F': 9.111,
            'I': 7.472,
            'H': 8.634,
            'K': 9.306,
            'M': 8.401,
            'L': 7.696,
            'N': 7.672,
            'Q': 8.601,
            'P': 7.321,
            'S': 6.76,
            'R': 10.123,
            'T': 6.971,
            'W': 9.867,
            'V': 6.972,
            'Y': 9.979,
        },
        'E': {
            'A': 7.124,
            'C': 7.449,
            'E': 9.863,
            'D': 8.945,
            'G': 7.036,
            'F': 9.403,
            'I': 7.949,
            'H': 9.454,
            'K': 9.842,
            'M': 8.87,
            'L': 7.696,
            'N': 7.672,
            'Q': 8.601,
            'P': 7.321,
            'S': 6.76,
            'R': 10.123,
            'T': 6.971,
            'W': 9.867,
            'V': 6.972,
            'Y': 9.979,
        },
        'F': {
            'A': 8.162,
            'C': 9.026,
            'E': 9.403,
            'D': 9.111,
            'G': 7.966,
            'F': 10.903,
            'I': 9.602,
            'H': 9.602,
            'K': 9.344,
            'M': 10.253,
            'L': 9.9,
            'N': 9.168,
            'Q': 9.506,
            'P': 8.895,
            'S': 8.694,
            'R': 10.577,
            'T': 9.03,
            'W': 11.758,
            'V': 9.057,
            'Y': 10.999,
        },
        'G': {
            'A': 5.201,
            'C': 5.777,
            'E': 7.036,
            'D': 6.135,
            'G': 4.467,
            'F': 7.966,
            'I': 6.413,
            'H': 7.472,
            'K': 8.216,
            'M': 7.383,
            'L': 6.554,
            'N': 6.321,
            'Q': 7.297,
            'P': 6.14,
            'S': 5.51,
            'R': 9.166,
            'T': 5.619,
            'W': 8.966,
            'V': 5.671,
            'Y': 9.098,
        },
        'H': {
            'A': 7.591,
            'C': 8.422,
            'E': 9.454,
            'D': 8.634,
            'G': 7.472,
            'F': 9.602,
            'I': 8.523,
            'H': 10.606,
            'K': 9.582,
            'M': 9.396,
            'L': 8.676,
            'N': 8.672,
            'Q': 9.391,
            'P': 8.537,
            'S': 8.051,
            'R': 10.879,
            'T': 8.221,
            'W': 10.661,
            'V': 8.179,
            'Y': 10.843,
        },
        'I': {
            'A': 6.587,
            'C': 7.476,
            'E': 7.949,
            'D': 7.472,
            'G': 6.413,
            'F': 9.602,
            'I': 8.096,
            'H': 8.523,
            'K': 8.329,
            'M': 8.874,
            'L': 8.342,
            'N': 7.631,
            'Q': 8.302,
            'P': 7.554,
            'S': 7.142,
            'R': 9.746,
            'T': 7.442,
            'W': 10.47,
            'V': 7.441,
            'Y': 9.719,
        },
        'K': {
            'A': 8.327,
            'C': 8.494,
            'E': 9.842,
            'D': 9.306,
            'G': 8.216,
            'F': 9.344,
            'I': 8.329,
            'H': 9.582,
            'K': 10.662,
            'M': 9.096,
            'L': 8.479,
            'N': 9.319,
            'Q': 9.667,
            'P': 9.198,
            'S': 8.792,
            'R': 11.322,
            'T': 8.715,
            'W': 10.136,
            'V': 8.077,
            'Y': 10.627,
        },
        'L': {
            'A': 6.707,
            'C': 7.685,
            'E': 8.077,
            'D': 7.696,
            'G': 6.554,
            'F': 9.9,
            'I': 8.342,
            'H': 8.676,
            'K': 8.479,
            'M': 9.122,
            'L': 8.522,
            'N': 7.889,
            'Q': 8.48,
            'P': 7.751,
            'S': 7.394,
            'R': 9.852,
            'T': 7.642,
            'W': 10.707,
            'V': 7.633,
            'Y': 9.889,
        },
        'M': {
            'A': 7.605,
            'C': 8.265,
            'E': 8.87,
            'D': 8.401,
            'G': 7.383,
            'F': 10.253,
            'I': 8.874,
            'H': 9.396,
            'K': 9.096,
            'M': 9.53,
            'L': 9.122,
            'N': 8.55,
            'Q': 9.102,
            'P': 8.247,
            'S': 8.01,
            'R': 10.25,
            'T': 8.397,
            'W': 11.11,
            'V': 8.335,
            'Y': 10.4,
        },
        'N': {
            'A': 6.766,
            'C': 7.205,
            'E': 8.485,
            'D': 7.672,
            'G': 6.321,
            'F': 9.168,
            'I': 7.631,
            'H': 8.672,
            'K': 9.319,
            'M': 8.55,
            'L': 7.889,
            'N': 7.682,
            'Q': 8.502,
            'P': 7.497,
            'S': 7.081,
            'R': 10.135,
            'T': 7.159,
            'W': 9.976,
            'V': 7.219,
            'Y': 10.039,
        },
        'P': {
            'A': 6.412,
            'C': 7.157,
            'E': 7.938,
            'D': 7.321,
            'G': 6.14,
            'F': 8.895,
            'I': 7.554,
            'H': 8.537,
            'K': 9.198,
            'M': 8.247,
            'L': 7.751,
            'N': 7.497,
            'Q': 8.308,
            'P': 7.288,
            'S': 6.937,
            'R': 10.266,
            'T': 7.062,
            'W': 9.719,
            'V': 7.063,
            'Y': 9.965,
        },
        'Q': {
            'A': 7.583,
            'C': 7.962,
            'E': 9.328,
            'D': 8.601,
            'G': 7.297,
            'F': 9.506,
            'I': 8.302,
            'H': 9.391,
            'K': 9.667,
            'M': 9.102,
            'L': 8.48,
            'N': 8.502,
            'Q': 9.074,
            'P': 8.308,
            'S': 7.807,
            'R': 10.61,
            'T': 8.055,
            'W': 10.429,
            'V': 8.008,
            'Y': 10.534,
        },
        'R': {
            'A': 9.365,
            'C': 9.46,
            'E': 10.713,
            'D': 10.123,
            'G': 9.166,
            'F': 10.577,
            'I': 9.746,
            'H': 10.879,
            'K': 11.322,
            'M': 10.25,
            'L': 9.852,
            'N': 10.135,
            'Q': 10.61,
            'P': 10.266,
            'S': 9.753,
            'R': 12.05,
            'T': 9.764,
            'W': 11.355,
            'V': 9.513,
            'Y': 11.615,
        },
        'S': {
            'A': 5.829,
            'C': 6.59,
            'E': 7.483,
            'D': 6.76,
            'G': 5.51,
            'F': 8.694,
            'I': 7.142,
            'H': 8.051,
            'K': 8.792,
            'M': 8.01,
            'L': 7.394,
            'N': 7.081,
            'Q': 7.807,
            'P': 6.937,
            'S': 6.19,
            'R': 9.753,
            'T': 6.45,
            'W': 9.77,
            'V': 6.567,
            'Y': 9.594,
        },
        'T': {
            'A': 5.982,
            'C': 6.801,
            'E': 7.628,
            'D': 6.971,
            'G': 5.619,
            'F': 9.03,
            'I': 7.442,
            'H': 8.221,
            'K': 8.715,
            'M': 8.397,
            'L': 7.642,
            'N': 7.159,
            'Q': 8.055,
            'P': 7.062,
            'S': 6.45,
            'R': 9.764,
            'T': 6.676,
            'W': 9.98,
            'V': 6.791,
            'Y': 9.813,
        },
        'V': {
            'A': 5.854,
            'C': 6.941,
            'E': 7.404,
            'D': 6.972,
            'G': 5.671,
            'F': 9.057,
            'I': 7.441,
            'H': 8.179,
            'K': 8.077,
            'M': 8.335,
            'L': 7.633,
            'N': 7.219,
            'Q': 8.008,
            'P': 7.063,
            'S': 6.567,
            'R': 9.513,
            'T': 6.791,
            'W': 10.021,
            'V': 6.759,
            'Y': 9.442,
        },
        'W': {
            'A': 9.252,
            'C': 9.752,
            'E': 10.303,
            'D': 9.867,
            'G': 8.966,
            'F': 11.758,
            'I': 10.47,
            'H': 10.661,
            'K': 10.136,
            'M': 11.11,
            'L': 10.707,
            'N': 9.976,
            'Q': 10.429,
            'P': 9.719,
            'S': 9.77,
            'R': 11.355,
            'T': 9.98,
            'W': 12.806,
            'V': 10.021,
            'Y': 11.807,
        },
        'Y': {
            'A': 9.121,
            'C': 9.362,
            'E': 10.544,
            'D': 9.979,
            'G': 9.098,
            'F': 10.999,
            'I': 9.719,
            'H': 10.843,
            'K': 10.627,
            'M': 10.4,
            'L': 9.889,
            'N': 10.039,
            'Q': 10.534,
            'P': 9.965,
            'S': 9.594,
            'R': 11.615,
            'T': 9.813,
            'W': 11.807,
            'V': 9.442,
            'Y': 11.536,
        },
    }

    _CB_CB_PERCENT = {
        'A': {
            'A': 0.262,
            'C': 0.394,
            'E': 0.34,
            'D': 0.289,
            'G': 0.269,
            'F': 0.26,
            'I': 0.214,
            'H': 0.38,
            'K': 0.55,
            'M': 0.394,
            'L': 0.25,
            'N': 0.349,
            'Q': 0.356,
            'P': 0.399,
            'S': 0.291,
            'R': 0.485,
            'T': 0.378,
            'W': 0.29,
            'V': 0.312,
            'Y': 0.443,
        },
        'C': {
            'A': 0.394,
            'C': 0.178,
            'E': 0.538,
            'D': 0.299,
            'G': 0.129,
            'F': 0.286,
            'I': 0.295,
            'H': 0.203,
            'K': 0.521,
            'M': 0.439,
            'L': 0.206,
            'N': 0.24,
            'Q': 0.347,
            'P': 0.259,
            'S': 0.24,
            'R': 0.491,
            'T': 0.181,
            'W': 0.417,
            'V': 0.173,
            'Y': 0.585,
        },
        'D': {
            'A': 0.289,
            'C': 0.299,
            'E': 0.354,
            'D': 0.392,
            'G': 0.193,
            'F': 0.351,
            'I': 0.341,
            'H': 0.325,
            'K': 0.343,
            'M': 0.361,
            'L': 0.348,
            'N': 0.337,
            'Q': 0.357,
            'P': 0.416,
            'S': 0.323,
            'R': 0.327,
            'T': 0.307,
            'W': 0.475,
            'V': 0.287,
            'Y': 0.676,
        },
        'E': {
            'A': 0.34,
            'C': 0.538,
            'E': 0.389,
            'D': 0.354,
            'G': 0.249,
            'F': 0.512,
            'I': 0.453,
            'H': 0.443,
            'K': 0.434,
            'M': 0.511,
            'L': 0.475,
            'N': 0.423,
            'Q': 0.45,
            'P': 0.475,
            'S': 0.446,
            'R': 0.363,
            'T': 0.409,
            'W': 0.493,
            'V': 0.51,
            'Y': 0.469,
        },
        'F': {
            'A': 0.26,
            'C': 0.286,
            'E': 0.512,
            'D': 0.351,
            'G': 0.219,
            'F': 0.46,
            'I': 0.347,
            'H': 0.542,
            'K': 0.441,
            'M': 0.377,
            'L': 0.26,
            'N': 0.393,
            'Q': 0.451,
            'P': 0.425,
            'S': 0.394,
            'R': 0.738,
            'T': 0.264,
            'W': 0.447,
            'V': 0.246,
            'Y': 0.767,
        },
        'G': {
            'A': 0.269,
            'C': 0.129,
            'E': 0.249,
            'D': 0.193,
            'G': 0.017,
            'F': 0.219,
            'I': 0.179,
            'H': 0.206,
            'K': 0.358,
            'M': 0.255,
            'L': 0.125,
            'N': 0.169,
            'Q': 0.216,
            'P': 0.245,
            'S': 0.153,
            'R': 0.334,
            'T': 0.12,
            'W': 0.239,
            'V': 0.107,
            'Y': 0.267,
        },
        'H': {
            'A': 0.38,
            'C': 0.203,
            'E': 0.443,
            'D': 0.325,
            'G': 0.206,
            'F': 0.542,
            'I': 0.379,
            'H': 0.333,
            'K': 0.714,
            'M': 0.342,
            'L': 0.401,
            'N': 0.289,
            'Q': 0.401,
            'P': 0.457,
            'S': 0.435,
            'R': 0.595,
            'T': 0.417,
            'W': 0.458,
            'V': 0.383,
            'Y': 0.554,
        },
        'I': {
            'A': 0.214,
            'C': 0.295,
            'E': 0.453,
            'D': 0.341,
            'G': 0.179,
            'F': 0.347,
            'I': 0.321,
            'H': 0.379,
            'K': 0.582,
            'M': 0.327,
            'L': 0.261,
            'N': 0.341,
            'Q': 0.406,
            'P': 0.336,
            'S': 0.342,
            'R': 0.557,
            'T': 0.259,
            'W': 0.397,
            'V': 0.242,
            'Y': 0.589,
        },
        'K': {
            'A': 0.55,
            'C': 0.521,
            'E': 0.434,
            'D': 0.343,
            'G': 0.358,
            'F': 0.441,
            'I': 0.582,
            'H': 0.714,
            'K': 0.738,
            'M': 0.611,
            'L': 0.591,
            'N': 0.398,
            'Q': 0.521,
            'P': 0.55,
            'S': 0.445,
            'R': 0.648,
            'T': 0.464,
            'W': 0.47,
            'V': 0.634,
            'Y': 0.704,
        },
        'L': {
            'A': 0.25,
            'C': 0.206,
            'E': 0.475,
            'D': 0.348,
            'G': 0.125,
            'F': 0.26,
            'I': 0.261,
            'H': 0.401,
            'K': 0.591,
            'M': 0.318,
            'L': 0.198,
            'N': 0.279,
            'Q': 0.411,
            'P': 0.317,
            'S': 0.287,
            'R': 0.578,
            'T': 0.19,
            'W': 0.331,
            'V': 0.179,
            'Y': 0.611,
        },
        'M': {
            'A': 0.394,
            'C': 0.439,
            'E': 0.511,
            'D': 0.361,
            'G': 0.255,
            'F': 0.377,
            'I': 0.327,
            'H': 0.342,
            'K': 0.611,
            'M': 0.457,
            'L': 0.318,
            'N': 0.31,
            'Q': 0.498,
            'P': 0.388,
            'S': 0.369,
            'R': 0.641,
            'T': 0.292,
            'W': 0.397,
            'V': 0.295,
            'Y': 0.661,
        },
        'N': {
            'A': 0.349,
            'C': 0.24,
            'E': 0.423,
            'D': 0.337,
            'G': 0.169,
            'F': 0.393,
            'I': 0.341,
            'H': 0.289,
            'K': 0.398,
            'M': 0.31,
            'L': 0.279,
            'N': 0.249,
            'Q': 0.373,
            'P': 0.334,
            'S': 0.305,
            'R': 0.372,
            'T': 0.262,
            'W': 0.458,
            'V': 0.232,
            'Y': 0.586,
        },
        'P': {
            'A': 0.399,
            'C': 0.259,
            'E': 0.475,
            'D': 0.416,
            'G': 0.245,
            'F': 0.425,
            'I': 0.336,
            'H': 0.457,
            'K': 0.55,
            'M': 0.388,
            'L': 0.317,
            'N': 0.334,
            'Q': 0.41,
            'P': 0.339,
            'S': 0.321,
            'R': 0.506,
            'T': 0.32,
            'W': 0.462,
            'V': 0.298,
            'Y': 0.506,
        },
        'Q': {
            'A': 0.356,
            'C': 0.347,
            'E': 0.45,
            'D': 0.357,
            'G': 0.216,
            'F': 0.451,
            'I': 0.406,
            'H': 0.401,
            'K': 0.521,
            'M': 0.498,
            'L': 0.411,
            'N': 0.373,
            'Q': 0.436,
            'P': 0.41,
            'S': 0.408,
            'R': 0.535,
            'T': 0.378,
            'W': 0.49,
            'V': 0.359,
            'Y': 0.547,
        },
        'R': {
            'A': 0.485,
            'C': 0.491,
            'E': 0.363,
            'D': 0.327,
            'G': 0.334,
            'F': 0.738,
            'I': 0.557,
            'H': 0.595,
            'K': 0.648,
            'M': 0.641,
            'L': 0.578,
            'N': 0.372,
            'Q': 0.535,
            'P': 0.506,
            'S': 0.483,
            'R': 0.704,
            'T': 0.477,
            'W': 0.889,
            'V': 0.514,
            'Y': 0.822,
        },
        'S': {
            'A': 0.291,
            'C': 0.24,
            'E': 0.446,
            'D': 0.323,
            'G': 0.153,
            'F': 0.394,
            'I': 0.342,
            'H': 0.435,
            'K': 0.445,
            'M': 0.369,
            'L': 0.287,
            'N': 0.305,
            'Q': 0.408,
            'P': 0.321,
            'S': 0.292,
            'R': 0.483,
            'T': 0.214,
            'W': 0.497,
            'V': 0.205,
            'Y': 0.467,
        },
        'T': {
            'A': 0.378,
            'C': 0.181,
            'E': 0.409,
            'D': 0.307,
            'G': 0.12,
            'F': 0.264,
            'I': 0.259,
            'H': 0.417,
            'K': 0.464,
            'M': 0.292,
            'L': 0.19,
            'N': 0.262,
            'Q': 0.378,
            'P': 0.32,
            'S': 0.214,
            'R': 0.477,
            'T': 0.188,
            'W': 0.315,
            'V': 0.138,
            'Y': 0.43,
        },
        'V': {
            'A': 0.312,
            'C': 0.173,
            'E': 0.51,
            'D': 0.287,
            'G': 0.107,
            'F': 0.246,
            'I': 0.242,
            'H': 0.383,
            'K': 0.634,
            'M': 0.295,
            'L': 0.179,
            'N': 0.232,
            'Q': 0.359,
            'P': 0.298,
            'S': 0.205,
            'R': 0.514,
            'T': 0.138,
            'W': 0.271,
            'V': 0.145,
            'Y': 0.535,
        },
        'W': {
            'A': 0.29,
            'C': 0.417,
            'E': 0.493,
            'D': 0.475,
            'G': 0.239,
            'F': 0.447,
            'I': 0.397,
            'H': 0.458,
            'K': 0.47,
            'M': 0.397,
            'L': 0.331,
            'N': 0.458,
            'Q': 0.49,
            'P': 0.462,
            'S': 0.497,
            'R': 0.889,
            'T': 0.315,
            'W': 0.473,
            'V': 0.271,
            'Y': 0.684,
        },
        'Y': {
            'A': 0.443,
            'C': 0.585,
            'E': 0.469,
            'D': 0.676,
            'G': 0.267,
            'F': 0.767,
            'I': 0.589,
            'H': 0.554,
            'K': 0.704,
            'M': 0.661,
            'L': 0.611,
            'N': 0.586,
            'Q': 0.547,
            'P': 0.506,
            'S': 0.467,
            'R': 0.822,
            'T': 0.43,
            'W': 0.684,
            'V': 0.535,
            'Y': 0.855,
        },
    }

    @classmethod
    def cutoff(cls, x, y):
        """Return the amino acid pair-specific cB-cB cutoff

        Parameters
        ----------
        x : str
           Single-letter amino acid
        y : str
           Single-letter amino acid

        """
        return cls._CB_CB_CUTOFF[x][y]

    @classmethod
    def percentile(cls, x, y):
        """Return 95-97 percentile data

        Parameters
        ----------
        x : str
           Single-letter amino acid
        y : str
           Single-letter amino acid

        """
        return cls._CB_CB_PERCENT[x][y]


class RosettaFunctionConstructs(object):
    """Storage for string formats of different Rosetta energy function constructs

    Description
    -----------
    For more information on the different energy functions, please refer to the
    corresponding references or the official `RosettaCommons documentation`_

    .. _RosettaCommons documentation: https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/constraint-file

    """

    __slots__ = ()

    _ATOMPAIR = "AtomPair {atom1: >2} {res1_seq: >4} {atom2: >2} {res2_seq: >4} "
    _SCALARWEIGHTED = "SCALARWEIGHTEDFUNC {scalar_score: .3f} "

    @property
    def BOUNDED_default(self):
        """Simple bounded energy function"""
        construct = RosettaFunctionConstructs._ATOMPAIR + "BOUNDED {lower_bound: >.3f} {upper_bound: >.3f} 1 0.5 #"
        return construct

    @property
    def BOUNDED_gremlin(self):
        """Energy function according to [2]_

        References
        ----------
        .. [2] Ovchinnekov et al. (2015). Large-scale determination of previously unsolved
           protein structures using evolutionary information. Elife 3(4), e09248.

        """
        construct = (
            RosettaFunctionConstructs._ATOMPAIR
            + RosettaFunctionConstructs._SCALARWEIGHTED
            + "BOUNDED 0 {lower_bound: >.3f} 1 0.5"
        )
        return construct

    @property
    def FADE(self):
        """Energy function according to [3]_ and [4]_

        References
        ----------
        .. [3] Simkovic et al. (2016). Residue contacts predicted by evolutionary covariance
           extend the application of ab initio molecular replacement to larger and more
           challenging protein folds. IUCrJ 3(Pt 4), 259-270.
        .. [4] Michel et al. (2014). PconsFold: improved contact predictions improve protein
           models. Bioinformatics 30(17), i482-i488

        """
        construct = RosettaFunctionConstructs._ATOMPAIR + "FADE -10 19 10 {energy_bonus: >5.2f} 0"
        return construct

    @property
    def FADE_default(self):
        """Energy function according to [4]_

        References
        ----------
        .. [4] Michel et al. (2014). PconsFold: improved contact predictions improve protein
           models. Bioinformatics 30(17), i482-i488

        """
        construct = RosettaFunctionConstructs._ATOMPAIR + "FADE -10 19 10 -15.00 0"
        return construct

    @property
    def FLAT_HARMONIC(self):
        """ROSETTA FLAT_HARMONIC energy function"""
        construct = RosettaFunctionConstructs._ATOMPAIR + "FLAT_HARMONIC {x0: >.2f} {stddev: >.2f} {tol: >.2f}"
        return construct

    @property
    def GAUSSIAN(self):
        """Simple Gaussian energy function"""
        construct = RosettaFunctionConstructs._ATOMPAIR + "GAUSSIANFUNC {mean: >.2f} {stddev: >.2f} TAG"
        return construct

    @property
    def SIGMOID_default(self):
        """Simple sigmoidal energy function"""
        construct = RosettaFunctionConstructs._ATOMPAIR + "SIGMOID 8.00 1.00 #ContactMap: {raw_score: >4.3f}"
        return construct

    @property
    def SIGMOID_gremlin(self):
        """Energy function according to [2]_

        References
        ----------
        .. [2] Ovchinnekov et al. (2015). Large-scale determination of previously unsolved
           protein structures using evolutionary information. Elife 4, e09248.

        """
        construct = (
            RosettaFunctionConstructs._ATOMPAIR
            + RosettaFunctionConstructs._SCALARWEIGHTED
            + "SUMFUNC 2 SIGMOID {sigmoid_cutoff: >6.3f} {sigmoid_slope: >6.3f} CONSTANTFUNC -0.5"
        )
        return construct


class Saint2FunctionConstructs(object):
    """Storage for string formats of different SAINT2 energy function constructs

    """

    __slots__ = ()

    @property
    def DEFAULT(self):
        """Default format"""
        # Actually a contact format but keep it as restraints style
        construct = "{res1_seq: >4} {res2_seq: >4} {raw_score: >4.3f}"
        return construct
