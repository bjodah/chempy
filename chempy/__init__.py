# -*- coding: utf-8 -*-
"""
ChemPy is a Python package useful for solving problems in chemistry.
"""

from __future__ import absolute_import, division, print_function

from ._url import __url__
from ._release import __version__
from .chemistry import Substance, Reaction, Equilibrium, Species, balance_stoichiometry, mass_fractions
from .reactionsystem import ReactionSystem
from .henry import Henry
from .util.periodic import atomic_number
from .kinetics import EyringParam, EyringHS, MassAction

from .util.pyutil import ChemPyDeprecationWarning

import sys

if sys.version_info < (3, 5, 0):
    import warnings
    warnings.warn("Use 'chempy<0.7' if using python versions < 3.5", ChemPyDeprecationWarning)
