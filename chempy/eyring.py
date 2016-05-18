# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from .util.pyutil import ChemPyDeprecationWarning

import warnings

from .kinetics.eyring import (
    eyring_equation, EyringParam, EyringParamWithUnits
)


warnings.warn("use .kinetics.eyring instead of .eyring", ChemPyDeprecationWarning)
