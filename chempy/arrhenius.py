# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from .util.pyutil import ChemPyDeprecationWarning

import warnings

from .kinetics.arrhenius import (
    arrhenius_equation, fit_arrhenius_equation, ArrheniusParam, ArrheniusParamWithUnits
)


warnings.warn("use .kinetics.arrhenius instead of .arrhenius", ChemPyDeprecationWarning)
