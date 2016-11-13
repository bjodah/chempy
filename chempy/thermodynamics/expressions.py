# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math

from ..util._expr import Expr


class EqExpr(Expr):
    """ Baseclass for equilibrium expressions """
    kw = {'eq': None, 'ref': None}


class GibbsEqConst(EqExpr):
    argument_names = ('dH_over_R', 'dS_over_R', 'dCp_over_R', 'Tref')
    argument_defaults = (0, 298.15)
    parameter_keys = ('temperature',)

    def __call__(self, variables, backend=math):
        dH_over_R, dS_over_R, dCp_over_R, Tref = self.all_args(variables, backend=backend)
        T, = self.all_params(variables, backend=backend)
        _dH_over_R = dH_over_R + dCp_over_R*(T-Tref)
        _dS_over_R = dS_over_R + dCp_over_R*backend.log(T/Tref)
        return backend.exp(_dS_over_R - _dH_over_R/T)
