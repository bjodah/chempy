# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math

from ..util._expr import Expr


class EqExpr(Expr):
    """ Baseclass for equilibrium expressions """
    kw = {'eq': None, 'ref': None}


class GibbsEnergyExpr(EqExpr):
    argument_names = ('dH_over_R', 'dS_over_R')
    parameter_keys = ('temperature',)

    def __call__(self, variables, backend=math):
        dH_over_R, dS_over_R = self.all_args(variables, backend=backend)
        T, = self.all_params(variables, backend=backend)
        return backend.exp(dS_over_R - dH_over_R/T)
