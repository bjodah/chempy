# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math

from ..util.pyutil import deprecated
from ..util._expr import Expr


class MassActionEq(Expr):

    argument_names = ('equilibrium_constant',)

    def active_conc_prod(self, variables, backend=math, equilibrium=None):
        result = None
        for exp_factor, stoichs in [(1, equilibrium.prod), (-1, equilibrium.reac)]:
            for k, v in stoichs.items():
                if result is None:
                    result = variables[k]**(exp_factor*v)
                else:
                    result *= variables[k]**(exp_factor*v)
        return result

    def eq_const(self, variables, backend=math, **kwargs):
        eq_c, = self.all_args(variables, backend=backend, **kwargs)
        return eq_c

    def __call__(self, *args, **kwargs):
        return self.eq_const(*args, **kwargs)

    def equilibrium_equation(self, variables, backend=math, equilibrium=None, **kwargs):
        return self.eq_const(variables, backend=backend, **kwargs) - self.active_conc_prod(
            variables, backend=backend, equilibrium=equilibrium, **kwargs)

    @classmethod
    def from_callback(cls, callback, attr='eq_const', **kwargs):
        return super(MassActionEq, cls).from_callback(callback, attr=attr, **kwargs)


@deprecated(use_instead=MassActionEq)
class EqExpr(Expr):
    """ Baseclass for equilibrium expressions """
    kw = {'eq': None, 'ref': None}


class GibbsEqConst(MassActionEq):
    argument_names = ('dH_over_R', 'dS_over_R')
    parameter_keys = ('temperature',)

    def eq_const(self, variables, backend=math, **kwargs):
        dH_over_R, dS_over_R = self.all_args(variables, backend=backend)
        T, = self.all_params(variables, backend=backend)
        return backend.exp(dS_over_R - dH_over_R/T)
