# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math
from ..util.expr import Expr, mk_Poly


class RateExpr(Expr):

    kw = {'rxn': None, 'ref': None}


class MassAction(RateExpr):
    """ Arguments: k """
    nargs = 1

    def rate_coeff(self, variables, args, backend):  # for subclasses
        return self.arg(variables, args, 0)

    def __call__(self, variables, args=None, backend=math):
        prod = self.rate_coeff(variables, args, backend)
        for k, v in self.rxn.reac.items():
            prod *= variables[k]**v
        return prod


class ArrheniusMassAction(MassAction):
    """ Arguments: A, Ea_over_R """
    parameter_keys = ('temperature',)
    nargs = 2

    def rate_coeff(self, variables, args, backend):
        A, Ea_over_R = self.all_args(variables, args)
        return A*backend.exp(-Ea_over_R/variables['temperature'])


class Radiolytic(RateExpr):
    """ Arguments: yield [amount/volume] """

    parameter_keys = ('doserate', 'density')
    nargs = 1

    def g_value(self, variables, args, backend):  # for subclasses
        return self.arg(variables, args, 0)

    def __call__(self, variables, args=None, backend=math):
        return self.g_value(variables, args, 0)*variables['doserate']*variables['density']

TPoly = mk_Poly('temperature')
RTPoly = mk_Poly('temperature', reciprocal=True)


class TPolyRadiolytic(TPoly, Radiolytic):
    nargs = None

    def g_value(self, variables, args, backend):
        return self.eval_poly(variables, args, backend)


class TPolyMassAction(TPoly, MassAction):
    """ Arguments: temperature_offset, c0, c1, ... """
    parameter_keys = ('temperature',)
    nargs = None

    def rate_coeff(self, variables, args, backend):
        return self.eval_poly(variables, args, backend)


class RTPolyMassAction(RTPoly, MassAction):
    """ Arguments: temperature_offset, c0, c1, ... """
    parameter_keys = ('temperature',)
    nargs = None

    def rate_coeff(self, variables, args, backend):
        return self.eval_poly(variables, args, backend)


class Log10TPolyMassAction(TPolyMassAction):
    """ Arguments: k_unit, temperature_offset, c0, c1, ... """
    nargs = None

    def rate_coeff(self, variables, args, backend):
        k_unit = self.arg(variables, args, 0)
        return 10**super(Log10TPolyMassAction, self).rate_coeff(
            variables, self.all_args(variables, args)[1:], backend)*k_unit


class TPolyInLog10MassAction(TPolyMassAction):
    """ Arguments: T_unit, temperature_offset, c0, c1, ... """
    nargs = None

    def __call__(self, variables, args=None, backend=math):
        T_unit = self.arg(variables, args, 0)
        new_vars = variables.copy()
        new_vars['temperature'] = backend.log10(variables['temperature'] / T_unit)
        return super(TPolyInLog10MassAction, self).__call__(
            new_vars, self.all_args(variables, args)[1:], backend=backend)


def law_of_mass_action_rates(conc, rsys, variables=None):
    """ Returns a generator of reaction rate expressions

    Rates from the law of mass action (:attr:`Reaction.inact_reac` ignored)
    from a :class:`ReactionSystem`.

    Parameters
    ----------
    conc : array_like
        concentrations (floats or symbolic objects)
    rsys : ReactionSystem instance
        See :class:`ReactionSystem`
    variables : dict (optional)
        to override parameters in the rate expressions of the reactions

    Examples
    --------
    >>> from chempy import ReactionSystem, Reaction
    >>> line, keys = 'H2O -> H+ + OH- ; 1e-4', 'H2O H+ OH-'
    >>> rsys = ReactionSystem([Reaction.from_string(line, keys)], keys)
    >>> next(law_of_mass_action_rates([55.4, 1e-7, 1e-7], rsys))
    0.00554

    """
    for idx_r, rxn in enumerate(rsys.rxns):
        if isinstance(rxn.param, RateExpr):
            if isinstance(rxn.param, MassAction):
                yield rxn.param(variables)
            else:
                raise ValueError("Not mass-action rate in reaction %d" % idx_r)
        else:
            rate = 1
            for substance_key, coeff in rxn.reac.items():
                s_idx = rsys.as_substance_index(substance_key)
                rate *= conc[s_idx]**coeff
            yield rate*rxn.param
