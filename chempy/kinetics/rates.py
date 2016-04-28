# -*- coding: utf-8 -*-
"""
This module collects object representing rate expressions. It is based
on the `chemp.util._expr` module. The API is somewhat cumbersome since
it tries to be compatible with pure python, SymPy and the underlying
units library of SymPy (quantities). Consider the API to be provisional.
"""

from __future__ import (absolute_import, division, print_function)

import math

from ..util._expr import Expr, mk_Poly, mk_PiecewisePoly


class RateExpr(Expr):
    """ Baseclass for rate expressions, see e.g. MassAction & Radiolytic. """

    kw = {'rxn': None, 'ref': None}

    @classmethod
    def subclass_from_callback(cls, cb, cls_attrs=None):
        """ Override RateExpr.__call__

        Parameters
        ----------
        cb : callback
            With signature (variables, all_args, backend) -> scalar
            where `variables` is a dict, `all_args` a tuple and `backend` a module.
        cls_attrs : dict, optional
            Attributes to set in subclass, e.g. parameter_keys, nargs

        Examples
        --------
        >>> from chempy import Reaction
        >>> rxn = Reaction({'O2': 1, 'H2': 1}, {'H2O2': 1})  # d[H2O2]/dt = p0*exp(-p1/T)*sqrt([O2])
        >>> def cb(variables, all_args, backend):
        ...     O2, T = variables['O2'], variables['temperature']
        ...     p0, p1 = all_args
        ...     return p0*backend.sqrt(O2)*backend.exp(-p1/T)
        >>> MyRateExpr = RateExpr.subclass_from_callback(cb, dict(parameter_keys=('temperature',),nargs=2))
        >>> k = MyRateExpr([1.3e9, 4317.2], rxn=rxn)
        >>> print('%.5g' % k({'temperature': 298.15, 'O2': 1.1e-3}))
        22.186

        """
        class _RateExpr(cls):

            def __call__(self, variables, backend=math):
                return cb(variables, self.all_args(variables), backend)
        for k, v in (cls_attrs or {}).items():
            setattr(_RateExpr, k, v)
        return _RateExpr


class MassAction(RateExpr):
    """ Arguments: k """
    nargs = 1

    def rate_coeff(self, variables, backend):  # for subclasses
        return self.arg(variables, 0)

    def __call__(self, variables, backend=math):
        prod = self.rate_coeff(variables, backend)
        for k, v in self.rxn.reac.items():
            prod *= variables[k]**v
        return prod

    @classmethod
    def subclass_from_callback(cls, cb, cls_attrs=None):
        """ Override MassAction.rate_coeff

        Parameters
        ----------
        cb : callback
            With signature (variables, all_args, backend) -> scalar
            where `variables` is a dict, `all_args` a tuple and `backend` a module.
        cls_attrs : dict, optional
            Attributes to set in subclass, e.g. parameter_keys, nargs

        Examples
        --------
        >>> from functools import reduce
        >>> from operator import add
        >>> from chempy import Reaction # d[H2]/dt = 10**(p0 + p1/T + p2/T**2)*[e-]**2
        >>> rxn = Reaction({'e-': 2}, {'OH-': 2, 'H2': 1}, None, {'H2O': 2})
        >>> def cb(variables, all_args, backend):
        ...     T = variables['temperature']
        ...     return 10**reduce(add, [p*T**-i for i, p in enumerate(all_args)])
        >>> MyMassAction = MassAction.subclass_from_callback(cb, dict(parameter_keys=('temperature',), nargs=None))
        >>> k = MyMassAction([9, 300, -75000], rxn=rxn)
        >>> print('%.5g' % k({'temperature': 293., 'e-': 1e-10}))
        1.4134e-11

        """
        class _MassAction(cls):

            def rate_coeff(self, variables, backend=math):
                return cb(variables, self.all_args(variables), backend)
        for k, v in (cls_attrs or {}).items():
            setattr(_MassAction, k, v)
        return _MassAction


class ArrheniusMassAction(MassAction):
    """ Arguments: A, Ea_over_R """
    parameter_keys = ('temperature',)
    nargs = 2

    def rate_coeff(self, variables, backend):
        A, Ea_over_R = self.all_args(variables)
        return A*backend.exp(-Ea_over_R/variables['temperature'])


class EyringMassAction(ArrheniusMassAction):
    """ Arguments: kB_h_times_exp_dS_R, dH_over_R """
    def rate_coeff(self, variables, backend):
        kB_h_times_exp_dS_R, dH_over_R = self.all_args(variables)
        T = variables['temperature']
        return T * kB_h_times_exp_dS_R * backend.exp(-dH_over_R/T)


class Radiolytic(RateExpr):
    """ Arguments: radiolytic_yield [amount/energy] """

    parameter_keys = ('doserate', 'density')
    nargs = 1

    def g_value(self, variables, backend):  # for subclasses
        return self.arg(variables, 0)

    def __call__(self, variables, backend=math):
        return self.g_value(variables, 0)*variables['doserate']*variables['density']


TPoly = mk_Poly('temperature')
RTPoly = mk_Poly('temperature', reciprocal=True)


class TPolyMassAction(TPoly, MassAction):
    """ Arguments: temperature_offset, c0, c1, ... """
    nargs = None
    parameter_keys = TPoly.parameter_keys

    def rate_coeff(self, variables, backend):
        return self.eval_poly(variables, backend)


TPiecewisePoly = mk_PiecewisePoly('temperature')
RTPiecewisePoly = mk_PiecewisePoly('temperature', reciprocal=True)


class PiecewiseTPolyMassAction(TPiecewisePoly, MassAction):
    nargs = None
    parameter_keys = TPoly.parameter_keys

    def rate_coeff(self, variables, backend):
        return self.eval_poly(variables, backend)


class TPolyRadiolytic(TPoly, Radiolytic):
    nargs = None
    parameter_keys = Radiolytic.parameter_keys + TPoly.parameter_keys

    def g_value(self, variables, backend):
        return self.eval_poly(variables, backend)


class RTPolyMassAction(RTPoly, MassAction):
    """ Arguments: temperature_offset, c0, c1, ... """
    parameter_keys = RTPoly.parameter_keys
    nargs = None

    def rate_coeff(self, variables, backend):
        return self.eval_poly(variables, backend)


class _Log10XPolyMassAction(MassAction):
    skip_poly = 1  # kunit

    def rate_coeff(self, variables, backend):
        k_unit = self.arg(variables, 0)
        return 10**self.eval_poly(variables, backend)*k_unit


class Log10TPolyMassAction(TPoly, _Log10XPolyMassAction):
    """ Arguments: k_unit, temperature_offset, c0, c1, ... """
    nargs = None
    parameter_keys = TPoly.parameter_keys
    skip_poly = 1  # kunit


class Log10RTPolyMassAction(RTPoly, _Log10XPolyMassAction):
    """ Arguments: k_unit, temperature_offset, c0, c1, ... """
    nargs = None
    parameter_keys = RTPoly.parameter_keys
    skip_poly = 1  # kunit


class Log10PiecewiseRTPolyMassAction(RTPiecewisePoly, _Log10XPolyMassAction):
    """ Arguments: k_unit, *args """
    nargs = None
    parameter_keys = RTPiecewisePoly.parameter_keys
    skip_poly = 1  # kunit


class TPolyInLog10MassAction(TPoly, MassAction):
    """ Arguments: T_unit, temperature_offset, c0, c1, ... """
    nargs = None
    parameter_keys = TPoly.parameter_keys
    skip_poly = 1  # T_unit

    def rate_coeff(self, variables, backend=math):
        T_u = self.arg(variables, 0)  # T_unit
        new_vars = variables.copy()
        new_vars['temperature'] = backend.log10(variables['temperature'] / T_u)
        return self.eval_poly(new_vars, backend=backend)


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
