# -*- coding: utf-8 -*-
"""
This module collects object representing rate expressions. It is based
on the ``chemp.util._expr`` module. The API is somewhat cumbersome since
it tries to be compatible with pure python, SymPy and the underlying
units library of ChemPy (``quantities``). Consider the API to be provisional.
"""

from __future__ import (absolute_import, division, print_function)

import math


from ..util.pyutil import memoize, deprecated
from ..util._expr import Expr


class RateExpr(Expr):
    """ Baseclass for rate expressions, see source code of e.g. MassAction & Radiolytic. """

    kw = {'rxn': None, 'ref': None}

    @property
    def rxn(self):
        return self._rxn

    @rxn.setter
    def rxn(self, value):
        self._rxn = value
        for arg in self.args:
            if isinstance(arg, RateExpr):
                arg.rxn = value

    def _recursive_as_RateExpr(self):
        new_args = []
        for arg in self.args:
            if isinstance(arg, Expr):
                new_args.append(arg)
            else:
                if hasattr(arg, '_as_RateExpr'):
                    new_args.append(arg._as_RateExpr(self.rxn))
                else:
                    new_args.append(arg)
        if self.kw is None:
            kw = {}
        else:
            kw = {k: getattr(self, k) for k in self.kw}
        return self.__class__(new_args, self.unique_keys, **kw)

    @classmethod
    @deprecated(use_instead=Expr.from_callback)
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
                return cb(variables, self.all_args(variables), backend=backend)
        for k, v in (cls_attrs or {}).items():
            setattr(_RateExpr, k, v)
        return _RateExpr


class RadiolyticBase(RateExpr):
    pass  # for isinstance checks


@memoize(1)
def mk_Radiolytic(doserate_name='doserate'):
    """ Create a Radiolytic rate expression

    Note that there is no mass-action dependence in the resulting
    class, i.e. the rates does not depend on any concentrations.

    Examples
    --------
    >>> RadiolyticAlpha = mk_Radiolytic('doserate_alpha')
    >>> RadiolyticGamma = mk_Radiolytic('doserate_gamma')
    >>> dihydrogen_alpha = RadiolyticAlpha([0.8e-7])
    >>> dihydrogen_gamma = RadiolyticGamma([0.45e-7])

    """
    class _Radiolytic(RadiolyticBase):
        argument_names = ('radiolytic_yield',)  # [amount/energy]
        parameter_keys = (doserate_name, 'density')
        print_name = 'Radiolytic' if doserate_name == 'doserate' else ('Radiolytic{'+doserate_name+'}')

        def g_value(self, variables, backend=math):  # for subclasses
            return self.arg(variables, 0, backend=backend)

        def __call__(self, variables, backend=math):
            return self.g_value(variables, 0)*variables[doserate_name]*variables['density']
    return _Radiolytic


Radiolytic = mk_Radiolytic()


class MassAction(RateExpr):
    """ Rate-expression of mass-action type

    Examples
    --------
    >>> ma = MassAction([3.14])
    >>> ma.rate_coeff({})
    3.14
    >>> from chempy import Reaction
    >>> r = Reaction.from_string('3 A -> B', param=ma)
    >>> r.rate({'A': 2}) == {'A': -75.36, 'B': 25.12}
    True

    """

    argument_names = ('rate_constant',)

    def rate_coeff(self, variables, backend=math):  # for subclasses
        return self.arg(variables, 0, backend=backend)

    def __call__(self, variables, backend=math):
        prod = self.rate_coeff(variables, backend=backend)
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
        >>> MyMassAction = MassAction.subclass_from_callback(cb, dict(parameter_keys=('temperature',), nargs=-1))
        >>> k = MyMassAction([9, 300, -75000], rxn=rxn)
        >>> print('%.5g' % k({'temperature': 293., 'e-': 1e-10}))
        1.4134e-11

        """
        class _MassAction(cls):

            def rate_coeff(self, variables, backend=math):
                return cb(variables, self.all_args(variables), backend=backend)
        for k, v in (cls_attrs or {}).items():
            setattr(_MassAction, k, v)
        return _MassAction

    def as_mass_action(self, variables, backend=math):
        return MassAction([self.rate_coeff(variables, backend=backend)], self.unique_keys, **self.kwargs)


class ArrheniusMassAction(MassAction):
    """ Rate expression for a Arrhenius-type of rate

    Examples
    --------
    >>> from math import exp
    >>> from chempy import Reaction
    >>> from chempy.units import allclose, default_units as u
    >>> A = 1e11 / u.second
    >>> Ea_over_R = 42e3/8.3145 * u.K**-1
    >>> ratex = ArrheniusMassAction([A, Ea_over_R])
    >>> rxn = Reaction({'R'}, {'P'}, ratex)
    >>> dRdt = rxn.rate({'R': 3*u.M, 'temperature': 298.15*u.K})['R']
    >>> allclose(dRdt, -3*1e11*exp(-42e3/8.3145/298.15)*u.M/u.s)
    True

    """
    argument_names = ('A', 'Ea_over_R')
    parameter_keys = ('temperature',)

    def rate_coeff(self, variables, backend=math):
        A, Ea_over_R = self.all_args(variables, backend=backend)
        return A*backend.exp(-Ea_over_R/variables['temperature'])


class EyringMassAction(ArrheniusMassAction):
    argument_names = ('kB_h_times_exp_dS_R', 'dH_over_R')

    def rate_coeff(self, variables, backend=math):
        kB_h_times_exp_dS_R, dH_over_R = self.all_args(variables, backend=backend)
        T = variables['temperature']
        return T * kB_h_times_exp_dS_R * backend.exp(-dH_over_R/T)
