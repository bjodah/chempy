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
        >>> k = MyRateExpr([1.3e9, 4317.2])
        >>> print('%.5g' % k({'temperature': 298.15, 'O2': 1.1e-3, 'rxn': rxn}))
        22.186

        """
        class _RateExpr(cls):

            def __call__(self, variables, backend=math, **kwargs):
                return cb(variables, self.all_args(variables), backend=backend, **kwargs)
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

    Notes
    -----
    The instance __call__ will require ``'density'`` and ``doserate_name``
    in variables.

    """
    class _Radiolytic(RadiolyticBase):
        argument_names = ('radiolytic_yield',)  # [amount/energy]
        parameter_keys = (doserate_name, 'density')

        def g_value(self, variables, backend=math, **kwargs):
            g_val, = self.all_args(variables, backend=backend, **kwargs)
            return g_val

        def __call__(self, variables, backend=math, reaction=None, **kwargs):
            g_value, = self.all_args(variables, backend=backend, **kwargs)
            return self.g_value(variables, backend=backend)*variables[doserate_name]*variables['density']

    _Radiolytic.__name__ = 'Radiolytic' if doserate_name == 'doserate' else ('Radiolytic{'+doserate_name+'}')
    return _Radiolytic


Radiolytic = mk_Radiolytic()


class MassAction(RateExpr):
    """ Rate-expression of mass-action type

    Notes
    -----
    :meth:`__call__` requires a :class:`Reaction` instance to be passed as ``reaction``
    keyword argument.

    Examples
    --------
    >>> ma = MassAction([3.14])
    >>> from chempy import Reaction
    >>> r = Reaction.from_string('3 A -> B', param=ma)
    >>> r.rate({'A': 2}) == {'A': -75.36, 'B': 25.12}
    True

    """

    argument_names = ('rate_constant',)

    def active_conc_prod(self, variables, backend=math, reaction=None):
        result = None
        for k, v in reaction.reac.items():
            if result is None:
                result = variables[k]**v
            else:
                result *= variables[k]**v
        return result

    def rate_coeff(self, variables, backend=math, **kwargs):
        rat_c, = self.all_args(variables, backend=backend, **kwargs)
        return rat_c

    def __call__(self, variables, backend=math, reaction=None, **kwargs):
        return self.rate_coeff(variables, backend=backend)*self.active_conc_prod(
            variables, backend=backend, reaction=reaction, **kwargs)

    @classmethod
    def from_callback(cls, callback, attr='rate_coeff', **kwargs):
        return super(MassAction, cls).from_callback(callback, attr=attr, **kwargs)

    @classmethod
    @deprecated(use_instead=Expr.from_callback)
    def subclass_from_callback(cls, cb, cls_attrs=None):
        """ Override MassAction.__call__ """
        _RateExpr = super(MassAction, cls).subclass_from_callback(cb, cls_attrs=cls_attrs)

        def wrapper(*args, **kwargs):
            obj = _RateExpr(*args, **kwargs)
            return cls(obj)
        return wrapper


class Arrhenius(Expr):
    """ Rate expression for a Arrhenius-type of rate: c0*exp(-c1/T)

    Examples
    --------
    >>> from math import exp
    >>> from chempy import Reaction
    >>> from chempy.units import allclose, default_units as u
    >>> A = 1e11 / u.second
    >>> Ea_over_R = 42e3/8.3145 * u.K**-1
    >>> ratex = MassAction(Arrhenius([A, Ea_over_R]))
    >>> rxn = Reaction({'R'}, {'P'}, ratex)
    >>> dRdt = rxn.rate({'R': 3*u.M, 'temperature': 298.15*u.K})['R']
    >>> allclose(dRdt, -3*1e11*exp(-42e3/8.3145/298.15)*u.M/u.s)
    True

    """
    argument_names = ('A', 'Ea_over_R')
    parameter_keys = ('temperature',)

    def __call__(self, variables, backend=math, **kwargs):
        A, Ea_over_R = self.all_args(variables, backend=backend, **kwargs)
        return A*backend.exp(-Ea_over_R/variables['temperature'])


class Eyring(Expr):
    """ Rate expression for Eyring eq: c0*T*exp(-c1/T) """

    argument_names = ('kB_h_times_exp_dS_R', 'dH_over_R')

    def __call__(self, variables, backend=math, **kwargs):
        c0, c1 = self.all_args(variables, backend=backend, **kwargs)
        T = variables['temperature']
        return c0*T*backend.exp(-c1/T)


class RampedTemp(Expr):
    """ Ramped temperature, pass as substitution to e.g. ``get_odesys`` """
    argument_names = ('T0', 'dTdt')
    parameter_keys = ('time',)  # consider e.g. a parameter such as 'init_time'

    def __call__(self, variables, backend=None, **kwargs):
        T0, dTdt = self.all_args(variables, backend=backend, **kwargs)
        return T0 + dTdt*variables['time']
