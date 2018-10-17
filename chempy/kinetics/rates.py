# -*- coding: utf-8 -*-
"""
This module collects object representing rate expressions. It is based
on the ``chemp.util._expr`` module. The API is somewhat cumbersome since
it tries to be compatible with pure python, SymPy and the underlying
units library of ChemPy (``quantities``). Consider the API to be provisional.
"""

from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict
from functools import reduce
import math
from operator import add

from ..units import get_derived_unit, default_units, energy, concentration
from ..util._dimensionality import dimension_codes, base_registry
from ..util.pyutil import memoize, deprecated
from ..util._expr import Expr


_molar = getattr(default_units, 'molar', 1)  # makes module importable.


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


@memoize(None)
def mk_Radiolytic(*doserate_names):
    """ Create a Radiolytic rate expression

    Note that there is no mass-action dependence in the resulting
    class, i.e. the rates does not depend on any concentrations.

    Parameters
    ----------
    \\*doserate_names : str instances
        Default: ('',)


    Examples
    --------
    >>> RadiolyticAlpha = mk_Radiolytic('alpha')
    >>> RadiolyticGamma = mk_Radiolytic('gamma')
    >>> dihydrogen_alpha = RadiolyticAlpha([0.8e-7])
    >>> dihydrogen_gamma = RadiolyticGamma([0.45e-7])
    >>> RadiolyticAB = mk_Radiolytic('alpha', 'beta')

    Notes
    -----
    The instance __call__ will require by default ``'density'`` and ``'doserate'``
    in variables.

    """
    if len(doserate_names) == 0:
        doserate_names = ('',)

    class _Radiolytic(RadiolyticBase):
        argument_names = tuple('radiolytic_yield{0}'.format('' if drn == '' else '_' + drn)
                               for drn in doserate_names)  # [amount/energy]
        parameter_keys = ('density',) + tuple('doserate{0}'.format('' if drn == '' else '_' + drn)
                                              for drn in doserate_names)

        def args_dimensionality(self, reaction):
            N = base_registry['amount']
            E = get_derived_unit(base_registry, 'energy')
            return (dict(zip(dimension_codes, N/E)),)*self.nargs

        def g_values(self, *args, **kwargs):
            return OrderedDict(zip(self.parameter_keys[1:], self.all_args(*args, **kwargs)))

        @deprecated(use_instead='Radiolytic.all_args')
        def g_value(self, variables, backend=math, **kwargs):
            g_val, = self.all_args(variables, backend=backend, **kwargs)
            return g_val

        def __call__(self, variables, backend=math, reaction=None, **kwargs):
            return variables['density']*reduce(add, [variables[k]*gval for k, gval in zip(
                self.parameter_keys[1:], self.all_args(variables, backend=backend, **kwargs))])

    _Radiolytic.__name__ = 'Radiolytic' if doserate_names == ('',) else ('Radiolytic_' + '_'.join(doserate_names))
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

    def args_dimensionality(self, reaction):
        order = reaction.order()
        return ({'time': -1, 'amount': 1-order, 'length': 3*(order - 1)},)

    def active_conc_prod(self, variables, backend=math, reaction=None):
        result = 1
        for k, v in reaction.reac.items():
            result *= variables[k]**v
        return result

    def rate_coeff(self, variables, backend=math, **kwargs):
        rat_c, = self.all_args(variables, backend=backend, **kwargs)
        return rat_c

    def __call__(self, variables, backend=math, reaction=None, **kwargs):
        return self.rate_coeff(variables, backend=backend, reaction=reaction)*self.active_conc_prod(
            variables, backend=backend, reaction=reaction, **kwargs)

    @classmethod
    def from_callback(cls, callback, attr='rate_coeff', **kwargs):
        return super(MassAction, cls).from_callback(callback, attr=attr, **kwargs)

    def string(self, *args, **kwargs):
        if self.args is None and len(self.unique_keys) == 1:
            return self.unique_keys[0]
        else:
            return super(MassAction, self).string(*args, **kwargs)

    @classmethod
    @deprecated(use_instead='MassAction.from_callback')
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

    def args_dimensionality(self, reaction):
        order = reaction.order()
        return (
            {'time': -1,
             'amount': 1-order, 'length': 3*(order - 1)},
            {'temperature': 1},
        )

    def __call__(self, variables, backend=math, **kwargs):
        A, Ea_over_R = self.all_args(variables, backend=backend, **kwargs)
        try:
            Ea_over_R = Ea_over_R.simplified
        except AttributeError:
            pass
        return A*backend.exp(-Ea_over_R/variables['temperature'])


class Eyring(Expr):
    """ Rate expression for Eyring eq: c0*T*exp(-c1/T)

    Note that choice of standard state (c^0) will matter for order > 1.
    """

    argument_names = ('kB_h_times_exp_dS_R', 'dH_over_R', 'conc0')
    argument_defaults = (1*_molar,)
    parameter_keys = ('temperature',)

    def args_dimensionality(self, reaction):
        order = reaction.order()
        return (
            {'time': -1, 'temperature': -1,
             'amount': 1-order, 'length': 3*(order - 1)},
            {'temperature': 1},
            concentration
        )

    def __call__(self, variables, backend=math, **kwargs):
        c0, c1, conc0 = self.all_args(variables, backend=backend, **kwargs)
        T = variables['temperature']
        return c0*T*backend.exp(-c1/T)*conc0**(1-kwargs['reaction'].order())


class EyringHS(Expr):
    argument_names = ('dH', 'dS', 'c0')
    argument_defaults = (1*_molar,)
    parameter_keys = ('temperature', 'molar_gas_constant',
                      'Boltzmann_constant', 'Planck_constant')

    def args_dimensionality(self, **kwargs):
        return (
            energy + {'amount': -1},
            energy + {'amount': -1, 'temperature': -1},
            concentration
        )

    def __call__(self, variables, backend=math, reaction=None, **kwargs):
        dH, dS, c0 = self.all_args(variables, backend=backend, **kwargs)
        T, R, kB, h = [variables[k] for k in self.parameter_keys]
        return kB/h*T*backend.exp(-(dH-T*dS)/(R*T))*c0**(1-reaction.order())


class RampedTemp(Expr):
    """ Ramped temperature, pass as substitution to e.g. ``get_odesys`` """
    argument_names = ('T0', 'dTdt')
    parameter_keys = ('time',)  # consider e.g. a parameter such as 'init_time'

    def args_dimensionality(self, **kwargs):
        return ({'temperature': 1}, {'temperature': 1, 'time': -1})

    def __call__(self, variables, backend=None, **kwargs):
        T0, dTdt = self.all_args(variables, backend=backend, **kwargs)
        return T0 + dTdt*variables['time']


class SinTemp(Expr):
    argument_names = ('Tbase', 'Tamp', 'angvel', 'phase')
    parameter_keys = ('time',)

    def args_dimensionality(self, **kwargs):
        return ({'temperature': 1}, {'temperature': 1}, {'time': -1}, {})

    def __call__(self, variables, backend=math, **kwargs):
        Tbase, Tamp, angvel, phase = self.all_args(variables, backend=backend, **kwargs)
        return Tbase + Tamp*backend.sin(angvel*variables['time'] + phase)
