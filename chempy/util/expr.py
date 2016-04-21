# -*- coding: utf-8 -*-
"""
This module provides a class :class:`Expr` to subclass from in order to
describe expressions. The value of the class is that it allows straightforward
interoperability between python packages handling symbolics (SymPy) and units
(quantities) as well as working without either of those. The price one has to
pay to allow for this is a somewhat contrived syntax
"""
from __future__ import (absolute_import, division, print_function)


class Expr(object):
    '''

    Examples
    --------
    >>> class HeatCapacity(Expr):
    ...     parameter_keys = ('temperature',)
    ...     kw = {'substance': None}
    ...
    >>> import math
    >>> class EinsteinSolid(HeatCapacity):
    ...     """ arguments: einstein temperature """
    ...     def __call__(self, variables, args=None, backend=math):
    ...         molar_mass = self.substance.mass
    ...         TE = self.arg(variables, args, 0)  # einstein_temperature
    ...         R = variables['R']
    ...         T = variables['temperature']
    ...         molar_c_v = 3*R*(TE/(2*T))**2 * backend.sinh(TE/(2*T))**-2
    ...         return molar_c_v/molar_mass
    ...
    >>> from chempy import Substance
    >>> Al = Substance.from_formula('Al', other_properties={'DebyeT': 428})
    >>> Be = Substance.from_formula('Be', other_properties={'DebyeT': 1440})
    >>> einT = lambda s: 0.806*s.other_properties['DebyeT']
    >>> cv = {s.name: EinsteinSolid([einT(s)], substance=s) for s in (Al, Be)}
    >>> print('%.4f' % cv['Al']({'temperature': 273.15, 'R': 8.3145}))  # J/(g*K)
    0.8108
    >>> import sympy
    >>> print(cv['Be']({'temperature': sympy.Symbol('T'), 'R': sympy.Symbol('R')}, backend=sympy))
    112105.346283965*R/(T**2*sinh(580.32/T)**2)

    '''

    parameter_keys = ()
    kw = None

    def __init__(self, args, arg_keys=None, **kwargs):
        self.args = args
        self.arg_keys = arg_keys
        for k, v in (self.kw or {}).items():
            setattr(self, k, kwargs.pop(k, v))
        if kwargs:
            raise ValueError("Unexpected keyword arguments %s" % kwargs)

    def __call__(self, variables, args=None, backend=None):
        raise NotImplementedError("Subclass and implement __call__")

    def arg(self, variables, args, index):
        args = args or self.args
        if self.arg_keys is None:
            return args[index]
        else:
            return variables.get(self.arg_keys[index], args[index])

    def all_args(self, variables, args):
        return [self.arg(variables, args, i) for i in range(len(self.args))]

    def _dedimensionalisation(self, unit_registry):
        from ..units import default_unit_in_registry, to_unitless
        units = [default_unit_in_registry(arg, unit_registry) for arg in self.args]
        unitless_args = [to_unitless(arg, unit) for arg, unit in zip(self.args, units)]
        return units, self.__class__(unitless_args, self.arg_keys, **{k: getattr(self, k) for k in self.kw})
