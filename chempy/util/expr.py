# -*- coding: utf-8 -*-
"""
This module provides a class :class:`Expr` to subclass from in order to
describe expressions. The value of the class is that it allows straightforward
interoperability between python packages handling symbolics (SymPy) and units
(quantities) as well as working without either of those. The price one has to
pay to allow for this is a somewhat contrived syntax
"""
from __future__ import (absolute_import, division, print_function)

from itertools import chain


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
    ...     """ arguments: einstein_temperature """
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
    nargs = None

    def __init__(self, args, arg_keys=None, **kwargs):
        if self.nargs is not None and len(args) != self.nargs:
            raise ValueError("Incorrect number of arguments: %d (expected %d)" % (len(args), self.nargs))
        if arg_keys is not None and self.nargs is not None and len(arg_keys) != self.nargs:
            raise ValueError("Incorrect number of arg_keys: %d (expected %d)" % (len(arg_keys), self.nargs))

        self.args = args
        self.arg_keys = arg_keys
        for k, v in (self.kw or {}).items():
            setattr(self, k, kwargs.pop(k, v))
        if kwargs:
            raise ValueError("Unexpected keyword arguments %s" % kwargs)

    def __call__(self, variables, args=None, backend=None):
        raise NotImplementedError("Subclass and implement __call__")

    def _str(self, arg_fmt, arg_keys_fmt=str, with_kw=False):
        if len(self.args) == 0:
            args_str = ''
        elif len(self.args) == 1:
            args_str = '%s,' % self.args[0]
        else:
            args_str = '%s' % ', '.join(map(arg_fmt, self.args))
        args_kwargs_strs = [', '.join(chain(
            ['(%s)' % args_str],
            [arg_keys_fmt(self.arg_keys)] if self.arg_keys is not None else []
        ))]
        print_kw = {k: getattr(self, k) for k in self.kw if getattr(self, k) != self.kw[k]}
        if with_kw and print_kw:
            args_kwargs_strs += [', '.join('{}={}'.format(k, v) for k, v in print_kw.items())]
        return "{}({})".format(self.__class__.__name__, ', '.join(args_kwargs_strs))

    def __repr__(self):
        return self._str(repr)

    def string(self, arg_fmt=str):
        return self._str(arg_fmt)

    def arg(self, variables, args, index):
        args = args or self.args
        if self.arg_keys is None:
            return args[index]
        else:
            return variables.get(self.arg_keys[index], args[index])

    def all_args(self, variables, args):
        return [self.arg(variables, args, i) for i in range(len(args or self.args))]

    def _dedimensionalisation(self, unit_registry):
        from ..units import default_unit_in_registry, to_unitless
        units = [default_unit_in_registry(arg, unit_registry) for arg in self.args]
        unitless_args = [to_unitless(arg, unit) for arg, unit in zip(self.args, units)]
        return units, self.__class__(unitless_args, self.arg_keys, **{k: getattr(self, k) for k in self.kw})


def mk_Poly(parameter, reciprocal=False):
    """ Classfactory of Expr subclass for (shifted) polynomial

    Parameters
    ----------
    parameter: str
        name of paramter
    reciprocal: bool
        whether the polynomial is in the reciprocal of the parameter

    Returns
    -------
    Expr subclass for a shifted polynomial with the args: offset, p0, p1, ...
    the class has the method "eval_poly" with same signature as __call__


    Examples
    --------
    >>> P = mk_Poly('x')
    >>> p = P([3, 5, 7, 2])
    >>> p.eval_poly({'x': 13}) == 5 + 7*(13-3) + 2*(13-3)**2
    True

    """
    class Poly(Expr):

        parameter_keys = (parameter,)

        def eval_poly(self, variables, args=None, backend=None):
            all_args = self.all_args(variables, args)
            offset, coeffs = all_args[0], all_args[1:]
            _x0 = variables[parameter] - offset
            _x = _x0/_x0
            k = None
            for coeff in coeffs:
                if k is None:
                    k = coeff*_x
                else:
                    k += coeff*_x

                if reciprocal is True:
                    _x /= _x0
                elif reciprocal is False:
                    _x *= _x0
                else:
                    raise NotImplementedError
            return k
    return Poly


class PiecewisePoly(Expr):

    def __init__(self, args, arg_keys=None, **kwargs):
        bounds, polys = args
        if not all(p.parameter_keys == polys[0].parameter_keys for p in polys):
            raise ValueError("Mixed parameter_keys")
        super(PiecewisePoly, self).__init__(args, arg_keys, **kwargs)
        self.parameter_keys = polys[0].parameter_keys
        self.bounds = bounds
        self.polys = polys

    def eval_poly(self, variables, args=None, backend=None):
        bounds, polys = self.all_args(variables, args)
        x = variables[self.parameter_keys[0]]
        if args is not None:
            raise NotImplementedError("args not flattened")  # nargs == None
        try:
            pw = backend.Piecewise
        except AttributeError:
            for (lower, upper), p in zip(bounds, polys):
                if lower <= x <= upper:
                    return p.eval_poly(variables, args, backend)
            else:
                raise ValueError("not within any bounds: %s" % str(x))
        else:
            return pw(*[(p.eval_poly(variables, args, backend), backend.And(l <= x, x <= u)) for (l, u), p
                        in zip(bounds, polys)])
