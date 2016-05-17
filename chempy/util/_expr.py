# -*- coding: utf-8 -*-
"""
This module provides a class :class:`Expr` to subclass from in order to
describe expressions. The value of the class is that it allows straightforward
interoperability between python packages handling symbolics (SymPy) and units
(quantities) as well as working without either of those. The price one has to
pay to allow for this is a somewhat contrived syntax
"""
from __future__ import (absolute_import, division, print_function)

import math
from functools import reduce
from itertools import chain
from operator import add


class Expr(object):
    ''' Baseclass for Expressions corresponding to physical quantitites.

    The design assumes that a large group of different Expr subclasses may
    be evaluated with some shared state (parameter_keys). The backend kwarg
    in call enables use of e.g. math, numpy or sympy interchangeably.

    Parameters
    ----------
    args : tuple/list of scalars or dict mapping name to scalar
        When dict it is converted to a list using self.argument_names or self.unique_keys
    unique_keys : iterable of strings
        Unique names (among all instances) for late overriding
    \*\*kwargs :
        keyword arguments intercepted in subclasses (directed by :attr:`kw`)

    Examples
    --------
    >>> class HeatCapacity(Expr):
    ...     parameter_keys = ('temperature',)
    ...     kw = {'substance': None}
    ...
    >>> import math
    >>> class EinsteinSolid(HeatCapacity):
    ...     argument_names = ('einstein_temperature',)
    ...     def __call__(self, variables, backend=math):
    ...         molar_mass = self.substance.mass
    ...         TE = self.arg(variables, 'einstein_temperature')  # einstein_temperature
    ...         R = variables['R']  # shared "state"
    ...         T = variables['temperature']  # shared state
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

    Attributes
    ----------
    argument_names : tuple of strings, optional
        For documentation and referencing positional arguments in self.args
        If set, and `nargs` is `None`: its length is used to set `nargs`
        (unless argument_names ends with an Ellipsis()).
    parameter_keys : tuple of strings
    kw : dict or None
        kwargs to be intercepted in __init__ and set as attributes
    nargs : int
        number of arguments (`None` signifies unset, -1 signifies any number)

    '''

    argument_names = None
    parameter_keys = ()
    kw = None
    nargs = None

    def __init__(self, args, unique_keys=None, **kwargs):
        if self.argument_names is not None and self.argument_names[-1] != Ellipsis and self.nargs is None:
            self.nargs = len(self.argument_names)
        if self.nargs not in (None, -1) and len(args) != self.nargs:
            raise ValueError("Incorrect number of arguments: %d (expected %d)" % (len(args), self.nargs))
        if unique_keys is not None and self.nargs is not None and len(unique_keys) != self.nargs:
            raise ValueError("Incorrect number of unique_keys: %d (expected %d)" % (len(unique_keys), self.nargs))
        self.unique_keys = unique_keys

        if isinstance(args, dict):
            args = [args[k] for k in self.argument_names or self.unique_keys]

        self.args = args
        for k, v in (self.kw or {}).items():
            setattr(self, k, kwargs.pop(k, v))
        if kwargs:
            raise ValueError("Unexpected keyword arguments %s" % kwargs)

    def __call__(self, variables, backend=None):
        raise NotImplementedError("Subclass and implement __call__")

    def _str(self, arg_fmt, unique_keys_fmt=str, with_kw=False):
        if len(self.args) == 0:
            args_str = ''
        elif len(self.args) == 1:
            args_str = '%s,' % self.args[0]
        else:
            args_str = '%s' % ', '.join(map(arg_fmt, self.args))
        args_kwargs_strs = [', '.join(chain(
            ['(%s)' % args_str],
            [unique_keys_fmt(self.unique_keys)] if self.unique_keys is not None else []
        ))]
        print_kw = {k: getattr(self, k) for k in self.kw if getattr(self, k) != self.kw[k]}
        if with_kw and print_kw:
            args_kwargs_strs += [', '.join('{}={}'.format(k, v) for k, v in print_kw.items())]
        return "{}({})".format(self.__class__.__name__, ', '.join(args_kwargs_strs))

    def __repr__(self):
        return self._str(repr)

    def string(self, arg_fmt=str):
        return self._str(arg_fmt)

    def arg(self, variables, index):
        if isinstance(index, str):
            index = self.argument_names.index(index)
        if self.unique_keys is None or len(self.unique_keys) <= index:
            return self.args[index]
        else:
            return variables.get(self.unique_keys[index], self.args[index])

    def all_args(self, variables):
        return [self.arg(variables, i) for i in range(len(self.args))]

    def _dedimensionalisation(self, unit_registry):
        from ..units import default_unit_in_registry, to_unitless
        units = [default_unit_in_registry(arg, unit_registry) for arg in self.args]
        unitless_args = [to_unitless(arg, unit) for arg, unit in zip(self.args, units)]
        return units, self.__class__(unitless_args, self.unique_keys, **{k: getattr(self, k) for k in self.kw})


def Expr_from_callback(callback, **kwargs):
    """ Factory of Expr subclasses

    Parameters
    ----------
    callback : callable
        signature: *args, backend=None
    argument_names : tuple of str, optional
    parameter_keys : tuple of str, optional,
    kw : dict, optional
    nargs : int, optional

    Examples
    --------
    >>> from operator import add; from functools import reduce
    >>> def poly(args, x, backend=None):
    ...     x0 = args[0]
    ...     return reduce(add, [c*(x-x0)**i for i, c in enumerate(args[1:])])
    ...
    >>> Poly = Expr_from_callback(poly, parameter_keys=('x',), argument_names=('x0', Ellipsis))
    >>> p = Poly([1, 3, 2, 5])
    >>> p({'x': 7}) == 3 + 2*(7-1) + 5*(7-1)**2
    True
    >>> q = Poly([1, 3, 2, 5], unique_keys=('x0_q',))
    >>> q({'x': 7, 'x0_q': 0}) == 3 + 2*7 + 5*7**2
    True

    """
    class Wrapper(Expr):
        def __call__(self, variables, backend=math):
            params = [variables[k] for k in self.parameter_keys]
            return callback(self.all_args(variables), *params, backend=backend)
    for k, v in kwargs.items():
        setattr(Wrapper, k, v)
    return Wrapper


def _eval_poly(x, offset, coeffs, reciprocal=False):
    _x0 = x - offset
    _x = _x0/_x0
    res = None
    for coeff in coeffs:
        if res is None:
            res = coeff*_x
        else:
            res += coeff*_x

        if reciprocal:
            _x /= _x0
        else:
            _x *= _x0
    return res


def mk_Poly(parameter, reciprocal=False):
    """ Class factory of Expr subclass for (shifted) polynomial

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
        """ Args: shift, p0, p1, ... """
        argument_names = ('shift', Ellipsis)
        parameter_keys = (parameter,)
        skip_poly = 0

        def eval_poly(self, variables, backend=None):
            all_args = self.all_args(variables)
            x = variables[parameter]
            offset, coeffs = all_args[self.skip_poly], all_args[self.skip_poly+1:]
            return _eval_poly(x, offset, coeffs, reciprocal)
    return Poly


def mk_PiecewisePoly(parameter, reciprocal=False):
    """ Class factory of Expr subclass for piecewise (shifted) polynomial """
    class PiecewisePoly(Expr):
        """ Args: npolys, ncoeff0, lower0, upper0, ncoeff1, ..., shift0, p0_0, p0_1, ... shiftn, p0_n, p1_n, ... """
        argument_names = ('npolys', Ellipsis)
        parameter_keys = (parameter,)
        skip_poly = 0

        def eval_poly(self, variables, backend=None):
            all_args = self.all_args(variables)[self.skip_poly:]
            npoly = all_args[0]
            arg_idx = 1
            poly_args = []
            meta = []
            for poly_idx in range(npoly):
                meta.append(all_args[arg_idx:arg_idx+3])  # nargs, lower, upper
                arg_idx += 3
            for poly_idx in range(npoly):
                narg = 1+meta[poly_idx][0]
                poly_args.append(all_args[arg_idx:arg_idx+narg])
                arg_idx += narg
            if arg_idx != len(all_args):
                raise Exception("Bug in PiecewisePoly.eval_poly")

            x = variables[parameter]
            try:
                pw = backend.Piecewise
            except AttributeError:
                for (ncoeff, lower, upper), args in zip(meta, poly_args):
                    if lower <= x <= upper:
                        return _eval_poly(x, args[0], args[1:], reciprocal)
                else:
                    raise ValueError("not within any bounds: %s" % str(x))
            else:
                return pw(*[(_eval_poly(x, a[0], a[1:], reciprocal),
                             backend.And(l <= x, x <= u)) for (n, l, u), a in zip(meta, poly_args)])

        @classmethod
        def from_polynomials(cls, bounds, polys, inject=[], **kwargs):
            if any(p.parameter_keys != (parameter,) for p in polys):
                raise ValueError("Mixed parameter_keys")
            npolys = len(polys)
            if len(bounds) != npolys:
                raise ValueError("Length mismatch")

            meta = reduce(add, [[len(p.args[p.skip_poly:]) - 1, l, u] for (l, u), p in zip(bounds, polys)])
            p_args = reduce(add, [p.args[p.skip_poly:] for p in polys])
            return cls(inject + [npolys] + meta + p_args, **kwargs)
    return PiecewisePoly
