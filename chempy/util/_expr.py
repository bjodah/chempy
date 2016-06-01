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
from operator import add, mul, truediv, sub
from .pyutil import defaultkeydict


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
        note that parameters in :attr:`kw` are not processed in e.g. dedimensionalisation.

    Examples
    --------
    >>> class HeatCapacity(Expr):
    ...     parameter_keys = ('temperature',)
    ...
    >>> import math
    >>> class EinsteinSolid(HeatCapacity):
    ...     parameter_keys = HeatCapacity.parameter_keys + ('molar_gas_constant',)
    ...     argument_names = ('einstein_temperature', 'molar_mass')
    ...
    ...     def __call__(self, variables, backend=math):
    ...         TE, molar_mass = self.all_args(variables, backend=backend)  # einstein_temperature
    ...         T, R = self.all_params(variables, backend=backend)
    ...         # Canonical ensemble:
    ...         molar_c_v = 3*R*(TE/(2*T))**2 * backend.sinh(TE/(2*T))**-2
    ...         return molar_c_v/molar_mass
    ...
    >>> from chempy import Substance
    >>> Al = Substance.from_formula('Al', data={'DebyeT': 428})
    >>> Be = Substance.from_formula('Be', data={'DebyeT': 1440})
    >>> einT = lambda s: 0.806*s.data['DebyeT']
    >>> cv = {s.name: EinsteinSolid([einT(s), s.mass]) for s in (Al, Be)}
    >>> print('%.4f' % cv['Al']({'temperature': 273.15, 'molar_gas_constant': 8.3145}))  # J/(g*K)
    0.8108
    >>> import sympy; from sympy import Symbol as Symb
    >>> print(cv['Be']({'temperature': Symb('T'), 'molar_gas_constant': Symb('R')}, backend=sympy))
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
    print_name : str
        Name of class

    '''

    argument_names = None
    parameter_keys = ()
    kw = None
    nargs = None
    print_name = None

    def __init__(self, args, unique_keys=None, **kwargs):
        if self.argument_names is not None and self.argument_names[-1] != Ellipsis and self.nargs is None:
            self.nargs = len(self.argument_names)
        if self.nargs == 1 and (isinstance(args, (float, int)) or getattr(args, 'ndim', -1) == 0):
            args = [args]
            nargs = 1
        else:
            nargs = len(args)
        if self.nargs not in (None, -1) and nargs != self.nargs:
            raise ValueError("Incorrect number of arguments: %d (expected %d)" % (nargs, self.nargs))
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

    @classmethod
    def from_callback(cls, callback, **kwargs):
        """ Factory of subclasses

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
        >>> Poly = Expr.from_callback(poly, parameter_keys=('x',), argument_names=('x0', Ellipsis))
        >>> p = Poly([1, 3, 2, 5])
        >>> p({'x': 7}) == 3 + 2*(7-1) + 5*(7-1)**2
        True
        >>> q = Poly([1, 3, 2, 5], unique_keys=('x0_q',))
        >>> q({'x': 7, 'x0_q': 0}) == 3 + 2*7 + 5*7**2
        True

        """
        class Wrapper(cls):
            def __call__(self, variables, backend=math):
                args = self.all_args(variables, backend=backend)
                params = self.all_params(variables, backend=backend)
                return callback(args, *params, backend=backend)
        for k, v in kwargs.items():
            setattr(Wrapper, k, v)
        return Wrapper

    def __call__(self, variables, backend=None):
        raise NotImplementedError("Subclass and implement __call__")

    @property
    def kwargs(self):
        return {k: getattr(self, k, v) for k, v in self.kw.items()}

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
        if self.kw is not None:
            print_kw = {k: getattr(self, k) for k in self.kw if getattr(self, k) != self.kw[k]}
        if with_kw and self.kw is not None and print_kw:
            args_kwargs_strs += [', '.join('{}={}'.format(k, v) for k, v in print_kw.items())]
        return "{}({})".format(self.print_name or self.__class__.__name__, ', '.join(args_kwargs_strs))

    def __repr__(self):
        return self._str(repr)

    def string(self, arg_fmt=str, **kwargs):
        return self._str(arg_fmt, **kwargs)

    def arg(self, variables, index, backend=None, evaluate=True):
        if isinstance(index, str):
            index = self.argument_names.index(index)
        if self.unique_keys is None or len(self.unique_keys) <= index:
            res = self.args[index]
        else:
            res = variables.get(self.unique_keys[index], self.args[index])
        if isinstance(res, Expr) and evaluate:
            return res(variables, backend=backend)
        else:
            return res

    def all_args(self, variables, backend=None, evaluate=True):
        return [self.arg(variables, i, backend, evaluate) for i in range(len(self.args))]

    def all_params(self, variables, backend=None):
        return [v(variables, backend=backend) if isinstance(v, Expr) else v for v
                in [variables[k] for k in self.parameter_keys]]

    def dedimensionalisation(self, unit_registry, variables={}, backend=math):
        """ Create an instance with consistent units

        Parameters
        ----------
        unit_registry : dict
        variables : dict
        backend : module

        Examples
        --------
        >>> class Pressure(Expr):
        ...     argument_names = ('n',)
        ...     parameter_keys = ('temperature', 'volume', 'R')
        ...     def __call__(self, variables, backend=None):
        ...         n, = self.all_args(variables, backend=backend)
        ...         T, V, R = self.all_params(variables, backend=backend)
        ...         return n*R*T/V
        ...
        >>> from chempy.units import SI_base_registry, default_units as u
        >>> p = Pressure([2*u.micromole])
        >>> units, d = p.dedimensionalisation(SI_base_registry)
        >>> units[0] == 1e6*u.micromole
        True
        >>> d.args[0] == 2e-6
        True


        Returns
        -------
        A new instance where all args have been (recursively) expressed in the unit system
        of unit_registry

        """
        from ..units import default_unit_in_registry, to_unitless
        units = [None if isinstance(arg, Expr) else default_unit_in_registry(arg, unit_registry) for arg
                 in self.all_args(variables, backend=backend, evaluate=False)]
        new_units, unitless_args = [], []
        for arg, unit in zip(self.all_args(variables, backend=backend, evaluate=False), units):
            if isinstance(arg, Expr):
                if unit is not None:
                    raise ValueError()
                _unit, _dedim = arg.dedimensionalisation(unit_registry, variables, backend=backend)
            else:
                _unit, _dedim = unit, to_unitless(arg, unit)
            new_units.append(_unit)
            unitless_args.append(_dedim)
        if self.kw is None:
            kw = {}
        else:
            kw = {k: getattr(self, k) for k in self.kw}
        return new_units, self.__class__(unitless_args, self.unique_keys, **kw)

    def _sympy_format(self, method, variables, backend, default):
        variables = variables or {}
        if backend is None:
            import sympy as backend
        variables = defaultkeydict(
            None if default is None else (lambda k: backend.Symbol(default(k))),
            {k: v if isinstance(v, Expr) else backend.Symbol(v) for k, v in variables.items()})
        expr = self(variables, backend=backend).simplify()
        if method == 'latex':
            return backend.latex(expr)
        elif method == 'unicode':
            return backend.pprint(expr, use_unicode=True)
        elif method == 'mathml':
            from sympy.printing.mathml import print_mathml
            return print_mathml(expr)
        else:
            raise NotImplementedError("Unknown method: %s" % method)

    def latex(self, variables=None, backend=None, default=None):
        r"""
        Parameters
        ----------
        variables : dict
        backend : module
        default : callable
            Format string based on missing key, signature: str -> str.

        Examples
        --------
        >>> def pressure(args, *params, **kw):
        ...     return args[0]*params[0]*params[1]/params[2]
        >>> Pressure = Expr.from_callback(pressure, parameter_keys='R temp vol'.split(), nargs=1)
        >>> p = Pressure([7])
        >>> p.latex({'R': 'R', 'temp': 'T', 'vol': 'V'})
        '\\frac{7 R}{V} T'

        Notes
        -----
        Requires SymPy

        """
        return self._sympy_format('latex', variables, backend=backend, default=default)

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        if len(self.args) != len(other.args):
            return False
        for arg0, arg1 in zip(self.args, other.args):
            if arg0 != arg1:
                return False
        if self.kw is not None:
            for k in self.kw:
                print(k)
                print(getattr(self, k), getattr(other, k))
                if getattr(self, k) != getattr(other, k):
                    return False
        return True

    def __add__(self, other):
        if other == other*0:
            return self
        return _AddExpr([self, other])

    def __sub__(self, other):
        if other == other*0:
            return self
        return _SubExpr([self, other])

    def __mul__(self, other):
        if other == 1:
            return self
        return _MulExpr([self, other])

    def __truediv__(self, other):
        if other == 1:
            return self
        return _DivExpr([self, other])

    def __neg__(self):
        if isinstance(self, _NegExpr):
            return self.args[0]
        return _NegExpr((self,))

    def __radd__(self, other):
        return self+other

    def __rmul__(self, other):
        return self*other

    def __rsub__(self, other):
        return (-self) + other

    def __rtruediv__(self, other):
        return _DivExpr([other, self])


class _NegExpr(Expr):

    def __call__(self, variables, backend=None):
        arg0, = self.all_args(variables, backend=backend)
        return -arg0


class _BinaryExpr(Expr):
    _op = None

    def __call__(self, variables, backend=None):
        arg0, arg1 = self.all_args(variables, backend=backend)
        return self._op(arg0, arg1)


class _AddExpr(_BinaryExpr):
    _op = add


class _SubExpr(_BinaryExpr):
    _op = sub


class _MulExpr(_BinaryExpr):
    _op = mul


class _DivExpr(_BinaryExpr):
    _op = truediv


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
            all_args = self.all_args(variables, backend=backend)
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
            all_args = self.all_args(variables, backend=backend)[self.skip_poly:]
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
