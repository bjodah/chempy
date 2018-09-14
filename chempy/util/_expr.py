# -*- coding: utf-8 -*-
"""
This module provides a class :class:`Expr` to subclass from in order to
describe expressions. The value of the class is that it allows straightforward
interoperability between python packages handling symbolics (SymPy) and units
(quantities) as well as working without either of those. The price one has to
pay to allow for this is a somewhat contrived syntax.

Note that this module is to be considered an implementation detail, and not
something that should be relied upon in external code.
"""
from __future__ import (absolute_import, division, print_function)

import math
from itertools import chain
from operator import add, mul, truediv, sub, pow
from .pyutil import defaultkeydict, deprecated

try:
    import sympy
except ImportError:
    sympy = None


def _implicit_conversion(obj):
    if isinstance(obj, (int, float)):
        return Constant(obj)
    elif isinstance(obj, Expr):
        return obj
    elif isinstance(obj, str):
        return Symbol(obj)

    if sympy is not None:
        if isinstance(obj, sympy.Mul):
            if len(obj.args) != 2:
                raise NotImplementedError("Did you use evaluate=False?")
            return _MulExpr([_implicit_conversion(obj.args[0]), _implicit_conversion(obj.args[1])])
        elif isinstance(obj, sympy.Add):
            if len(obj.args) != 2:
                raise NotImplementedError("Did you use evaluate=False?")
            return _AddExpr([_implicit_conversion(obj.args[0]), _implicit_conversion(obj.args[1])])
        elif isinstance(obj, sympy.Pow):
            return _PowExpr(_implicit_conversion(obj.base), _implicit_conversion(obj.exp))
        elif isinstance(obj, sympy.Float):
            return Constant(float(obj))
        elif isinstance(obj, sympy.Symbol):
            return Symbol(obj.name)

    raise NotImplementedError(
        "Don't know how to convert %s (of type %s)" % (obj, type(obj)))


class Expr(object):
    ''' Baseclass for Expressions corresponding to physical quantitites.

    The design assumes that a large group of different Expr subclasses may
    be evaluated with some shared state (parameter_keys). The backend kwarg
    in call enables use of e.g. math, numpy or sympy interchangeably.

    Parameters
    ----------
    args : tuple/list of scalars or dict mapping name to scalar
        When dict: it is converted to a list using ``self.argument_names`` or
        ``self.unique_keys``.
    unique_keys : iterable of strings
        Unique names (among all instances) for late overriding, aligned with beginning of
        ``args``.

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
    argument_defaults : tuple of floats, optional
        Default values for arguments, aligned from the end of argument names.
    parameter_keys : tuple of strings
    nargs : int
        number of arguments (`None` signifies unset, -1 signifies any number)
    '''

    argument_names = None
    argument_defaults = None
    parameter_keys = ()
    nargs = None

    def __init__(self, args=None, unique_keys=None):
        if isinstance(args, str):
            args = (args,)
        if self.argument_names is not None and self.argument_names[-1] != Ellipsis and self.nargs is None:
            self.nargs = len(self.argument_names)
        if self.argument_defaults is not None:
            if self.nargs == -1:
                raise ValueError("Cannot have defaults when number of arguments is unbounded.")
            if len(self.argument_defaults) > len(self.argument_names):
                raise ValueError("Cannot have more defaults than actual arguments")
            if args is not None:
                n_missing = self.nargs - len(args)
                if n_missing > 0:
                    args = tuple(chain(args, self.argument_defaults[-n_missing:]))

        if self.nargs == 1 and (isinstance(args, (float, int)) or
                                getattr(args, 'ndim', -1) == 0 or
                                isinstance(args, Expr)):
            args = [args]
            nargs = 1
        elif args is None:
            nargs = None
        else:
            nargs = len(args)

        if self.nargs not in (None, -1) and nargs is not None and nargs != self.nargs:
            raise ValueError("Incorrect number of arguments: %d (expected %d)" % (nargs, self.nargs))
        if unique_keys is not None and self.nargs is not None and len(unique_keys) > self.nargs:
            raise ValueError("Incorrect number of unique_keys: %d (expected %d or less)" % (
                len(unique_keys), self.nargs))
        self.unique_keys = None if unique_keys is None else tuple(unique_keys)

        if isinstance(args, dict):
            args = [args[k] for k in self.argument_names or self.unique_keys]

        self.args = args

    @classmethod
    def fk(cls, *args):
        """ Alternative constructor "from keys", \\*args is used as ``unique_keys``. """
        return cls(unique_keys=args)

    @classmethod
    def from_callback(cls, callback, attr='__call__', **kwargs):
        """ Factory of subclasses

        Parameters
        ----------
        callback : callable
            signature: *args, backend=math
        attr : str
            What attribute to override
        argument_names : tuple of str, optional
        argument_defaults : tuple of floats, optional
        parameter_keys : tuple of str, optional,
        nargs : int, optional

        Examples
        --------
        >>> from operator import add; from functools import reduce
        >>> def poly(args, x, backend=math):
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
        def body(self, variables, backend=math, **kw):
            args = self.all_args(variables, backend=backend)
            params = self.all_params(variables, backend=backend)
            return callback(args, *params, backend=backend, **kw)

        class Wrapper(cls):
            pass

        setattr(Wrapper, attr, body)
        Wrapper.__name__ = callback.__name__
        for k, v in kwargs.items():
            setattr(Wrapper, k, v)
        return Wrapper

    def __call__(self, variables, backend=math, **kwargs):
        raise NotImplementedError("Subclass and implement __call__")

    def _all_keys(self, attr):
        _keys = getattr(self, attr)
        _all = set() if _keys is None else set(_keys)
        if self.args is not None:
            for arg in self.args:
                if isinstance(arg, Expr):
                    _all = _all.union(arg._all_keys(attr))
        return _all

    def all_parameter_keys(self):
        return self._all_keys('parameter_keys')

    def all_unique_keys(self):
        return self._all_keys('unique_keys')

    def _str(self, arg_fmt, unique_keys_fmt=str):
        if self.args is None or len(self.args) == 0:
            args_str = ''
        elif len(self.args) == 1:
            args_str = '%s,' % self.args[0]
        else:
            args_str = '%s' % ', '.join(map(arg_fmt, self.args))
        args_strs = [', '.join(chain(
            ['(%s)' % args_str],
            [unique_keys_fmt(self.unique_keys)] if self.unique_keys is not None else []
        ))]
        return "{}({})".format(self.__class__.__name__, ', '.join(args_strs))

    def __repr__(self):
        return self._str(repr)

    def string(self, arg_fmt=str, **kwargs):
        return self._str(arg_fmt, **kwargs)

    def arg(self, variables, index, backend=math, evaluate=True, **kwargs):
        """
        Parameters
        ----------
        variables : container
        index : int or str
            When str: index from ``self.argument_names``.
        backend : module
        evaluate : bool

        Notes
        -----
        Priority:
            1. unique_keys
            2. variables[k] for k in argument_names
        """
        if isinstance(index, str):
            index = self.argument_names.index(index)

        if self.unique_keys is None:
            res = self.args[index]
        elif index < len(self.unique_keys):
            uk = self.unique_keys[index]
            try:
                res = variables[uk]
            except KeyError:
                if self.args is None:
                    raise KeyError("Unique key missing: %s" % uk)
                else:
                    res = self.args[index]
        else:
            if self.args is None or index > len(self.args):
                res = self.argument_defaults[index - self.nargs + len(self.argument_defaults)]
            else:
                res = self.args[index]

        if isinstance(res, str):
            res = variables[res]
        elif isinstance(res, Symbol):
            res = variables[res.args[0]]

        if isinstance(res, Expr) and evaluate:
            return res(variables, backend=backend, **kwargs)
        else:
            return res

    def all_args(self, variables, backend=math, evaluate=True, **kwargs):
        if self.nargs is None or self.nargs == -1:
            nargs = len(self.args)
        else:
            nargs = self.nargs
        return [self.arg(variables, i, backend, evaluate, **kwargs) for i in range(nargs)]

    def all_params(self, variables, backend=math):
        return [v(variables, backend=backend) if isinstance(v, Expr) else v for v
                in [variables[k] for k in self.parameter_keys]]

    def args_dimensionality(self, **kwargs):
        """ return tuple of dicts mapping str to int ('length', 'mass', 'time', 'current',
        'temperature', 'luminous_intensity', 'amount') """
        raise NotImplementedError("method not implemented in subclass.")

    def dedimensionalisation(self, unit_registry, variables={}, backend=math):
        """ Create an instance with consistent units from a unit_registry

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
        ...     def __call__(self, variables, backend=math, **kwargs):
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
        new_units: list of units of the dedimensionalised args.
        self.__class__ instance: with dedimensioanlised arguments

        """
        from ..units import default_unit_in_registry, to_unitless, unitless_in_registry
        new_units = []
        if self.args is None:
            unitless_args = None
        else:
            unitless_args = []
            units = [None if isinstance(arg, Expr) else default_unit_in_registry(arg, unit_registry) for arg
                     in self.all_args(variables, backend=backend, evaluate=False)]
            for arg, unit in zip(self.all_args(variables, backend=backend, evaluate=False), units):
                if isinstance(arg, Expr):
                    if unit is not None:
                        raise ValueError()
                    _unit, _dedim = arg.dedimensionalisation(unit_registry, variables, backend=backend)
                else:
                    _unit, _dedim = unit, to_unitless(arg, unit)
                new_units.append(_unit)
                unitless_args.append(_dedim)
        instance = self.__class__(unitless_args, self.unique_keys)
        if self.argument_defaults is not None:
            instance.argument_defaults = tuple(unitless_in_registry(arg, unit_registry)
                                               for arg in self.argument_defaults)

        return new_units, instance

    def _sympy_format(self, method, variables, backend, default, **kwargs):
        variables = variables or {}
        if backend in (None, math):
            backend = sympy
        variables = defaultkeydict(
            None if default is None else (lambda k: backend.Symbol(default(k))),
            {k: v if isinstance(v, Expr) else (backend.Symbol(v) if isinstance(v, str) else backend.Float(v))
             for k, v in variables.items()})
        expr = self(variables, backend=backend, **kwargs).simplify()
        if method == 'latex':
            return backend.latex(expr)
        elif method == 'str':
            return str(expr)
        elif method == 'unicode':
            return backend.pretty(expr, use_unicode=True)
        elif method == 'mathml':
            from sympy.printing.mathml import mathml
            return mathml(expr)
        else:
            raise NotImplementedError("Unknown method: %s" % method)

    def latex(self, variables=None, backend=math, default=None):
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
        >>> p.latex({'R': 'R', 'temp': 'T', 'vol': 'V'})  # doctest: +SKIP
        '\\frac{7 R T}{V}'

        Notes
        -----
        Requires SymPy

        """
        return self._sympy_format('latex', variables, backend=backend, default=default)

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        if self.args is None and other.args is None:
            for uk1, uk2 in zip(self.unique_keys, other.unique_keys):
                if uk1 != uk2:
                    return False
            return True
        if self.args is None or other.args is None:
            return False
        if len(self.args) != len(other.args):
            return False
        from ..units import compare_equality
        for arg0, arg1 in zip(self.args, other.args):
            if not compare_equality(arg0, arg1):
                return False
        return True

    def __add__(self, other):
        if other == other*0:
            return self
        return _AddExpr([self, _implicit_conversion(other)])

    def __sub__(self, other):
        if other == other*0:
            return self
        return _SubExpr([self, _implicit_conversion(other)])

    def __mul__(self, other):
        if other == 1:
            return self
        return _MulExpr([self, _implicit_conversion(other)])

    def __truediv__(self, other):
        if other == 1:
            return self
        return _DivExpr([self, _implicit_conversion(other)])

    def __rtruediv(self, other):
        return _DivExpr([_implicit_conversion(other), self])

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
        return _DivExpr([_implicit_conversion(other), self])

    def __pow__(self, other):
        return _PowExpr([self, _implicit_conversion(other)])

    def __rpow__(self, other):
        return _PowExpr([_implicit_conversion(other), self])


class _NegExpr(Expr):

    def __call__(self, variables, backend=math, **kwargs):
        arg0, = self.all_args(variables, backend=backend, **kwargs)
        return -arg0

    def rate_coeff(self, *args, **kwargs):
        return -self.args[0].rate_coeff(*args, **kwargs),


class _BinaryExpr(Expr):
    _op = None

    def _str(self, *args, **kwargs):
        return ("({0} %s {1})" % self._op_str).format(*[arg._str(*args, **kwargs) for arg in self.args])

    def __call__(self, variables, backend=math, **kwargs):
        arg0, arg1 = self.all_args(variables, backend=backend, **kwargs)
        return self._op(arg0, arg1)

    def rate_coeff(self, *args, **kwargs):
        return self._op(self.args[0].rate_coeff(*args, **kwargs),
                        self.args[1].rate_coeff(*args, **kwargs))


class _AddExpr(_BinaryExpr):
    _op = add
    _op_str = '+'


class _SubExpr(_BinaryExpr):
    _op = sub
    _op_str = '-'


class _MulExpr(_BinaryExpr):
    _op = mul
    _op_str = '*'


class _DivExpr(_BinaryExpr):
    _op = truediv
    _op_str = '/'


class _PowExpr(_BinaryExpr):
    _op = pow
    _op_str = '**'


class Constant(Expr):
    nargs = 1

    def __call__(self, variables, backend=None, **kwargs):
        return self.args[0]

    def rate_coeff(self, *args, **kwargs):
        return self.args[0]


class Symbol(Expr):
    nargs = 1

    def __call__(self, variables, backend=None, **kwargs):
        return variables[self.args[0]]


class Function(Expr):
    pass


class UnaryFunction(Function):
    nargs = 1
    _func_name = None

    def __call__(self, variables, backend=math, **kwargs):
        arg, = self.all_args(variables, backend=backend, **kwargs)
        return getattr(backend, self._func_name)(arg)

    def rate_coeff(self, *args, **kwargs):
        return getattr(kwargs.get('backend', math), self._func_name)(self.args[0].rate_coeff(*args, **kwargs))


class BinaryFunction(Function):
    nargs = 2
    _func_name = None


class Log10(UnaryFunction):
    _func_name = 'log10'


class Exp(UnaryFunction):
    _func_name = 'exp'


def create_Piecewise(parameter_name, nan_fallback=False):
    """
    Examples
    --------
    >>> Power = Expr.from_callback(lambda args, x, backend=None: args[0]*x**args[1],
    ...     argument_names=('scale', 'pow'), parameter_keys=('x',))
    >>> minus_x = Power([-1, 1])
    >>> cube = Power([1, 3])
    >>> PW = create_Piecewise('x')
    >>> pw = PW([-float('inf'), minus_x, 0, cube, float('inf')])
    >>> pw({'x': -5}) == 5
    True
    >>> pw({'x': 2}) == 8
    True

    """
    def _pw(bounds_exprs, x, backend=math, **kwargs):
        if len(bounds_exprs) < 3:
            raise ValueError("Need at least 3 args")
        if len(bounds_exprs) % 2 != 1:
            raise ValueError("Need an odd number of bounds/exprs")
        n_exprs = (len(bounds_exprs) - 1) // 2
        lower = [bounds_exprs[2*(i+0)] for i in range(n_exprs)]
        upper = [bounds_exprs[2*(i+1)] for i in range(n_exprs)]
        exprs = [bounds_exprs[2*i + 1] for i in range(n_exprs)]

        try:
            pw = backend.Piecewise
        except AttributeError:
            for lo, up, ex in zip(lower, upper, exprs):
                if lo <= x <= up:
                    return ex
            else:
                raise ValueError("not within any bounds: %s" % x)
        else:
            _NAN = backend.Symbol('NAN')
            return pw(*([(ex, backend.And(lo <= x, x <= up)) for lo, up, ex in zip(lower, upper, exprs)] +
                        ([(_NAN, True)] if nan_fallback else [])))

    return Expr.from_callback(_pw, parameter_keys=(parameter_name,))


def create_Poly(parameter_name, reciprocal=False, shift=None, name=None):
    """
    Examples
    --------
    >>> Poly = create_Poly('x')
    >>> p1 = Poly([3, 4, 5])
    >>> p1({'x': 7}) == 3 + 4*7 + 5*49
    True
    >>> RPoly = create_Poly('T', reciprocal=True)
    >>> p2 = RPoly([64, 32, 16, 8])
    >>> p2({'T': 2}) == 64 + 16 + 4 + 1
    True
    >>> SPoly = create_Poly('z', shift=True)
    >>> p3 = SPoly([7, 2, 3, 5], unique_keys=('z0',))
    >>> p3({'z': 9}) == 2 + 3*(9-7) + 5*(9-7)**2
    True
    >>> p3({'z': 9, 'z0': 6}) == 2 + 3*(9-6) + 5*(9-6)**2
    True

    """
    if shift is True:
        shift = 'shift'

    def _poly(args, x, backend=math, **kwargs):

        if shift is None:
            coeffs = args
            x0 = x
        else:
            coeffs = args[1:]
            x_shift = args[0]
            x0 = x - x_shift

        cur = 1
        res = None
        for coeff in coeffs:
            if res is None:
                res = coeff*cur
            else:
                res += coeff*cur

            if reciprocal:
                cur /= x0
            else:
                cur *= x0
        return res

    if shift is None:
        argument_names = None
    else:
        argument_names = (shift, Ellipsis)
    if name is not None:
        _poly.__name__ = name
    return Expr.from_callback(_poly, parameter_keys=(parameter_name,), argument_names=argument_names)

from ._expr_deprecated import _mk_PiecewisePoly, _mk_Poly  # noqa

mk_PiecewisePoly = deprecated(use_instead=create_Piecewise)(_mk_PiecewisePoly)
mk_Poly = deprecated(use_instead=create_Poly)(_mk_Poly)
