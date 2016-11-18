# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from functools import reduce
import math
from operator import add

from ._expr import Expr


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


def _mk_Poly(parameter_name, reciprocal=False, shift_name='shift'):
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
    >>> P = _mk_Poly('x')
    >>> p = P([3, 5, 7, 2])
    >>> p.eval_poly({'x': 13}) == 5 + 7*(13-3) + 2*(13-3)**2
    True

    """
    class Poly(Expr):
        """ Args: shift, p0, p1, ... """
        argument_names = (shift_name, Ellipsis)
        parameter_keys = (parameter_name,)
        skip_poly = 0

        def eval_poly(self, variables, backend=math):
            all_args = self.all_args(variables, backend=backend)
            x = variables[parameter_name]
            offset, coeffs = all_args[self.skip_poly], all_args[self.skip_poly+1:]
            return _eval_poly(x, offset, coeffs, reciprocal)
    return Poly


def _mk_PiecewisePoly(parameter, reciprocal=False):
    """ Class factory of Expr subclass for piecewise (shifted) polynomial """
    class PiecewisePoly(Expr):
        """ Args: npolys, ncoeff0, lower0, upper0, ncoeff1, ..., shift0, p0_0, p0_1, ... shiftn, p0_n, p1_n, ... """
        argument_names = ('npolys', Ellipsis)
        parameter_keys = (parameter,)
        skip_poly = 0

        def eval_poly(self, variables, backend=math):
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
