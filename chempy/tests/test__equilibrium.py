# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from ..util.testing import requires

try:
    import numpy as np
except ImportError:
    np, NumSysLin, NumSysLog = [None]*3
else:
    from scipy.optimize import fsolve

from .._equilibrium import equilibrium_residual, solve_equilibrium, _get_rc_interval
from .._util import prodpow


@requires('numpy')
def test_equilibrium_residual():
    c0 = np.array((13., 11, 17))
    stoich = np.array((-2, 3, -4))
    K = 3.14
    assert abs(equilibrium_residual(0.1, c0, stoich, K) -
               (K - (13-0.2)**-2*(11 + 0.3)**3*(17 - 0.4)**-4)) < 1e-14


@requires('numpy')
def test_get__rc_interval():
    c = np.array((13., 11, 17))
    stoich = np.array((-2, 3, -4))
    limits = _get_rc_interval(stoich, c)
    lower = -11/3.
    upper = 17./4
    assert abs(limits[0] - lower) < 1e-14
    assert abs(limits[1] - upper) < 1e-14


@requires('numpy')
def test_solve_equilibrium_1():
    c = np.array((13., 11, 17))
    stoich = np.array((-2, 3, -4))
    K = 3.14

    def f(x):
        return (13 - 2*x)**-2 * (
            11 + 3*x)**3 * (17 - 4*x)**-4 - K
    assert np.allclose(solve_equilibrium(c, stoich, K),
                       c + stoich*fsolve(f, 3.48))


@requires('numpy')
def test_solve_equilibrium_2():
    c = np.array([1.7e-03, 3.0e+06, 3.0e+06, 9.7e+07, 5.55e+09])
    stoich = (1, 1, 0, 0, -1)
    K = 55*1e-6

    def f(x):
        return prodpow(c+x*stoich, stoich) - K
    solution = solve_equilibrium(c, stoich, K)
    assert np.allclose(solution, c + stoich*fsolve(f, 0.1))
