import numpy as np
from scipy.optimize import fsolve

from ..equilibria import (
    equilibrium_quotient, equilibrium_residual, get_rc_interval,
    solve_equilibrium, EqSystem, prodexp
)


def test_equilibrium_quotient():
    assert abs(equilibrium_quotient([2.3, 3.7, 5.1], (-1, -1, 1)) -
               5.1/2.3/3.7) < 1e-14


def test_equilibrium_residual():
    c0 = np.array((13., 11, 17))
    stoich = np.array((-2, 3, -4))
    K = 3.14
    assert abs(equilibrium_residual(0.1, c0, stoich, K) -
               (K - (13-0.2)**-2*(11 + 0.3)**3*(17 - 0.4)**-4)) < 1e-14


def test_get_rc_interval():
    c = np.array((13., 11, 17))
    stoich = np.array((-2, 3, -4))
    limits = get_rc_interval(stoich, c)
    lower = -11/3.
    upper = 17./4
    assert abs(limits[0] - lower) < 1e-14
    assert abs(limits[1] - upper) < 1e-14


def test_solve_equilibrium_1():
    c = np.array((13., 11, 17))
    stoich = np.array((-2, 3, -4))
    K = 3.14

    def f(x):
        return (13 - 2*x)**-2 * (
            11 + 3*x)**3 * (17 - 4*x)**-4 - K
    assert np.allclose(solve_equilibrium(c, stoich, K),
                       c + stoich*fsolve(f, 3.48))

def test_solve_equilibrium_2():
    c = np.array([  1.7000e-03,   3.0000e+06,   3.0000e+06,   9.7000e+07,   5.5500e+09])
    stoich = (1, 1, 0, 0, -1)
    K = 1e-14

    def f(x):
        return prodexp(c, stoich) - K
    assert np.allclose(solve_equilibrium(c, stoich, K),
                       c + stoich*fsolve(f, 3.48))


def test_EqSystem():
    pass
