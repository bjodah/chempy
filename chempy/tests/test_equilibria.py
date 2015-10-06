import pytest
import numpy as np
from scipy.optimize import fsolve

from ..chemistry import (
    Solute, Substance, Reaction
)

from ..equilibria import (
    equilibrium_quotient, equilibrium_residual, get_rc_interval,
    solve_equilibrium, EqSystemBase, prodpow, _solve_equilibrium_coord,
    Equilibrium, EqSystemLin
)

from .ammonical_cupric_solution import get_ammonical_cupric_eqsys


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
    c = np.array([1.7e-03, 3.0e+06, 3.0e+06, 9.7e+07, 5.55e+09])
    stoich = (1, 1, 0, 0, -1)
    K = 55*1e-6

    def f(x):
        return prodpow(c+x*stoich, stoich) - K
    solution = solve_equilibrium(c, stoich, K)
    assert np.allclose(solution, c + stoich*fsolve(f, 0.1))


def test_EqSystemBase():
    a, b = sbstncs = Substance('a'), Substance('b')
    rxns = [Reaction({a: 1}, {b: 1})]
    es = EqSystemBase(rxns, sbstncs)
    assert es.stoichs.tolist() == [[-1, 1]]


def _get_es1():
    a, b = sbstncs = Solute('a'), Solute('b')
    rxns = [Equilibrium({a: 1}, {b: 1}, 3)]
    return EqSystemBase(rxns, sbstncs)


def _get_es_water(EqSys=EqSystemBase):
    H2O = Solute('H2O', charge=0, composition={1: 2, 8: 1})
    OHm = Solute('OH-', charge=-1, composition={1: 1, 8: 1})
    Hp = Solute('H+', charge=1, composition={1: 1})
    Kw = 1e-14/55.5
    w_auto_p = Equilibrium({H2O: 1}, {Hp: 1, OHm: 1}, Kw)
    return EqSys([w_auto_p], [H2O, OHm, Hp])


def test_EqSystemBase_1():
    es = _get_es1()
    assert es.stoichs.tolist() == [[-1, 1]]
    assert es.eq_constants() == [3]


@pytest.mark.xfail  # look into this later
def test_EqSystemLin_2():
    wes = _get_es_water(EqSystemLin)
    C = [55.5, 1e-7, 1e-7]
    fval, elim, elim_cbs = wes.f(C, C, rref=False)
    assert np.all(np.abs(fval) < 1e-16)
    fval, elim, elim_cbs = wes.f(1000*np.array(C), C,
                                 scaling=1000, rref=False)
    assert np.all(np.abs(fval) < 1e-16)


@pytest.mark.xfail
def test__solve_equilibrium_coord():
    c = np.array([-1.2882e-14, 3.1156e-10, 3.2099e-10, 9.679e-09, 5.5469e-07])
    stoich = np.array([1, 1, 0, 0, -1])
    K = 1e-22
    # K = (c[0] + r)*(c[1] + r)/(c[4] - r)
    # Kc[4] - Kr = c[0]c[1] + r(c[0] + c[1]) + r**2
    # {p = (c[0] + c[1] + K)/2}
    # {q = c[0]c[1] - Kc[4]}
    # r = p +/- sqrt(p*p/4 - q)
    # ... scrap that, there's a neg. conc
    _solve_equilibrium_coord(c, stoich, K)


def test_Equilibria_arithmetics():
    es1 = _get_es1()
    e, = es1.rxns
    e2 = 2*e
    sum2 = e + e
    assert sum2 == e2


def test_Equilibria_root():
    eqsys, c0 = get_ammonical_cupric_eqsys()
    x, sol = eqsys.root(c0)
    assert sol.success
