import pytest
import numpy as np
from scipy.optimize import fsolve

from ..chemistry import (
    Solute, Substance, Reaction
)

from ..equilibria import (
    equilibrium_quotient, equilibrium_residual, get_rc_interval,
    solve_equilibrium, EqSystemBase, prodpow, _solve_equilibrium_coord,
    Equilibrium, EqSystemLin, EqSystemLog
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
    assert es.net_stoichs().tolist() == [[-1, 1]]


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
    assert es.stoichs().tolist() == [[-1, 1]]
    assert es.eq_constants() == [3]


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


@pytest.mark.parametrize('EqSys', [EqSystemLin])
def test_solid(EqSys):
    Na_p, Cl_m, NaCl = sbstncs = (
        Solute('Na+', 1, composition={11: 1}),
        Solute('Cl-', -1, composition={17: 1}),
        Solute('NaCl', composition={11: 1, 17: 1}, solid=True)
    )
    sp = Equilibrium({NaCl: 1}, {Na_p: 1, Cl_m: 1}, 4.0)
    eqsys = EqSys([sp], sbstncs)

    x, sol = eqsys.root({Na_p: 1.0, Cl_m: 1.0, NaCl: 0.0})
    assert np.allclose(x, [1, 1, 0])

    x, sol = eqsys.root({Na_p: 0, Cl_m: 0, NaCl: 7.0})
    assert np.allclose(x, [2, 2, 5])

    x, sol = eqsys.root({Na_p: 0.5, Cl_m: 0.5, NaCl: 0.5})
    assert np.allclose(x, [1, 1, 0])

    x, sol = eqsys.root({Na_p: 0, Cl_m: 0, NaCl: 1.0})
    assert np.allclose(x, [1, 1, 0])

    # x, sol = eqsys.root({Na_p: 0, Cl_m: 0, NaCl: 3.5})
    # assert np.allclose(x, [2.5, 2.5, 1])
