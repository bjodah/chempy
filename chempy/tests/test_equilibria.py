# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import collections

import pytest
import numpy as np
from scipy.optimize import fsolve

from ..chemistry import (
    Substance, Reaction, Equilibrium, Species, Solute,
)

from ..equilibria import (
    equilibrium_quotient, equilibrium_residual, get_rc_interval,
    solve_equilibrium, prodpow, EqSystem, NumSysLin, NumSysLog
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


def test_EqSystem():
    a, b = sbstncs = Substance('a'), Substance('b')
    rxns = [Reaction({'a': 1}, {'b': 1})]
    es = EqSystem(rxns, [(s.name, s) for s in sbstncs])
    assert es.net_stoichs().tolist() == [[-1, 1]]


def _get_es1():
    a, b = sbstncs = Solute('a'), Solute('b')
    rxns = [Equilibrium({'a': 1}, {'b': 1}, 3)]
    return EqSystem(rxns, sbstncs)


def _get_es_water(EqSys=EqSystem):
    H2O = Substance('H2O', charge=0, composition={1: 2, 8: 1})
    OHm = Substance('OH-', charge=-1, composition={1: 1, 8: 1})
    Hp = Substance('H+', charge=1, composition={1: 1})
    Kw = 1e-14/55.5
    w_auto_p = Equilibrium({'H2O': 1}, {'Hp': 1, 'OHm': 1}, Kw)
    return EqSys([w_auto_p], [H2O, OHm, Hp])


def test_EqSystem_1():
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
    x, sol, sane = eqsys.root(c0, NumSys=(NumSysLog,))
    assert sol['success'] and sane


def _solutes(Cls):
    return (
        Cls('H2O', 0, composition={1: 2, 8: 1}),
        Cls('H+', 1, composition={1: 1}),
        Cls('OH-', -1, composition={1: 1, 8: 1}),
        Cls('NH4+', 1, composition={1: 4, 7: 1}),
        Cls('NH3', 0, composition={1: 3, 7: 1})
    )


@pytest.mark.parametrize('Cls', [Solute, Species])
def test_Equilibria_root_simple(Cls):
    solutes = water, hydronium, hydroxide, ammonium, ammonia = _solutes(Cls)

    water_auto_protolysis = Equilibrium(
        {water.name: 1}, {hydronium.name: 1, hydroxide.name: 1}, 1e-14/55.5)
    ammonia_protolysis = Equilibrium(
        {ammonium.name: 1}, {hydronium.name: 1, ammonia.name: 1},
        10**-9.26/55.5
    )
    eqsys = EqSystem([water_auto_protolysis, ammonia_protolysis], solutes)

    init_concs = collections.defaultdict(float, {
        water.name: 55.5, ammonia.name: 1e-3})
    x, sol1, sane1 = eqsys.root(init_concs)
    x, sol2, sane2 = eqsys.root(init_concs, x, NumSys=(NumSysLog,))
    assert sane2
    ref = eqsys.as_per_substance_array({
        water.name: 55.5,
        ammonia.name: 1e-3 - 6.2e-4,
        ammonium.name: 6.2e-4,
        hydronium.name: 1.6e-11,
        hydroxide.name: 6.2e-4
    })
    assert np.allclose(x, ref, rtol=0.02, atol=1e-16)


def _get_NaCl(Cls, **precipitate_kwargs):
    Na_p, Cl_m, NaCl = sbstncs = (
        Cls('Na+', 1, composition={11: 1}),
        Cls('Cl-', -1, composition={17: 1}),
        Cls('NaCl', composition={11: 1, 17: 1}, **precipitate_kwargs)
    )
    sp = Equilibrium({'NaCl': 1}, {'Na+': 1, 'Cl-': 1}, 4.0)
    eqsys = EqSystem([sp], sbstncs)
    cases = [
        [(0, 0, .1), (.1, .1, 0)],
        [(.5, .5, .4), (.9, .9, 0)],
        [(1, 1, 1), (2, 2, 0)],
        [(0, 0, 2), (2, 2, 0)],
        [(3, 3, 3), (2, 2, 4)],
        [(3, 3, 0), (2, 2, 1)],
        [(0, 0, 3), (2, 2, 1)],
    ]
    return eqsys, [s.name for s in sbstncs], cases


@pytest.mark.parametrize('kwargs', [
    dict(Cls=Solute, precipitate=True),
    dict(Cls=Species, phase_idx=1)
])
def test_EqSystem_dissolved(kwargs):
    eqsys, names, _ = _get_NaCl(**kwargs)
    inp = eqsys.as_per_substance_array({'Na+': 1, 'Cl-': 2, 'NaCl': 4})
    result = eqsys.dissolved(inp)
    ref = eqsys.as_per_substance_array({'Na+': 5, 'Cl-': 6, 'NaCl': 0})
    assert np.allclose(result, ref)


@pytest.mark.parametrize('NumSys', [(NumSysLin,), (NumSysLog,),
                                    (NumSysLog, NumSysLin)])
def test_precipitate(NumSys):
    eqsys, species, cases = _get_NaCl(Solute, precipitate=True)

    for init, final in cases:
        x, sol, sane = eqsys.root(dict(zip(species, init)),
                                  NumSys=NumSys, rref_preserv=True)
        assert sol['success'] and sane
        assert x is not None
        assert np.allclose(x, np.asarray(final))
