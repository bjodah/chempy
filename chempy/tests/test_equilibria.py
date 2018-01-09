# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import collections

import pytest
try:
    import numpy as np
except ImportError:
    np, NumSysLin, NumSysLog = [None]*3
else:
    from ..equilibria import EqSystem, NumSysLin, NumSysLog
    from .ammonical_cupric_solution import get_ammonical_cupric_eqsys

from ..util.testing import requires
from ..chemistry import (
    Substance, Reaction, Equilibrium, Species
)


@requires('numpy')
def test_EqSystem():
    a, b = sbstncs = Substance('a'), Substance('b')
    rxns = [Reaction({'a': 1}, {'b': 1})]
    es = EqSystem(rxns, collections.OrderedDict(
        [(s.name, s) for s in sbstncs]))
    assert es.net_stoichs().tolist() == [[-1, 1]]


def _get_es1():
    a, b = sbstncs = Species('a'), Species('b')
    rxns = [Equilibrium({'a': 1}, {'b': 1}, 3)]
    return EqSystem(rxns, sbstncs)


def _get_es_water(EqSys=None):
    if EqSys is None:
        EqSys = EqSystem
    H2O = Substance('H2O', charge=0, composition={1: 2, 8: 1})
    OHm = Substance('OH-', charge=-1, composition={1: 1, 8: 1})
    Hp = Substance('H+', charge=1, composition={1: 1})
    Kw = 1e-14/55.5
    w_auto_p = Equilibrium({'H2O': 1}, {'Hp': 1, 'OHm': 1}, Kw)
    return EqSys([w_auto_p], [H2O, OHm, Hp])


@requires('numpy')
def test_EqSystem_1():
    es = _get_es1()
    assert es.stoichs().tolist() == [[-1, 1]]
    assert es.eq_constants() == [3]


@requires('numpy')
def test_Equilibria_arithmetics():
    es1 = _get_es1()
    e, = es1.rxns
    e2 = 2*e
    sum2 = e + e
    assert sum2 == e2


@pytest.mark.slow
@requires('numpy')
def test_Equilibria_root():
    eqsys, c0 = get_ammonical_cupric_eqsys()
    x, sol, sane = eqsys.root(c0, NumSys=(NumSysLog,))
    assert sol['success'] and sane


def _species(Cls):
    return (
        Cls('H2O', 0, composition={1: 2, 8: 1}),
        Cls('H+', 1, composition={1: 1}),
        Cls('OH-', -1, composition={1: 1, 8: 1}),
        Cls('NH4+', 1, composition={1: 4, 7: 1}),
        Cls('NH3', 0, composition={1: 3, 7: 1})
    )


@requires('numpy')
def test_Equilibria_root_simple():
    species = water, hydronium, hydroxide, ammonium, ammonia = _species(Species)

    water_auto_protolysis = Equilibrium(
        {water.name: 1}, {hydronium.name: 1, hydroxide.name: 1}, 1e-14/55.5)
    ammonia_protolysis = Equilibrium(
        {ammonium.name: 1}, {hydronium.name: 1, ammonia.name: 1},
        10**-9.26/55.5
    )
    eqsys = EqSystem([water_auto_protolysis, ammonia_protolysis], species)

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


@requires('numpy')
def test_EqSystem_dissolved():
    eqsys, names, _ = _get_NaCl(Cls=Species, phase_idx=1)
    inp = eqsys.as_per_substance_array({'Na+': 1, 'Cl-': 2, 'NaCl': 4})
    result = eqsys.dissolved(inp)
    ref = eqsys.as_per_substance_array({'Na+': 5, 'Cl-': 6, 'NaCl': 0})
    assert np.allclose(result, ref)


@requires('numpy')
@pytest.mark.parametrize('NumSys', [(NumSysLin,), (NumSysLog,),
                                    (NumSysLog, NumSysLin)])
def test_precipitate(NumSys):
    eqsys, species, cases = _get_NaCl(Species, phase_idx=1)

    for init, final in cases:
        x, sol, sane = eqsys.root(dict(zip(species, init)),
                                  NumSys=NumSys, rref_preserv=True, tol=1e-12)
        assert sol['success'] and sane
        assert x is not None
        assert np.allclose(x, np.asarray(final))
