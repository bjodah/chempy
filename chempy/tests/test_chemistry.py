# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math
from operator import attrgetter

import pytest

from ..units import default_units
from ..chemistry import (
    Substance, Species, Reaction, ReactionSystem, ArrheniusRate,
    ArrheniusRateWithUnits, Equilibrium
)


def test_Substance():
    s = Substance.from_formula('H+')
    assert s.composition == {0: 1, 1: 1}
    assert s.charge == 1
    assert abs(s.mass - 1.008) < 1e-3


def test_Substance__2():
    H2O = Substance(name='H2O',  charge=0, latex_name=r'$\mathrm{H_{2}O}$',
                    other_properties={'pKa': 14})
    OH_m = Substance(name='OH-',  charge=-1, latex_name=r'$\mathrm{OH^{-}}$')
    assert sorted([OH_m, H2O], key=attrgetter('name')) == [H2O, OH_m]


def test_Species():
    s = Species.from_formula('H2O')
    assert s.phase_idx == 0
    mapping = {'(aq)': 0, '(s)': 1, '(g)': 2}
    assert Species.from_formula('CO2(g)').phase_idx == 3
    assert Species.from_formula('CO2(g)', mapping).phase_idx == 2
    assert Species.from_formula('CO2(aq)', mapping).phase_idx == 0
    assert Species.from_formula('NaCl(s)').phase_idx == 1
    assert Species.from_formula('NaCl(s)', phase_idx=7).phase_idx == 7
    assert Species.from_formula('CO2(aq)', mapping, phase_idx=7).phase_idx == 7


def test_Solute():
    from ..chemistry import Solute
    from ..util.pyutil import ChemPyDeprecationWarning
    with pytest.warns(ChemPyDeprecationWarning):
        w = Solute('H2O')
    assert w.name == 'H2O'


def test_Reaction():
    substances = s_Hp, s_OHm, s_H2O = (
        Substance('H+', composition={0: 1, 1: 1}),
        Substance('OH-', composition={0: -1, 1: 1, 8: 1}),
        Substance('H2O', composition={0: 0, 1: 2, 8: 1}),
    )
    substance_names = Hp, OHm, H2O = [s.name for s in substances]
    substance_dict = {n: s for n, s in zip(substance_names, substances)}
    r1 = Reaction({Hp: 1, OHm: 1}, {H2O: 1})
    assert sum(r1.composition_violation(substance_dict)) == 0
    assert r1.charge_neutrality_violation(substance_dict) == 0

    r2 = Reaction({Hp: 1, OHm: 1}, {H2O: 2})
    assert sum(r2.composition_violation(substance_dict)) != 0
    assert r2.charge_neutrality_violation(substance_dict) == 0

    r3 = Reaction({Hp: 2, OHm: 1}, {H2O: 2})
    assert sum(r3.composition_violation(substance_dict)) != 0
    assert r3.charge_neutrality_violation(substance_dict) != 0


def test_ReactionSystem__as_per_substance_array():
    mol = default_units.mol
    m = default_units.metre
    M = default_units.molar
    rs = ReactionSystem([], [Substance('H2O')])
    c = rs.as_per_substance_array({'H2O': 1*M},
                                  unit=M)
    assert c.dimensionality == M.dimensionality
    assert abs(c[0]/(1000*mol/m**3) - 1) < 1e-16

    c = rs.as_per_substance_array({'H2O': 1})
    with pytest.raises(KeyError):
        c = rs.as_per_substance_array({'H': 1})


def test_ArrheniusRate():
    k = ArrheniusRate(1e10, 42e3)(273.15)
    ref = 1e10 * math.exp(-42e3/(8.3145*273.15))
    assert abs((k - ref)/ref) < 1e-4


def test_ArrheniusRateWithUnits():
    s = default_units.second
    mol = default_units.mol
    J = default_units.joule
    K = default_units.kelvin
    k = ArrheniusRateWithUnits(1e10/s, 42e3 * J/mol)(273.15*K)
    ref = 1e10/s * math.exp(-42e3/(8.3145*273.15))
    assert abs((k - ref)/ref) < 1e-4


def test_Equilibrium__as_reactions():
    s = default_units.second
    M = default_units.molar
    H2O, Hp, OHm = map(Substance, 'H2O H+ OH-'.split())
    eq = Equilibrium({'H2O': 1}, {'H+': 1, 'OH-': 1}, 1e-14)
    rate = 1.31e11/M/s
    fw, bw = eq.as_reactions(kb=rate, units=default_units)
    assert abs((bw.param - rate)/rate) < 1e-15
    assert abs((fw.param / M)/bw.param - 1e-14)/1e-14 < 1e-15


def test_Reaction__from_string():
    r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", 'H2O H+ OH-'.split())
    assert r.reac == {'H2O': 1} and r.prod == {'H+': 1, 'OH-': 1}

    with pytest.raises(ValueError):
        Reaction.from_string("H2O -> H+ + OH-; 1e-4", 'H2O H OH-'.split())
