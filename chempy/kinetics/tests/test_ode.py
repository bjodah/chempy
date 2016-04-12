# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import numpy as np

from chempy.chemistry import Substance, Reaction, ReactionSystem
from chempy.units import (
    default_units, SI_base_registry, get_derived_unit, allclose
)
from ..ode import get_odesys
from ..integrated import dimerization_irrev

from .test_rates import _get_h2_br2_rsys


def test_get_odesys_1():
    k = .2
    a = Substance('A')
    b = Substance('B')
    r = Reaction({'A': 1}, {'B': 1}, param=k)
    rsys = ReactionSystem([r], [a, b])
    odesys = get_odesys(rsys, include_params=True)
    c0 = {
        'A': 1.0,
        'B': 3.0,
    }
    t = np.linspace(0.0, 10.0)
    xout, yout, info = odesys.integrate(t, rsys.as_per_substance_array(c0))
    yref = np.zeros((t.size, 2))
    yref[:, 0] = np.exp(-k*t)
    yref[:, 1] = 4 - np.exp(-k*t)
    assert np.allclose(yout, yref)


def test_get_odesys__with_units():
    a = Substance('A')
    b = Substance('B')
    molar = default_units.molar
    second = default_units.second
    r = Reaction({'A': 2}, {'B': 1}, param=1e-3/molar/second)
    rsys = ReactionSystem([r], [a, b])
    odesys = get_odesys(rsys, include_params=True,
                        unit_registry=SI_base_registry)
    c0 = {
        'A': 13*default_units.mol / default_units.metre**3,
        'B': .2 * default_units.molar
    }
    conc_unit = get_derived_unit(SI_base_registry, 'concentration')
    t = np.linspace(0, 10)*default_units.hour
    xout, yout, info = odesys.integrate(
        t, rsys.as_per_substance_array(c0, unit=conc_unit))
    Aref = dimerization_irrev(3600, 1e-6, 13.0)
    yref = np.zeros((xout.size, 2))
    yref[:, 0] = Aref
    yref[:, 1] = .2e-3 + 2*Aref
    assert allclose(yout, yref*conc_unit)


def test_get_odesys_2():
    M = default_units.molar
    s = default_units.second
    substances = list(map(Substance, 'H2O H+ OH-'.split()))
    dissociation = Reaction({'H2O': 1}, {'H+': 1, 'OH-': 1}, 2.47e-5/s)
    recombination = Reaction({'H+': 1, 'OH-': 1}, {'H2O': 1}, 1.37e11/M/s)
    rsys = ReactionSystem([dissociation, recombination], substances)
    get_odesys(rsys, include_params=True,
               unit_registry=SI_base_registry, output_conc_unit=M)


def test_h2_br2():
    k, kprime = 3.142, 2.718
    rsys = _get_h2_br2_rsys(k, kprime)

    odesys = get_odesys(rsys, include_params=True)
    c0 = {'H2': 13, 'Br2': 17, 'HBr': 19}
    r = k*c0['H2']*c0['Br2']**(3/2)/(c0['Br2'] + kprime*c0['HBr'])
    ref = rsys.as_per_substance_array({'H2': -r, 'Br2': -r, 'HBr': 2*r})
    res = odesys.f_cb(0, rsys.as_per_substance_array(c0))
    assert np.allclose(res, ref)


def test_h2_br2_with_units():
    u = default_units
    k, kprime = 3.142 * u.s**-1 * u.molar**-0.5, 2.718
    rsys = _get_h2_br2_rsys(k, kprime)
    odesys = get_odesys(rsys, include_params=True,
                        unit_registry=SI_base_registry)
    c0 = {'H2': 13*u.molar, 'Br2': 16*u.molar, 'HBr': 19*u.molar}
    r = k*c0['H2']*c0['Br2']**(3/2)/(c0['Br2'] + kprime*c0['HBr'])
    ref = rsys.as_per_substance_array({'H2': -r, 'Br2': -r, 'HBr': 2*r})
    res = odesys.f_cb(0, rsys.as_per_substance_array(c0, unit=u.molar))
    assert allclose(ref, res)
