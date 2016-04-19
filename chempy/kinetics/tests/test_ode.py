# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

try:
    import numpy as np
except ImportError:
    np = None

from chempy.chemistry import Substance, Reaction, ReactionSystem
from chempy.units import (
    default_units, SI_base_registry, get_derived_unit, allclose, units_library,
    to_unitless
)
from chempy.util.testing import requires
from .test_rates import _get_SpecialFraction_rsys
from ..rates import ArrheniusMassAction
from ..ode import get_odesys
from ..integrated import dimerization_irrev


@requires('numpy', 'pyodesys')
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


@requires(units_library, 'pyodesys')
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


@requires(units_library, 'pyodesys')
def test_get_odesys_2():
    M = default_units.molar
    s = default_units.second
    mol = default_units.mol
    m = default_units.metre
    substances = list(map(Substance, 'H2O H+ OH-'.split()))
    dissociation = Reaction({'H2O': 1}, {'H+': 1, 'OH-': 1}, 2.47e-5/s)
    recombination = Reaction({'H+': 1, 'OH-': 1}, {'H2O': 1}, 1.37e11/M/s)
    rsys = ReactionSystem([dissociation, recombination], substances)
    odesys = get_odesys(
        rsys, include_params=True, unit_registry=SI_base_registry,
        output_conc_unit=M)
    c0 = {'H2O': 55.4*M, 'H+': 1e-7*M, 'OH-': 1e-4*mol/m**3}
    x, y, p = odesys.pre_process(-42*default_units.second,
                                 rsys.as_per_substance_array(c0, unit=M))
    fout = odesys.f_cb(x, y, p)

    time_unit = get_derived_unit(SI_base_registry, 'time')
    conc_unit = get_derived_unit(SI_base_registry, 'concentration')

    r1 = to_unitless(55.4*2.47e-5*M/s, conc_unit/time_unit)
    r2 = to_unitless(1e-14*1.37e11*M/s, conc_unit/time_unit)
    assert abs(fout[0] - r2 + r1) < 1e-10
    assert abs(fout[1] - r1 + r2) < 1e-10
    assert abs(fout[2] - r1 + r2) < 1e-10


@requires('numpy', 'pyodesys')
def test_SpecialFraction():
    k, kprime = 3.142, 2.718
    rsys = _get_SpecialFraction_rsys(k, kprime)

    odesys = get_odesys(rsys, include_params=True)
    c0 = {'H2': 13, 'Br2': 17, 'HBr': 19}
    r = k*c0['H2']*c0['Br2']**(3/2)/(c0['Br2'] + kprime*c0['HBr'])
    ref = rsys.as_per_substance_array({'H2': -r, 'Br2': -r, 'HBr': 2*r})
    res = odesys.f_cb(0, rsys.as_per_substance_array(c0))
    assert np.allclose(res, ref)


@requires(units_library, 'pyodesys')
def test_SpecialFraction_with_units():
    u = default_units
    k, kprime = 3.142 * u.s**-1 * u.molar**-0.5, 2.718
    rsys = _get_SpecialFraction_rsys(k, kprime)
    odesys = get_odesys(rsys, include_params=True,
                        unit_registry=SI_base_registry)
    c0 = {'H2': 13*u.molar, 'Br2': 16*u.molar, 'HBr': 19*u.molar}
    r = k*c0['H2']*c0['Br2']**(3/2)/(c0['Br2'] + kprime*c0['HBr'])
    ref = rsys.as_per_substance_array({'H2': -r, 'Br2': -r, 'HBr': 2*r})
    res = odesys.f_cb(0, rsys.as_per_substance_array(c0, unit=u.molar))
    assert allclose(ref, res)


@requires('pyodesys')
def test_ode_with_global_parameters():
    ratex = ArrheniusMassAction([1e10, 40e3/8.3145])
    rxn = Reaction({'A': 1}, {'B': 1}, ratex)
    rsys = ReactionSystem([rxn], 'A B')
    odesys = get_odesys(rsys, include_params=True)
    conc = {'A': 3, 'B': 5}
    x, y, p = odesys.pre_process(-37, conc, {'temperature': 298.15})
    fout = odesys.f_cb(x, y, p)
    ref = 3*1e10*np.exp(-40e3/8.3145/298.15)
    assert abs((fout[0] + ref)/ref) < 1e-14
    assert abs((fout[1] - ref)/ref) < 1e-14
