# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

try:
    import numpy as np
except ImportError:
    np = None

from chempy.arrhenius import ArrheniusParam
from chempy.chemistry import Substance, Reaction, ReactionSystem
from chempy.units import (
    SI_base_registry, get_derived_unit, allclose, units_library,
    to_unitless, default_units as u
)
from chempy.util._expr import Expr
from chempy.util.testing import requires
from .test_rates import _get_SpecialFraction_rsys
from ..rates import ArrheniusMassAction, Radiolytic
from .._rates import TPolyMassAction
from ..ode import get_odesys
from ..integrated import dimerization_irrev


@requires('numpy', 'pyodesys')
def test_get_odesys_1():
    k = .2
    a = Substance('A')
    b = Substance('B')
    r = Reaction({'A': 1}, {'B': 1}, param=k)
    rsys = ReactionSystem([r], [a, b])
    odesys = get_odesys(rsys, include_params=True)[0]
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


@requires('numpy', 'pyodesys')
def test_get_odesys_2():
    g = Radiolytic([3.14])
    a = Substance('A')
    b = Substance('B')
    r = Reaction({'A': 1}, {'B': 1}, param=g)
    rsys = ReactionSystem([r], [a, b])
    odesys = get_odesys(rsys, include_params=True)[0]
    c0 = {
        'A': 1.0,
        'B': 3.0,
    }
    t = np.linspace(0.0, .1)
    xout, yout, info = odesys.integrate(t, rsys.as_per_substance_array(c0), {'doserate': 2.72, 'density': .998})
    yref = np.zeros((t.size, 2))
    k = 3.14*2.72*.998
    yref[:, 0] = 1 - k*t
    yref[:, 1] = 3 + k*t
    assert np.allclose(yout, yref)


@requires(units_library, 'pyodesys')
def test_get_odesys_3():
    M = u.molar
    s = u.second
    mol = u.mol
    m = u.metre
    substances = list(map(Substance, 'H2O H+ OH-'.split()))
    dissociation = Reaction({'H2O': 1}, {'H+': 1, 'OH-': 1}, 2.47e-5/s)
    recombination = Reaction({'H+': 1, 'OH-': 1}, {'H2O': 1}, 1.37e11/M/s)
    rsys = ReactionSystem([dissociation, recombination], substances)
    odesys = get_odesys(
        rsys, include_params=True, unit_registry=SI_base_registry,
        output_conc_unit=M)[0]
    c0 = {'H2O': 55.4*M, 'H+': 1e-7*M, 'OH-': 1e-4*mol/m**3}
    x, y, p = odesys.pre_process(-42*u.second,
                                 rsys.as_per_substance_array(c0, unit=M))
    fout = odesys.f_cb(x, y, p)

    time_unit = get_derived_unit(SI_base_registry, 'time')
    conc_unit = get_derived_unit(SI_base_registry, 'concentration')

    r1 = to_unitless(55.4*2.47e-5*M/s, conc_unit/time_unit)
    r2 = to_unitless(1e-14*1.37e11*M/s, conc_unit/time_unit)
    assert abs(fout[0] - r2 + r1) < 1e-10
    assert abs(fout[1] - r1 + r2) < 1e-10
    assert abs(fout[2] - r1 + r2) < 1e-10


@requires(units_library, 'pyodesys')
def test_get_odesys__with_units():
    a = Substance('A')
    b = Substance('B')
    molar = u.molar
    second = u.second
    r = Reaction({'A': 2}, {'B': 1}, param=1e-3/molar/second)
    rsys = ReactionSystem([r], [a, b])
    odesys = get_odesys(rsys, include_params=True,
                        unit_registry=SI_base_registry)[0]
    c0 = {
        'A': 13*u.mol / u.metre**3,
        'B': .2 * u.molar
    }
    conc_unit = get_derived_unit(SI_base_registry, 'concentration')
    t = np.linspace(0, 10)*u.hour
    xout, yout, info = odesys.integrate(
        t, rsys.as_per_substance_array(c0, unit=conc_unit),
        atol=1e-10, rtol=1e-12)

    t_unitless = to_unitless(xout, u.second)
    Aref = dimerization_irrev(t_unitless, 1e-6, 13.0)
    # Aref = 1/(1/13 + 2*1e-6*t_unitless)
    yref = np.zeros((xout.size, 2))
    yref[:, 0] = Aref
    yref[:, 1] = 200 + (13-Aref)/2
    print((yout - yref*conc_unit)/yout)
    assert allclose(yout, yref*conc_unit)


@requires('numpy', 'pyodesys')
def test_SpecialFraction():
    k, kprime = 3.142, 2.718
    rsys = _get_SpecialFraction_rsys(k, kprime)

    odesys = get_odesys(rsys, include_params=True)[0]
    c0 = {'H2': 13, 'Br2': 17, 'HBr': 19}
    r = k*c0['H2']*c0['Br2']**(3/2)/(c0['Br2'] + kprime*c0['HBr'])
    ref = rsys.as_per_substance_array({'H2': -r, 'Br2': -r, 'HBr': 2*r})
    res = odesys.f_cb(0, rsys.as_per_substance_array(c0))
    assert np.allclose(res, ref)


@requires(units_library, 'pyodesys')
def test_SpecialFraction_with_units():
    k, kprime = 3.142 * u.s**-1 * u.molar**-0.5, 2.718
    rsys = _get_SpecialFraction_rsys(k, kprime)
    odesys = get_odesys(rsys, include_params=True,
                        unit_registry=SI_base_registry)[0]
    c0 = {'H2': 13*u.molar, 'Br2': 16*u.molar, 'HBr': 19*u.molar}
    r = k*c0['H2']*c0['Br2']**(3/2)/(c0['Br2'] + kprime*c0['HBr'])
    conc_unit = u.mol/u.metre**3
    rate_unit = conc_unit/u.second
    ref = rsys.as_per_substance_array({'H2': -r, 'Br2': -r, 'HBr': 2*r}, unit=rate_unit)
    res = odesys.f_cb(0, rsys.as_per_substance_array(c0, unit=conc_unit))
    assert allclose(to_unitless(ref, rate_unit), res)


@requires('pyodesys')
def test_ode_with_global_parameters():
    ratex = ArrheniusMassAction([1e10, 40e3/8.3145])
    rxn = Reaction({'A': 1}, {'B': 1}, ratex)
    rsys = ReactionSystem([rxn], 'A B')
    odesys = get_odesys(rsys, include_params=True)[0]
    conc = {'A': 3, 'B': 5}
    x, y, p = odesys.pre_process(-37, conc, {'temperature': 298.15})
    fout = odesys.f_cb(x, y, p)
    ref = 3*1e10*np.exp(-40e3/8.3145/298.15)
    assert abs((fout[0] + ref)/ref) < 1e-14
    assert abs((fout[1] - ref)/ref) < 1e-14


@requires('pyodesys')
def test_get_ode__ArrheniusParam():
    rxn = Reaction({'A': 1}, {'B': 1}, None)
    rxn.param = ArrheniusParam(1e10, 40e3)
    rsys = ReactionSystem([rxn], 'A B')
    odesys = get_odesys(rsys, include_params=True)[0]
    conc = {'A': 3, 'B': 5}
    x, y, p = odesys.pre_process(-37, conc, {'temperature': 200})
    fout = odesys.f_cb(x, y, p)
    ref = 3*1e10*np.exp(-40e3/8.314472/200)
    assert abs((fout[0] + ref)/ref) < 1e-14
    assert abs((fout[1] - ref)/ref) < 1e-14


@requires('pyodesys')
def test_get_ode__Radiolytic():
    rad = Radiolytic([2.4e-7])
    rxn = Reaction({'A': 4, 'B': 1}, {'C': 3, 'D': 2}, rad)
    rsys = ReactionSystem([rxn], 'A B C D')
    odesys = get_odesys(rsys, include_params=True)[0]
    c = {'A': 3, 'B': 5, 'C': 11, 'D': 13}
    x, y, p = odesys.pre_process(-37, c, {'doserate': 0.4, 'density': 0.998})
    fout = odesys.f_cb(x, y, p)
    r = 2.4e-7*0.4*0.998
    ref = [-4*r, -r, 3*r, 2*r]
    assert np.all(abs((fout - ref)/ref) < 1e-14)


@requires('pyodesys', units_library)
def test_get_ode__Radiolytic__units():
    rad = Radiolytic([2.4e-7*u.mol/u.joule])
    rxn = Reaction({'A': 4, 'B': 1}, {'C': 3, 'D': 2}, rad)
    rsys = ReactionSystem([rxn], 'A B C D')
    odesys = get_odesys(rsys, include_params=True,
                        unit_registry=SI_base_registry)[0]
    conc = {'A': 3*u.molar, 'B': 5*u.molar, 'C': 11*u.molar, 'D': 13*u.molar}
    x, y, p = odesys.pre_process(-37*u.second, conc, {
        'doserate': 0.4*u.gray/u.second,
        'density': 0.998*u.kg/u.decimetre**3
    })
    fout = odesys.f_cb(x, y, p)  # f_cb does not carry any units
    r = 2.4e-7*0.4*0.998*1e3  # mol/m3
    ref = [-4*r, -r, 3*r, 2*r]
    assert np.all(abs((fout - ref)/ref) < 1e-14)


class Density(Expr):
    """ Arguments: rho0 drhodT T0 """
    parameter_keys = ('temperature',)
    kw = {}

    def __call__(self, variables, backend=None):
        rho0, drhodT, T0 = self.all_args(variables)
        return rho0 + drhodT*(variables['temperature'] - T0)


@requires('pyodesys')
def test_get_ode__Radiolytic__substitutions():
    rad = Radiolytic([2.4e-7])
    rxn = Reaction({'A': 4, 'B': 1}, {'C': 3, 'D': 2}, rad)
    rsys = ReactionSystem([rxn], 'A B C D')
    substance_rho = Density([1, -1e-3, 273.15])
    odesys = get_odesys(rsys, include_params=True,
                        substitutions={'density': substance_rho})[0]
    conc = {'A': 3, 'B': 5, 'C': 11, 'D': 13}
    state = {'doserate': 0.4, 'temperature': 298.15}
    x, y, p = odesys.pre_process(-37, conc, state)
    fout = odesys.f_cb(x, y, p)
    r = 2.4e-7*0.4*substance_rho({'temperature': 298.15})
    ref = [-4*r, -r, 3*r, 2*r]
    assert np.all(abs((fout - ref)/ref) < 1e-14)


@requires('pyodesys', units_library)
def test_get_ode__Radiolytic__substitutions__units():
    rad = Radiolytic([2.4e-7*u.mol/u.joule])
    rxn = Reaction({'A': 4, 'B': 1}, {'C': 3, 'D': 2}, rad)
    rsys = ReactionSystem([rxn], 'A B C D')
    g_dm3 = u.gram / u.decimetre**3
    kg_dm3 = u.kg / u.decimetre**3
    substance_rho = Density([1*kg_dm3, -1*g_dm3/u.kelvin, 273.15*u.kelvin])
    odesys = get_odesys(rsys, include_params=True, unit_registry=SI_base_registry,
                        substitutions={'density': substance_rho})[0]
    conc = {'A': 3*u.molar, 'B': 5*u.molar, 'C': 11*u.molar, 'D': 13*u.molar}
    x, y, p = odesys.pre_process(-37*u.second, conc, {'doserate': 0.4*u.gray/u.second, 'temperature': 298.15*u.kelvin})
    fout = odesys.f_cb(x, y, p)
    r = 2.4e-7*0.4*0.975 * 1e3  # mol/m3/s
    ref = [-4*r, -r, 3*r, 2*r]
    assert np.all(abs((fout - ref)/ref) < 1e-14)


@requires('pyodesys', units_library)
def test_get_ode__TPoly():
    rate = TPolyMassAction([273.15*u.K, 10/u.molar/u.s, 2/u.molar/u.s/u.K])
    rxn = Reaction({'A': 1, 'B': 1}, {'C': 3, 'D': 2}, rate, {'A': 3})
    rsys = ReactionSystem([rxn], 'A B C D')
    odesys = get_odesys(rsys, unit_registry=SI_base_registry)[0]
    conc = {'A': 3*u.molar, 'B': 5*u.molar, 'C': 11*u.molar, 'D': 13*u.molar}
    x, y, p = odesys.pre_process(-37*u.second, conc, {'temperature': 298.15*u.kelvin})
    fout = odesys.f_cb(x, y, p)
    r = 3*5*(10+2*25)*1000  # mol/m3/s
    ref = [-4*r, -r, 3*r, 2*r]
    assert np.all(abs((fout - ref)/ref) < 1e-14)
