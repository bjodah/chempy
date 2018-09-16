# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from collections import defaultdict, OrderedDict
from itertools import permutations
import math

import pytest
try:
    import numpy as np
except ImportError:
    np = None


from chempy import Equilibrium, Reaction, ReactionSystem, Substance
from chempy.thermodynamics.expressions import MassActionEq
from chempy.units import (
    SI_base_registry, get_derived_unit, allclose, units_library,
    to_unitless, default_constants as const, default_units as u
)
from chempy.util._expr import Expr
from chempy.util.testing import requires
from .test_rates import _get_SpecialFraction_rsys
from ..arrhenius import ArrheniusParam
from ..rates import Arrhenius, MassAction, Radiolytic, RampedTemp
from .._rates import ShiftedTPoly
from ..ode import get_odesys, chained_parameter_variation, _create_odesys as create_odesys
from ..integrated import dimerization_irrev, binary_rev


@requires('numpy', 'pyodesys')
def test_get_odesys_1():
    k = .2
    a = Substance('A')
    b = Substance('B')
    r = Reaction({'A': 1}, {'B': 1}, param=k)
    rsys = ReactionSystem([r], [a, b])
    assert sorted(rsys.substances.keys()) == ['A', 'B']
    odesys = get_odesys(rsys, include_params=True)[0]
    c0 = {
        'A': 1.0,
        'B': 3.0,
    }
    t = np.linspace(0.0, 10.0)
    xout, yout, info = odesys.integrate(t, c0)
    yref = np.zeros((t.size, 2))
    yref[:, 0] = np.exp(-k*t)
    yref[:, 1] = 4 - np.exp(-k*t)
    assert np.allclose(yout, yref)


@requires('numpy', 'pyodesys', 'sympy')
def test_get_odesys__rate_exprs_cb():
    k = .2
    a = Substance('A')
    b = Substance('B')
    r = Reaction({'A': 1}, {'B': 1}, param=k)
    rsys = ReactionSystem([r], [a, b])
    assert sorted(rsys.substances.keys()) == ['A', 'B']
    odesys, extra = get_odesys(rsys)
    c0 = {'A': 1.0, 'B': 3.0}
    t = np.linspace(0.0, 10.0)
    res = odesys.integrate(t, c0)
    yref = np.zeros((t.size, 2))
    yref[:, 0] = np.exp(-k*t)
    yref[:, 1] = 4 - np.exp(-k*t)
    assert np.allclose(res.yout, yref)
    rate = extra['rate_exprs_cb'](res.xout, res.yout, res.params)
    assert np.allclose(rate[:, 0], k*yref[:, 0])


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
    xout, yout, info = odesys.integrate(t, rsys.as_per_substance_array(c0),
                                        {'doserate': 2.72, 'density': .998})
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
    x, y, p = odesys.to_arrays(-42*u.second,
                               rsys.as_per_substance_array(c0, unit=M), ())
    fout = odesys.f_cb(x, y, p)

    time_unit = get_derived_unit(SI_base_registry, 'time')
    conc_unit = get_derived_unit(SI_base_registry, 'concentration')

    r1 = to_unitless(55.4*2.47e-5*M/s, conc_unit/time_unit)
    r2 = to_unitless(1e-14*1.37e11*M/s, conc_unit/time_unit)
    assert np.all(abs(fout[:, 0] - r2 + r1)) < 1e-10
    assert np.all(abs(fout[:, 1] - r1 + r2)) < 1e-10
    assert np.all(abs(fout[:, 2] - r1 + r2)) < 1e-10


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
    ratex = MassAction(Arrhenius([1e10, 40e3/8.3145]))
    rxn = Reaction({'A': 1}, {'B': 1}, ratex)
    rsys = ReactionSystem([rxn], 'A B')
    odesys, extra = get_odesys(rsys, include_params=False)
    param_keys, unique_keys, p_units = map(extra.get, 'param_keys unique p_units'.split())
    conc = {'A': 3, 'B': 5}
    x, y, p = odesys.to_arrays(-37, conc, {'temperature': 298.15})
    fout = odesys.f_cb(x, y, p)
    ref = 3*1e10*np.exp(-40e3/8.3145/298.15)
    assert np.all(abs((fout[:, 0] + ref)/ref) < 1e-14)
    assert np.all(abs((fout[:, 1] - ref)/ref) < 1e-14)


@requires('pyodesys')
def test_get_ode__ArrheniusParam():
    rxn = Reaction({'A': 1}, {'B': 1}, None)
    rxn.param = ArrheniusParam(1e10, 40e3)
    rsys = ReactionSystem([rxn], 'A B')
    odesys = get_odesys(rsys, include_params=True)[0]
    conc = {'A': 3, 'B': 5}
    x, y, p = odesys.to_arrays(-37, conc, {'temperature': 200})
    fout = odesys.f_cb(x, y, p)
    ref = 3*1e10*np.exp(-40e3/8.314472/200)
    assert np.all(abs((fout[:, 0] + ref)/ref) < 1e-14)
    assert np.all(abs((fout[:, 1] - ref)/ref) < 1e-14)


@requires('pyodesys')
def test_get_ode__Radiolytic():
    rad = Radiolytic([2.4e-7])
    rxn = Reaction({'A': 4, 'B': 1}, {'C': 3, 'D': 2}, rad)
    rsys = ReactionSystem([rxn], 'A B C D')
    odesys = get_odesys(rsys, include_params=True)[0]
    c = {'A': 3, 'B': 5, 'C': 11, 'D': 13}
    x, y, p = odesys.to_arrays(-37, c, {'doserate': 0.4, 'density': 0.998})
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
    x, y, p = odesys.to_arrays(-37*u.second, conc, {
        'doserate': 0.4*u.gray/u.second,
        'density': 0.998*u.kg/u.decimetre**3
    })
    fout = odesys.f_cb(x, y, p)  # f_cb does not carry any units
    r = 2.4e-7*0.4*0.998*1e3  # mol/m3
    ref = [-4*r, -r, 3*r, 2*r]
    assert np.all(abs((fout - ref)/ref) < 1e-14)


@requires('pyodesys', units_library)
def test_get_ode__Radiolytic__units__multi():
    rad = Radiolytic([2.4e-7*u.mol/u.joule])
    rxn = Reaction({'A': 4, 'B': 1}, {'C': 3, 'D': 2}, rad)
    rsys = ReactionSystem([rxn], 'A B C D')
    odesys = get_odesys(rsys, include_params=True,
                        unit_registry=SI_base_registry)[0]
    conc = {'A': 3*u.molar, 'B': 5*u.molar, 'C': 11*u.molar, 'D': 13*u.molar}
    doserates = [dr*u.gray/u.second for dr in [.1, .2, .3, .4]]
    results = odesys.integrate(37*u.second, conc, {
        'doserate': doserates,
        'density': 0.998*u.kg/u.decimetre**3
    })
    assert len(results) == 4
    for i, r in enumerate(results):
        dr = r.params[odesys.param_names.index('doserate')]
        assert dr.ndim == 0 or len(dr) == 1
        assert dr == doserates[i]


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
    x, y, p = odesys.to_arrays(-37, conc, state)
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
    x, y, p = odesys.to_arrays(-37*u.second, conc, {'doserate': 0.4*u.gray/u.second, 'temperature': 298.15*u.kelvin})
    fout = odesys.f_cb(x, y, p)
    r = 2.4e-7*0.4*0.975 * 1e3  # mol/m3/s
    ref = [-4*r, -r, 3*r, 2*r]
    assert np.all(abs((fout - ref)/ref) < 1e-14)


@requires('pyodesys', units_library)
def test_get_ode__TPoly():
    rate = MassAction(ShiftedTPoly([273.15*u.K, 10/u.molar/u.s, 2/u.molar/u.s/u.K]))
    rxn = Reaction({'A': 1, 'B': 1}, {'C': 3, 'D': 2}, rate, {'A': 3})
    rsys = ReactionSystem([rxn], 'A B C D')
    odesys = get_odesys(rsys, unit_registry=SI_base_registry)[0]
    conc = {'A': 3*u.molar, 'B': 5*u.molar, 'C': 11*u.molar, 'D': 13*u.molar}
    x, y, p = odesys.to_arrays(-37*u.second, conc, {'temperature': 298.15*u.kelvin})
    fout = odesys.f_cb(x, y, p)
    r = 3*5*(10+2*25)*1000  # mol/m3/s
    ref = [-4*r, -r, 3*r, 2*r]
    assert np.all(abs((fout - ref)/ref) < 1e-14)


@requires('pyodesys', units_library)
def test_get_odesys__time_dep_rate():

    class RampedRate(Expr):
        argument_names = ('rate_constant', 'ramping_rate')

        def __call__(self, variables, reaction, backend=math):
            rate_constant, ramping_rate = self.all_args(variables, backend=backend)
            return rate_constant * ramping_rate * variables['time']

    rate = MassAction(RampedRate([7, 2]))
    rxn = Reaction({'A': 1}, {'B': 3}, rate)
    rsys = ReactionSystem([rxn], 'A B')
    odesys = get_odesys(rsys)[0]
    conc = {'A': 3, 'B': 11}
    x, y, p = odesys.to_arrays([5, 13, 17], conc, ())
    fout = odesys.f_cb(x, y, p)
    r = 2*7*3
    ref = np.array([
        [-r*5, -r*13, -r*17],
        [r*5*3, r*13*3, r*17*3]
    ]).T
    assert np.allclose(fout, ref)


@pytest.mark.slow
@requires('pyodesys', units_library)
def test_get_odesys__time_dep_temperature():
    import sympy as sp

    def refA(t, A0, A, Ea_over_R, T0, dTdt):
        T = (T0 + dTdt*t)
        d_Ei = sp.Ei(-Ea_over_R/T0).n(100).round(90) - sp.Ei(-Ea_over_R/T).n(100).round(90)
        d_Texp = T0*sp.exp(-Ea_over_R/T0) - T*sp.exp(-Ea_over_R/T)
        return A0*sp.exp(A/dTdt*(Ea_over_R*d_Ei + d_Texp)).n(30)

    params = A0, A, Ea_over_R, T0, dTdt = [13, 1e10, 56e3/8, 273, 2]
    B0 = 11
    rate = MassAction(Arrhenius([A, Ea_over_R]))
    rxn = Reaction({'A': 1}, {'B': 3}, rate)
    rsys = ReactionSystem([rxn], 'A B')
    rt = RampedTemp([T0, dTdt], ('init_temp', 'ramp_rate'))
    odesys, extra = get_odesys(rsys, False, substitutions={'temperature': rt})
    all_pk, unique, p_units = map(extra.get, 'param_keys unique p_units'.split())
    conc = {'A': A0, 'B': B0}
    tout = [2, 5, 10]

    for ramp_rate in [2, 3, 4]:
        unique['ramp_rate'] = ramp_rate
        xout, yout, info = odesys.integrate(10, conc, unique, atol=1e-10, rtol=1e-12)
        params[-1] = ramp_rate
        Aref = np.array([float(refA(t, *params)) for t in xout])
        # Aref = 1/(1/13 + 2*1e-6*t_unitless)
        yref = np.zeros((xout.size, 2))
        yref[:, 0] = Aref
        yref[:, 1] = B0 + 3*(A0-Aref)
        assert allclose(yout, yref)

    unique['ramp_rate'] = 2
    x, y, p = odesys.to_arrays(tout, conc, unique)
    fout = odesys.f_cb(x, y, p)

    def r(t):
        return A*np.exp(-Ea_over_R/(T0 + dTdt*t))*A0  # refA(t, *params)

    ref = np.array([
        [-r(2), -r(5), -r(10)],
        [3*r(2), 3*r(5), 3*r(10)]
    ]).T
    assert np.allclose(fout, ref)


@requires('numpy', 'pyodesys')
def test_get_odesys__late_binding():
    def _gibbs(args, T, R, backend, **kwargs):
        H, S = args
        return backend.exp(-(H - T*S)/(R*T))

    def _eyring(args, T, R, k_B, h, backend, **kwargs):
        H, S = args
        return k_B/h*T*backend.exp(-(H - T*S)/(R*T))

    gibbs_pk = ('temperature', 'molar_gas_constant')
    eyring_pk = gibbs_pk + ('Boltzmann_constant', 'Planck_constant')

    GibbsEC = MassActionEq.from_callback(_gibbs, argument_names=('H', 'S'), parameter_keys=gibbs_pk)
    EyringMA = MassAction.from_callback(_eyring, argument_names=('H', 'S'), parameter_keys=eyring_pk)

    uk_equil = ('He_assoc', 'Se_assoc')
    beta = GibbsEC(unique_keys=uk_equil)  # equilibrium parameters

    uk_kinet = ('Ha_assoc', 'Sa_assoc')
    bimol_barrier = EyringMA(unique_keys=uk_kinet)  # activation parameters

    eq = Equilibrium({'Fe+3', 'SCN-'}, {'FeSCN+2'}, beta)
    rsys = ReactionSystem(eq.as_reactions(kf=bimol_barrier))
    odesys, extra = get_odesys(rsys, include_params=False)
    pk, unique, p_units = map(extra.get, 'param_keys unique p_units'.split())
    assert sorted(unique) == sorted(uk_equil + uk_kinet)
    assert sorted(pk) == sorted(eyring_pk)


@requires('numpy', 'pyodesys')
def test_get_odesys__ScaledSys():
    from pyodesys.symbolic import ScaledSys
    k = .2
    a = Substance('A')
    b = Substance('B')
    r = Reaction({'A': 1}, {'B': 1}, param=k)
    rsys = ReactionSystem([r], [a, b])
    assert sorted(rsys.substances.keys()) == ['A', 'B']
    odesys = get_odesys(rsys, include_params=True, SymbolicSys=ScaledSys)[0]
    c0 = {
        'A': 1.0,
        'B': 3.0,
    }
    t = np.linspace(0.0, 10.0)
    xout, yout, info = odesys.integrate(t, c0)
    yref = np.zeros((t.size, 2))
    yref[:, 0] = np.exp(-k*t)
    yref[:, 1] = 4 - np.exp(-k*t)
    assert np.allclose(yout, yref)


@requires('numpy', 'pyodesys', 'sympy')
def test_get_odesys__max_euler_step_cb():
    rsys = ReactionSystem.from_string('\n'.join(['H2O -> H+ + OH-; 1e-4', 'OH- + H+ -> H2O; 1e10']))
    odesys, extra = get_odesys(rsys)
    r1 = 1.01e-4
    r2 = 6e-4
    dH2Odt = r2 - r1
    euler_ref = 2e-7/dH2Odt
    assert abs(extra['max_euler_step_cb'](0, {'H2O': 1.01, 'H+': 2e-7, 'OH-': 3e-7}) - euler_ref)/euler_ref < 1e-8


@requires('numpy', 'pyodesys', 'sympy')
@pytest.mark.parametrize('substances', permutations(['H2O', 'H+', 'OH-']))
def test_get_odesys__linear_dependencies__preferred(substances):
    rsys = ReactionSystem.from_string('\n'.join(['H2O -> H+ + OH-; 1e-4', 'OH- + H+ -> H2O; 1e10']), substances)
    assert isinstance(rsys.substances, OrderedDict)
    odesys, extra = get_odesys(rsys)

    af_H2O_H = extra['linear_dependencies'](['H+', 'H2O'])
    import sympy
    y0 = {k: sympy.Symbol(k+'0') for k in rsys.substances}
    af_H2O_H(None, {odesys[k]: v for k, v in y0.items()}, None, sympy)  # ensure idempotent
    exprs_H2O_H = af_H2O_H(None, {odesys[k]: v for k, v in y0.items()}, None, sympy)
    ref_H2O_H = {
        'H2O': y0['H2O'] + y0['OH-'] - odesys['OH-'],  # oxygen
        'H+': 2*y0['H2O'] + y0['H+'] + y0['OH-'] - odesys['OH-'] - 2*(
            y0['H2O'] + y0['OH-'] - odesys['OH-'])  # hydrogen
    }
    for k, v in ref_H2O_H.items():
        assert (exprs_H2O_H[odesys[k]] - v) == 0


@requires('numpy', 'pyodesys', 'sympy', 'pycvodes')
@pytest.mark.parametrize('preferred', [None, ['H+', 'OH-'], ['H2O', 'H+'], ['H2O', 'OH-']])
def test_get_odesys__linear_dependencies__PartiallySolvedSystem(preferred):
    import sympy
    from pyodesys.symbolic import PartiallySolvedSystem
    rsys = ReactionSystem.from_string('\n'.join(['H2O -> H+ + OH-; 1e-4', 'OH- + H+ -> H2O; 1e10']))
    odesys, extra = get_odesys(rsys)
    c0 = {'H2O': 0, 'H+': 2e-7, 'OH-': 3e-7}
    h0max = extra['max_euler_step_cb'](0, c0)
    analytic_factory = extra['linear_dependencies']()
    y0 = {k: sympy.Symbol(k+'0') for k in rsys.substances}
    analytic_factory(None, {odesys[k]: v for k, v in y0.items()}, None, sympy)
    psys = PartiallySolvedSystem(odesys, analytic_factory)
    xout, yout, info = psys.integrate(1, c0, atol=1e-12, rtol=1e-10, first_step=h0max*1e-12,
                                      integrator='cvode')
    c_reac = c0['H+'], c0['OH-']
    H2O_ref = binary_rev(xout, 1e10, 1e-4, c0['H2O'], max(c_reac), min(c_reac))
    assert np.allclose(yout[:, psys.names.index('H2O')], H2O_ref)
    assert np.allclose(yout[:, psys.names.index('H+')], c0['H+'] + c0['H2O'] - H2O_ref)
    assert np.allclose(yout[:, psys.names.index('OH-')], c0['OH-'] + c0['H2O'] - H2O_ref)


@requires('numpy', 'pyodesys', 'sympy', 'pycvodes')
def test_get_odesys__Equilibrium_as_reactions():
    from chempy import Equilibrium, ReactionSystem
    eq = Equilibrium({'Fe+3', 'SCN-'}, {'FeSCN+2'}, 10**2)
    substances = 'Fe+3 SCN- FeSCN+2'.split()
    rsys = ReactionSystem(eq.as_reactions(kf=3.0), substances)
    odesys, extra = get_odesys(rsys)
    init_conc = {'Fe+3': 1.0, 'SCN-': .3, 'FeSCN+2': 0}
    tout, Cout, info = odesys.integrate(5, init_conc, integrator='cvode', atol=1e-11, rtol=1e-12)
    cmplx_ref = binary_rev(tout, 3, 3.0/100, init_conc['FeSCN+2'], init_conc['Fe+3'], init_conc['SCN-'])
    assert np.allclose(Cout[:, 2], cmplx_ref)


@requires('numpy', 'pyodesys', 'sympy', 'pycvodes')
def test_get_odesys__Expr_as_param():
    def _eyring_pe(args, T, backend=math, **kwargs):
        freq, = args
        return freq*T

    EyringPreExp = Expr.from_callback(_eyring_pe, argument_names=('freq',),
                                      parameter_keys=('temperature',))

    def _k(args, T, backend=math, **kwargs):
        A, H, S = args
        return A*backend.exp(-(H - T*S)/(8.314511*T))

    EyringMA = MassAction.from_callback(_k, parameter_keys=('temperature',),
                                        argument_names=('Aa', 'Ha', 'Sa'))
    kb_h = 2.08e10
    rxn = Reaction({'A'}, {'B'}, EyringMA(unique_keys=('A_u', 'H_u', 'S_u')))
    rsys = ReactionSystem([rxn], ['A', 'B'])
    odesys, extra = get_odesys(rsys, include_params=False, substitutions={'A_u': EyringPreExp(kb_h)})
    y0 = defaultdict(float, {'A': 7.0})
    rt = 293.15
    xout, yout, info = odesys.integrate(5, y0, {'H_u': 117e3, 'S_u': 150, 'temperature': rt},
                                        integrator='cvode', atol=1e-12, rtol=1e-10, nsteps=1000)
    kref = kb_h*rt*np.exp(-(117e3 - rt*150)/(8.314511*rt))
    ref = y0['A']*np.exp(-kref*xout)
    assert np.allclose(yout[:, 0], ref)
    assert np.allclose(yout[:, 1], y0['A'] - ref)


@requires('numpy', 'pyodesys', 'sympy', 'pycvodes')
def test_get_odesys__Expr_as_param__unique_as_param():
    def _eyring_pe_coupled(args, T, S, backend=math, **kwargs):
        freq, = args
        return freq*T/S

    EyringPreExpCoupled = Expr.from_callback(_eyring_pe_coupled, argument_names=('freq',),
                                             parameter_keys=('temperature', 'S_u'))

    def _k(args, T, backend=math, **kwargs):
        A, H, S = args
        return A*backend.exp(-(H - T*S)/(8.314511*T))

    EyringMA = MassAction.from_callback(_k, parameter_keys=('temperature',),
                                        argument_names=('Aa', 'Ha', 'Sa'))
    kb_h = 2.08e10
    rxn = Reaction({'A'}, {'B'}, EyringMA(unique_keys=('A_u', 'H_u', 'S_u')))
    rsys = ReactionSystem([rxn], ['A', 'B'])
    odesys2, extra2 = get_odesys(rsys, include_params=False,
                                 substitutions={'A_u': EyringPreExpCoupled(kb_h)})
    y0 = defaultdict(float, {'A': 7.0})
    rt = 293.15
    xout2, yout2, info2 = odesys2.integrate(5, y0, {'H_u': 107e3, 'S_u': 150, 'temperature': rt},
                                            integrator='cvode', atol=1e-12, rtol=1e-10, nsteps=1000)
    kref2 = kb_h*rt*np.exp(-(107e3 - rt*150)/(8.314511*rt))/150
    ref2 = y0['A']*np.exp(-kref2*xout2)
    assert np.allclose(yout2[:, 0], ref2)
    assert np.allclose(yout2[:, 1], y0['A'] - ref2)


@requires('pyodesys', 'pycvodes')
def test_chained_parameter_variation():
    ratex = MassAction(Arrhenius([1e10, 63e3/8.3145]))
    rxn = Reaction({'A': 1}, {'B': 1}, ratex)
    rsys = ReactionSystem([rxn], 'A B')
    odesys, extra = get_odesys(rsys, include_params=False)
    param_keys, unique_keys, p_units = map(extra.get, 'param_keys unique p_units'.split())
    conc = {'A': 3.17, 'B': 5.03}
    Ts = (294, 304, 317)
    times = [3.1, 2.1, 5.3]
    kw = dict(integrator='cvode', atol=1e-12, rtol=1e-13, first_step=1e-14)
    tout, cout, info = chained_parameter_variation(
        odesys, times, conc, {'temperature': Ts}, {}, integrate_kwargs=kw)
    assert len(info['nfev']) == 3
    assert info['nfev'][0] > 2
    assert info['nfev'][1] > 2
    assert info['nfev'][2] > 2
    assert np.all(np.diff(tout) > 0)
    tout1 = tout[tout <= times[0]]
    tout23 = tout[tout > times[0]]
    tout2 = tout23[tout23 <= times[0] + times[1]]
    tout3 = tout23[tout23 > times[0] + times[1]]

    def _ref(y0, x, T, x0):
        k = 1e10*np.exp(-63e3/8.3145/T)
        return y0*np.exp(-k*(x-x0))

    Aref1 = _ref(conc['A'], tout1, Ts[0], tout1[0])
    Bref1 = conc['B'] + conc['A'] - Aref1

    Aref2 = _ref(Aref1[-1], tout2, Ts[1], tout1[-1])
    Bref2 = Bref1[-1] + Aref1[-1] - Aref2

    Aref3 = _ref(Aref2[-1], tout3, Ts[2], tout2[-1])
    Bref3 = Bref2[-1] + Aref2[-1] - Aref3

    cref = np.concatenate([np.vstack((a, b)).T for a, b in [(Aref1, Bref1), (Aref2, Bref2), (Aref3, Bref3)]])
    forgive = 27*1.1
    assert np.allclose(cref, cout, atol=kw['atol']*forgive, rtol=kw['rtol']*forgive)


@requires('pyodesys', 'scipy', 'sym')
def test_get_odesys__cstr():
    rsys = ReactionSystem.from_string("2 H2O2 -> O2 + 2 H2O; 5")
    odesys, extra = get_odesys(rsys, cstr=True)
    fr, fc = extra['cstr_fr_fc']
    tout, c0 = np.linspace(0, .13, 7), {'H2O2': 2, 'O2': 4, 'H2O': 3}
    params = {fr: 13, fc['H2O2']: 11, fc['O2']: 43, fc['H2O']: 45}
    res = odesys.integrate(tout, c0, params)
    from chempy.kinetics.integrated import binary_irrev_cstr

    def get_analytic(result, k, n):
        ref = binary_irrev_cstr(result.xout, 5, result.named_dep('H2O2')[0],
                                result.named_dep(k)[0], result.named_param(fc['H2O2']), result.named_param(fc[k]),
                                result.named_param(fr), n)
        return np.array(ref).T

    ref_O2 = get_analytic(res, 'O2', 1)
    ref_H2O = get_analytic(res, 'H2O', 2)
    assert np.allclose(res.named_dep('H2O2'), ref_O2[:, 0])
    assert np.allclose(res.named_dep('H2O2'), ref_H2O[:, 0])
    assert np.allclose(res.named_dep('O2'), ref_O2[:, 1])
    assert np.allclose(res.named_dep('H2O'), ref_H2O[:, 1])


@requires('pygslodeiv2', 'sym', units_library)
def test_get_odesys_rsys_with_units():
    rsys = ReactionSystem.from_string("""
    A -> B; 0.096/s
    B + C -> P; 4e3/M/s
    """, substance_factory=Substance)
    with pytest.raises(Exception):
        get_odesys(rsys)  # not a strict test, SI_base_registry could be made the default

    odesys, extra = get_odesys(rsys, unit_registry=SI_base_registry)
    tend = 10
    tend_units = tend*u.s
    c0 = {'A': 1e-6, 'B': 0, 'C': 1, 'P': 0}
    c0_units = {k: v*u.molar for k, v in c0.items()}
    result1 = odesys.integrate(tend_units, c0_units, integrator='gsl')
    assert result1.info['success']

    with pytest.raises(Exception):
        odesys.integrate(tend, c0, integrator='gsl')


@requires('pyodeint', 'sym', units_library)
def test_get_odesys_rsys_with_units__named_params():
    rsys = ReactionSystem.from_string("""
    A -> B; 'k1'
    B + C -> P; 'k2'
    """, substance_factory=Substance)
    odesys, extra = get_odesys(rsys, include_params=False, unit_registry=SI_base_registry)
    tend = 10
    tend_units = tend*u.s
    c0 = {'A': 1e-6, 'B': 0, 'C': 1, 'P': 0}
    p = {'k1': 3, 'k2': 4}
    p_units = {'k1': 3/u.s, 'k2': 4/u.M/u.s}
    c0_units = {k: v*u.molar for k, v in c0.items()}
    result1 = odesys.integrate(tend_units, c0_units, p_units, integrator='odeint')
    assert result1.info['success']

    with pytest.raises(Exception):
        odesys.integrate(tend, c0, p, integrator='odeint')


@requires('pycvodes', 'sym', units_library)
def test_get_odesys__Eyring():
    R = 8.314472
    T_K = 300
    dH = 80e3
    dS = 10
    rsys1 = ReactionSystem.from_string("""
    NOBr -> NO + Br; EyringParam(dH={dH}*J/mol, dS={dS}*J/K/mol)
    """.format(dH=dH, dS=dS), substances='NOBr NO Br'.split())
    kref = 20836643994.118652*T_K*np.exp(-(dH - T_K*dS)/(R*T_K))
    NOBr0_M = 0.7
    init_cond = dict(
        NOBr=NOBr0_M*u.M,
        NO=0*u.M,
        Br=0*u.M
    )
    t = 5*u.second
    params = dict(
        temperature=T_K*u.K
    )

    def check(rsys):
        odesys, extra = get_odesys(rsys, unit_registry=SI_base_registry, constants=const)
        res = odesys.integrate(t, init_cond, params, integrator='cvode')
        NOBr_ref = NOBr0_M*np.exp(-kref*to_unitless(res.xout, u.second))
        ref = np.array([NOBr_ref] + [NOBr0_M - NOBr_ref]*2).T
        cmp = to_unitless(res.yout, u.M)
        assert np.allclose(cmp, ref)
    check(rsys1)
    rsys2 = ReactionSystem.from_string("""
    NOBr -> NO + Br; MassAction(EyringHS([{dH}*J/mol, {dS}*J/K/mol]))
    """.format(dH=dH, dS=dS), substances='NOBr NO Br'.split())
    check(rsys2)


@requires('pycvodes', 'sym', units_library)
def test_get_odesys__Eyring_2nd_order():
    R = 8.314472
    T_K = 300
    dH = 80e3
    dS = 10
    rsys1b = ReactionSystem.from_string("""
    NO + Br -> NOBr; EyringParam(dH={dH}*J/mol, dS={dS}*J/K/mol)
    """.format(dH=dH, dS=dS))
    c0 = 1  # mol/dm3 === 1000 mol/m3
    kbref = 20836643994.118652*T_K*np.exp(-(dH - T_K*dS)/(R*T_K))/c0
    NO0_M = 1.5
    Br0_M = 0.7
    init_cond = dict(
        NOBr=0*u.M,
        NO=NO0_M*u.M,
        Br=Br0_M*u.M
    )
    t = 5*u.second
    params = dict(
        temperature=T_K*u.K
    )

    def analytic_b(t):
        U, V = NO0_M, Br0_M
        d = U - V
        return (U*(1 - np.exp(-kbref*t*d)))/(U/V - np.exp(-kbref*t*d))

    def check(rsys):
        odesys, extra = get_odesys(rsys, unit_registry=SI_base_registry, constants=const)
        res = odesys.integrate(t, init_cond, params, integrator='cvode')
        t_sec = to_unitless(res.xout, u.second)
        NOBr_ref = analytic_b(t_sec)
        cmp = to_unitless(res.yout, u.M)
        ref = np.empty_like(cmp)
        ref[:, odesys.names.index('NOBr')] = NOBr_ref
        ref[:, odesys.names.index('Br')] = Br0_M - NOBr_ref
        ref[:, odesys.names.index('NO')] = NO0_M - NOBr_ref
        assert np.allclose(cmp, ref)

    check(rsys1b)
    rsys2b = ReactionSystem.from_string("""
    NO + Br -> NOBr; MassAction(EyringHS([{dH}*J/mol, {dS}*J/K/mol]))
    """.format(dH=dH, dS=dS))
    check(rsys2b)


@requires('pycvodes', 'sym', 'scipy', units_library)
def test_get_odesys__Eyring_1st_order_linearly_ramped_temperature():
    from scipy.special import expi

    def analytic_unit0(t, T0, dH, dS):
        R = 8.314472
        kB = 1.3806504e-23
        h = 6.62606896e-34
        A = kB/h*np.exp(dS/R)
        B = dH/R
        return np.exp(A*(
            (-B**2*np.exp(B/T0)*expi(-B/T0) - T0*(B - T0))*np.exp(-B/T0) +
            (B**2*np.exp(B/(t + T0))*expi(-B/(t + T0)) - (t + T0)*(-B + t + T0))*np.exp(-B/(t + T0))
        )/2)

    T_K = 290
    dH = 80e3
    dS = 10
    rsys1 = ReactionSystem.from_string("""
    NOBr -> NO + Br; EyringParam(dH={dH}*J/mol, dS={dS}*J/K/mol)
    """.format(dH=dH, dS=dS))

    NOBr0_M = 0.7
    init_cond = dict(
        NOBr=NOBr0_M*u.M,
        NO=0*u.M,
        Br=0*u.M
    )
    t = 20*u.second

    def check(rsys):
        odes, extra = get_odesys(rsys, unit_registry=SI_base_registry, constants=const, substitutions={
            'temperature': RampedTemp([T_K*u.K, 1*u.K/u.s])})
        for odesys in [odes, odes.as_autonomous()]:
            res = odesys.integrate(t, init_cond, integrator='cvode')
            t_sec = to_unitless(res.xout, u.second)
            NOBr_ref = NOBr0_M*analytic_unit0(t_sec, T_K, dH, dS)
            cmp = to_unitless(res.yout, u.M)
            ref = np.empty_like(cmp)
            ref[:, odesys.names.index('NOBr')] = NOBr_ref
            ref[:, odesys.names.index('Br')] = NOBr0_M - NOBr_ref
            ref[:, odesys.names.index('NO')] = NOBr0_M - NOBr_ref
            assert np.allclose(cmp, ref)

    check(rsys1)
    rsys2 = ReactionSystem.from_string("""
    NOBr -> NO + Br; MassAction(EyringHS([{dH}*J/mol, {dS}*J/K/mol]))
    """.format(dH=dH, dS=dS))
    check(rsys2)


@requires('pycvodes', 'sym', 'scipy', units_library)
def test_get_odesys__Eyring_2nd_order_linearly_ramped_temperature():
    from scipy.special import expi

    def analytic_unit0(t, k, m, dH, dS):
        R = 8.314472
        kB = 1.3806504e-23
        h = 6.62606896e-34
        A = kB/h*np.exp(dS/R)
        B = dH/R
        return k*np.exp(B*(k*t + 2*m)/(m*(k*t + m)))/(
            A*(-B**2*np.exp(B/(k*t + m))*expi(-B/(k*t + m)) - B*k*t - B*m + k**2*t**2 + 2*k*m*t + m**2)*np.exp(B/m) +
            (A*B**2*np.exp(B/m)*expi(-B/m) - A*m*(-B + m) + k*np.exp(B/m))*np.exp(B/(k*t + m))
        )

    T_K = 290
    dTdt_Ks = 3
    dH = 80e3
    dS = 10
    rsys1 = ReactionSystem.from_string("""
    2 NO2 -> N2O4; EyringParam(dH={dH}*J/mol, dS={dS}*J/K/mol)
    """.format(dH=dH, dS=dS))

    NO2_M = 1.0
    init_cond = dict(
        NO2=NO2_M*u.M,
        N2O4=0*u.M
    )
    t = 20*u.second

    def check(rsys):
        odes, extra = get_odesys(rsys, unit_registry=SI_base_registry, constants=const, substitutions={
            'temperature': RampedTemp([T_K*u.K, dTdt_Ks*u.K/u.s])})
        for odesys in [odes, odes.as_autonomous()]:
            res = odesys.integrate(t, init_cond, integrator='cvode')
            t_sec = to_unitless(res.xout, u.second)
            NO2_ref = analytic_unit0(t_sec, dTdt_Ks, T_K, dH, dS)
            cmp = to_unitless(res.yout, u.M)
            ref = np.empty_like(cmp)
            ref[:, odesys.names.index('NO2')] = NO2_ref
            ref[:, odesys.names.index('N2O4')] = (NO2_M - NO2_ref)/2
            assert np.allclose(cmp, ref)

    check(rsys1)
    rsys2 = ReactionSystem.from_string("""
    2 NO2 -> N2O4; MassAction(EyringHS([{dH}*J/mol, {dS}*J/K/mol]))
    """.format(dH=dH, dS=dS))
    check(rsys2)


@requires('pycvodes', 'sym', units_library)
def test_get_odesys__Eyring_2nd_order_reversible():
    R = 8.314472
    T_K = 273.15 + 20  # 20 degree celsius
    kB = 1.3806504e-23
    h = 6.62606896e-34

    dHf = 74e3
    dSf = R*np.log(h/kB/T_K*1e16)

    dHb = 79e3
    dSb = dSf - 23

    rsys1 = ReactionSystem.from_string("""
    Fe+3 + SCN- -> FeSCN+2; EyringParam(dH={dHf}*J/mol, dS={dSf}*J/K/mol)
    FeSCN+2 -> Fe+3 + SCN-; EyringParam(dH={dHb}*J/mol, dS={dSb}*J/K/mol)
    """.format(dHf=dHf, dSf=dSf, dHb=dHb, dSb=dSb))
    kf_ref = 20836643994.118652*T_K*np.exp(-(dHf - T_K*dSf)/(R*T_K))
    kb_ref = 20836643994.118652*T_K*np.exp(-(dHb - T_K*dSb)/(R*T_K))
    Fe0 = 6e-3
    SCN0 = 2e-3
    init_cond = {
        'Fe+3': Fe0*u.M,
        'SCN-': SCN0*u.M,
        'FeSCN+2': 0*u.M
    }
    t = 3*u.second

    def check(rsys, params):
        odes, extra = get_odesys(rsys, include_params=False,
                                 unit_registry=SI_base_registry, constants=const)
        for odesys in [odes, odes.as_autonomous()]:
            res = odesys.integrate(t, init_cond, params, integrator='cvode')
            t_sec = to_unitless(res.xout, u.second)
            FeSCN_ref = binary_rev(t_sec, kf_ref, kb_ref, 0, Fe0, SCN0)
            cmp = to_unitless(res.yout, u.M)
            ref = np.empty_like(cmp)
            ref[:, odesys.names.index('FeSCN+2')] = FeSCN_ref
            ref[:, odesys.names.index('Fe+3')] = Fe0 - FeSCN_ref
            ref[:, odesys.names.index('SCN-')] = SCN0 - FeSCN_ref
            assert np.allclose(cmp, ref)

    check(rsys1, {'temperature': T_K*u.K})
    rsys2 = ReactionSystem.from_string("""
    Fe+3 + SCN- -> FeSCN+2; MassAction(EyringHS([{dHf}*J/mol, {dSf}*J/K/mol]))
    FeSCN+2 -> Fe+3 + SCN-; MassAction(EyringHS([{dHb}*J/mol, {dSb}*J/K/mol]))
    """.format(dHf=dHf, dSf=dSf, dHb=dHb, dSb=dSb))
    check(rsys2, {'temperature': T_K*u.K})
    rsys3 = ReactionSystem.from_string("""
    Fe+3 + SCN- -> FeSCN+2; MassAction(EyringHS.fk('dHf', 'dSf'))
    FeSCN+2 -> Fe+3 + SCN-; MassAction(EyringHS.fk('dHb', 'dSb'))
    """)
    check(rsys3, dict(temperature=T_K*u.K,
                      dHf=dHf*u.J/u.mol, dSf=dSf*u.J/u.mol/u.K,
                      dHb=dHb*u.J/u.mol, dSb=dSb*u.J/u.mol/u.K))


@requires('numpy', 'pyodesys', 'sympy', 'pycvodes')
def test_create_odesys():
    rsys = ReactionSystem.from_string("""
    A -> B; 'k1'
    B + C -> P; 'k2'
    """, substance_factory=Substance)

    odesys, odesys_extra = create_odesys(rsys, unit_registry=SI_base_registry)

    tend = 10
    c0 = {'A': 1e-6, 'B': 0, 'C': 1, 'P': 0}
    p = dict(k1=3, k2=4)

    tend_units = tend*u.s
    p_units = {'k1': p['k1']/u.s, 'k2': p['k2']/u.M/u.s}
    c0_units = {k: v*u.molar for k, v in c0.items()}

    validation = odesys_extra['validate'](dict(c0_units, **p_units))
    P, = validation['not_seen']
    assert P.name == 'P'

    ref_rates = {'A': -p_units['k1']*c0_units['A'], 'P': p_units['k2']*c0_units['B']*c0_units['C']}
    ref_rates['B'] = -ref_rates['A'] - ref_rates['P']
    ref_rates['C'] = -ref_rates['P']
    assert validation['rates'] == ref_rates

    result1, result1_extra = odesys_extra['unit_aware_solve'](
        tend_units, c0_units, p_units, integrator='cvode')
    assert result1.info['success']

    result2 = odesys.integrate(tend, c0, p, integrator='cvode')
    assert np.allclose(result2.yout[-1, :], to_unitless(result1.yout[-1, :], u.molar))


@requires('pycvodes', 'sym', units_library)
def test_create_odesys__Radiolytic():
    rsys1 = ReactionSystem.from_string("""
    -> e-(aq); Radiolytic.fk('g_emaq')
    """, checks=())
    ic1 = {'e-(aq)': 0.0}
    t1 = 5
    p1 = dict(
        g_emaq=42.0,
        doserate=17.0,
        density=5.0
    )
    odesys1, odesys_extra = create_odesys(rsys1)
    result1 = odesys1.integrate(t1, ic1, p1)
    yref1 = result1.xout*p1['g_emaq']*p1['doserate']*p1['density']
    assert np.allclose(yref1, result1.yout.squeeze())
