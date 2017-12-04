# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from collections import defaultdict
from functools import reduce
from operator import mul

try:
    import numpy as np
except ImportError:
    np = None

import pytest

from chempy import ReactionSystem
from chempy.units import (
    allclose, logspace_from_lin, unitless_in_registry, to_unitless, SI_base_registry,
    rescale, default_units as u
)
from chempy.util.testing import requires
from ..integrated import binary_rev
from ..ode import get_odesys
from .._native import get_native


decay_analytic = {
    0: lambda y0, k, t: (
        y0[0] * np.exp(-k[0]*t)),
    1: lambda y0, k, t: (
        y0[1] * np.exp(-k[1] * t) + y0[0] * k[0] / (k[1] - k[0]) *
        (np.exp(-k[0]*t) - np.exp(-k[1]*t))),
    2: lambda y0, k, t: (
        y0[2] * np.exp(-k[2] * t) + y0[1] * k[1] / (k[2] - k[1]) *
        (np.exp(-k[1]*t) - np.exp(-k[2]*t)) +
        k[1] * k[0] * y0[0] / (k[1] - k[0]) *
        (1 / (k[2] - k[0]) * (np.exp(-k[0]*t) - np.exp(-k[2]*t)) -
         1 / (k[2] - k[1]) * (np.exp(-k[1]*t) - np.exp(-k[2]*t))))
}


def decay_get_Cref(k, y0, tout):
    coeffs = list(k) + [0]*(3-len(k))

    return np.column_stack([
        decay_analytic[i](y0, coeffs, tout) for i in range(
            min(3, len(k)+1))])


@pytest.mark.veryslow
@requires('pycvodes', 'pyodesys')
@pytest.mark.parametrize('solve', [(), ('CNO',)])
def test_get_native__first_step(solve):
    integrator = 'cvode'
    forgive = 20

    def k(num):
        return "MassAction(unique_keys=('k%d',))" % num

    lines = [  # fictitious isomerization
        "CNO -> ONC; %s" % k(1),
        "ONC -> NCO; %s" % k(2),
        "NCO -> CON; %s" % k(3)
    ]
    rsys = ReactionSystem.from_string('\n'.join(lines), 'CNO ONC NCO CON')
    odesys, extra = get_odesys(rsys, include_params=False)
    if len(solve) > 0:
        from pyodesys.symbolic import PartiallySolvedSystem
        odesys = PartiallySolvedSystem(odesys, extra['linear_dependencies'](solve))
    c0 = defaultdict(float, {'CNO': .7})
    rate_coeffs = (1e78, 2, 3.)
    args = (5, c0, dict(zip('k1 k2 k3'.split(), rate_coeffs)))
    kwargs = dict(integrator=integrator, atol=1e-9, rtol=1e-9, nsteps=1000)
    native = get_native(rsys, odesys, integrator)

    h0 = extra['max_euler_step_cb'](0, *args[1:])
    xout1, yout1, info1 = odesys.integrate(*args, first_step=h0, **kwargs)
    xout2, yout2, info2 = native.integrate(*args, **kwargs)
    ref1 = decay_get_Cref(rate_coeffs, [c0[key] for key in native.names], xout1)
    ref2 = decay_get_Cref(rate_coeffs, [c0[key] for key in native.names], xout2)
    allclose_kw = dict(atol=kwargs['atol']*forgive, rtol=kwargs['rtol']*forgive)

    assert np.allclose(yout1[:, :3], ref1, **allclose_kw)

    assert info2['success']
    assert info2['nfev'] > 10 and info2['nfev'] > 1 and info2['time_cpu'] < 10 and info2['time_wall'] < 10
    assert np.allclose(yout2[:, :3], ref2, **allclose_kw)


@pytest.mark.veryslow
@requires('pygslodeiv2', 'pyodesys')
@pytest.mark.parametrize('solve', [(), ('H2O',)])
def test_get_native__a_substance_no_composition(solve):
    rsys = ReactionSystem.from_string('\n'.join(['H2O -> H2O+ + e-(aq); 1e-8', 'e-(aq) + H2O+ -> H2O; 1e10']))
    odesys, extra = get_odesys(rsys)
    c0 = {'H2O': 0, 'H2O+': 2e-9, 'e-(aq)': 3e-9}
    if len(solve) > 0:
        from pyodesys.symbolic import PartiallySolvedSystem
        odesys = PartiallySolvedSystem(odesys, extra['linear_dependencies'](solve))
    odesys = get_native(rsys, odesys, 'gsl')
    xout, yout, info = odesys.integrate(1, c0, atol=1e-15, rtol=1e-15, integrator='gsl')
    c_reac = c0['H2O+'], c0['e-(aq)']
    H2O_ref = binary_rev(xout, 1e10, 1e-4, c0['H2O'], max(c_reac), min(c_reac))
    assert np.allclose(yout[:, odesys.names.index('H2O')], H2O_ref)
    assert np.allclose(yout[:, odesys.names.index('H2O+')], c0['H2O+'] + c0['H2O'] - H2O_ref)
    assert np.allclose(yout[:, odesys.names.index('e-(aq)')], c0['e-(aq)'] + c0['H2O'] - H2O_ref)


@pytest.mark.veryslow
@requires('pycvodes', 'pyodesys')
@pytest.mark.parametrize('dep_scaling', [1, 763])
def test_get_native__named_parameter__units(dep_scaling):
    rsys = ReactionSystem.from_string("""
    -> H; 'p'
    H + H -> H2; 'k2'
    """, checks=('substance_keys', 'duplicate', 'duplicate_names'))
    from pyodesys.symbolic import ScaledSys
    kwargs = {} if dep_scaling == 1 else dict(SymbolicSys=ScaledSys, dep_scaling=dep_scaling)
    odesys, extra = get_odesys(rsys, include_params=False, unit_registry=SI_base_registry, **kwargs)
    c0 = {'H': 42e-6*u.molar, 'H2': 17*u.micromolar}
    native = get_native(rsys, odesys, 'cvode')
    tend = 7*60*u.minute
    g_rho_Ddot = g, rho, Ddot = 2*u.per100eV, 998*u.g/u.dm3, 314*u.Gy/u.hour
    params = {
        'p': reduce(mul, g_rho_Ddot),
        'k2': 53/u.molar/u.minute
    }
    result = native.integrate(tend, c0, params, atol=1e-15, rtol=1e-15, integrator='cvode', nsteps=16000)
    assert result.info['success']

    def analytic_H(t, p, k, H0):
        # dH/dt = p - k2*H**2
        x0 = np.sqrt(2)*np.sqrt(p)
        x1 = x0
        x2 = np.sqrt(k)
        x3 = t*x1*x2
        x4 = H0*x2
        x5 = np.sqrt(x0 + 2*x4)
        x6 = np.sqrt(-1/(2*H0*x2 - x0))
        x7 = x5*x6*np.exp(x3)
        x8 = np.exp(-x3)/(x5*x6)
        return x1*(x7 - x8)/(2*x2*(x7 + x8))

    t_ul = to_unitless(result.xout, u.s)
    p_ul = to_unitless(params['p'], u.micromolar/u.s)
    ref_H_uM = analytic_H(
        t_ul,
        p_ul,
        to_unitless(params['k2'], 1/u.micromolar/u.s),
        to_unitless(c0['H'], u.micromolar)
    )
    ref_H2_uM = to_unitless(c0['H2'], u.micromolar) + to_unitless(c0['H'], u.micromolar)/2 + t_ul*p_ul/2 - ref_H_uM/2
    assert np.allclose(to_unitless(result.named_dep('H'), u.micromolar), ref_H_uM)
    assert np.allclose(to_unitless(result.named_dep('H2'), u.micromolar), ref_H2_uM)


@pytest.mark.veryslow
@requires('pycvodes', 'pyodesys')
def test_get_native__conc_roots():
    M, s = u.M, u.s
    rsys = ReactionSystem.from_string("2 O3 -> 3 O2; 'k2'")
    u_reg = SI_base_registry.copy()
    odesys, extra = get_odesys(rsys, include_params=False, unit_registry=u_reg)
    c0 = {'O3': 4.2e-3*M, 'O2': 0*M}
    cr = ['O2']
    native = get_native(rsys, odesys, 'cvode', conc_roots=cr)
    tend = 1e5*u.s
    params = {'k2': logspace_from_lin(1e-3/M/s, 1e3/M/s, 14)}
    tgt_O2 = 1e-3*M
    results = native.integrate(tend, c0, params, integrator='native', return_on_root=True,
                               special_settings=[unitless_in_registry(tgt_O2, u_reg)])
    assert len(results) == params['k2'].size
    # dydt = -p*y**2
    # 1/y0 - 1/y = -2*pt
    # t = 1/2/p*(1/y - 1/y0)
    tgt_O3 = c0['O3'] - 2/3 * tgt_O2
    for r in results:
        ref = rescale(1/2/r.named_param('k2')*(1/tgt_O3 - 1/c0['O3']), u.s)
        assert allclose(r.xout[-1], ref, rtol=1e-6)


@pytest.mark.veryslow
@requires('pycvodes', 'pyodesys')
@pytest.mark.parametrize('scaling_density', [(1, False), (763, False), (1, True)])
def test_get_native__Radiolytic__named_parameter__units(scaling_density):
    scaling, density = scaling_density
    rsys = ReactionSystem.from_string("""
    -> H; Radiolytic(2*per100eV)
    H + H -> H2; 'k2'
    """, checks=('substance_keys', 'duplicate', 'duplicate_names'))
    gval = 2*u.per100eV

    from pyodesys.symbolic import ScaledSys
    kwargs = {} if scaling == 1 else dict(SymbolicSys=ScaledSys, dep_scaling=scaling)
    dens = {'density': 998*u.g/u.dm3}
    odesys, extra = get_odesys(rsys, include_params=False, substitutions=dens if density else {},
                               unit_registry=SI_base_registry, **kwargs)
    c0 = {'H': 42e-6*u.molar, 'H2': 17e3*u.nanomolar}
    native = get_native(rsys, odesys, 'cvode')
    tend = 7*60*u.minute
    params = {'doserate': 314*u.Gy/u.hour, 'k2': 53/u.molar/u.minute}
    if not density:
        params.update(dens)
    result = native.integrate(tend, c0, params, atol=1e-15, rtol=1e-15, integrator='cvode', nsteps=8000)
    assert result.info['success']

    def analytic_H(t, p, k, H0):
        # dH/dt = p - k2*H**2
        x0 = np.sqrt(2)*np.sqrt(p)
        x1 = x0
        x2 = np.sqrt(k)
        x3 = t*x1*x2
        x4 = H0*x2
        x5 = np.sqrt(x0 + 2*x4)
        x6 = np.sqrt(-1/(2*H0*x2 - x0))
        x7 = x5*x6*np.exp(x3)
        x8 = np.exp(-x3)/(x5*x6)
        return x1*(x7 - x8)/(2*x2*(x7 + x8))

    t_ul = to_unitless(result.xout, u.s)
    p_ul = to_unitless(params['doserate']*dens['density']*gval, u.micromolar/u.s)
    ref_H_uM = analytic_H(
        t_ul,
        p_ul,
        to_unitless(params['k2'], 1/u.micromolar/u.s),
        to_unitless(c0['H'], u.micromolar)
    )
    ref_H2_uM = to_unitless(c0['H2'], u.micromolar) + to_unitless(c0['H'], u.micromolar)/2 + t_ul*p_ul/2 - ref_H_uM/2
    assert np.allclose(to_unitless(result.named_dep('H'), u.micromolar), ref_H_uM)
    assert np.allclose(to_unitless(result.named_dep('H2'), u.micromolar), ref_H2_uM)
