# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math

import pytest

from chempy import Reaction, ReactionSystem, Substance
from chempy.units import (
    allclose, Backend, to_unitless, units_library, default_constants, default_units as u
)
from chempy.util._expr import Expr
from chempy.util.parsing import parsing_library
from chempy.util.testing import requires
from ..rates import RateExpr, MassAction, Arrhenius, Radiolytic, mk_Radiolytic, Eyring


class SpecialFraction(RateExpr):
    """
    Example from Atkins, De Paula, Physical Chemistry
    H2 + Br2 -> 2HBr
    ½ dHBr/dt = k[H2][Br2]**(3/2) / ([Br2] + k'[HBr])
    """
    def __call__(self, variables, backend=math, reaction=None):
        two = 2 * backend.pi**0
        k = self.arg(variables, 0)
        kp = self.arg(variables, 1)
        H2, Br2, HBr = variables['H2'], variables['Br2'], variables['HBr']
        return k*H2*Br2**(3/two) / (Br2 + kp*HBr)


def _get_SpecialFraction_rsys(k, kprime):
    ratex = SpecialFraction([k, kprime], ['k_HBr', 'kprime_HBr'])
    rxn = Reaction({'H2': 1, 'Br2': 1}, {'HBr': 2}, ratex)
    return ReactionSystem([rxn], 'H2 Br2 HBr')


@requires('numpy')
def test_specialfraction_rateexpr():
    rsys = _get_SpecialFraction_rsys(11, 13)
    r = rsys.rxns[0]

    conc = {'H2': 2, 'Br2': 3, 'HBr': 5}

    def _check(k, kprime, c):
        ref = k*c['H2']*c['Br2']**1.5/(c['Br2'] + kprime*c['HBr'])
        assert abs(r.rate_expr()(c) - ref) < 1e-15

    _check(11, 13, conc)


def test_RateExpr__subclass_from_callback():
    SF = RateExpr.subclass_from_callback(
        lambda v, a, backend: a[0]*v['H2']*v['Br2']**(3/2) / (v['Br2'] + a[1]*v['HBr'])
    )
    ratex = SF([11, 13], ['k_HBr', 'kprime_HBr'])
    r = Reaction({'H2': 1, 'Br2': 1}, {'HBr': 2}, ratex)

    res1 = r.rate_expr()({'H2': 5, 'Br2': 7, 'HBr': 15})
    ref1 = 11*5*7**1.5/(7+13*15)
    assert abs((res1-ref1)/ref1) < 1e-14

    res2 = r.rate_expr()({'H2': 5, 'Br2': 7, 'HBr': 15, 'k_HBr': 23, 'kprime_HBr': 42})
    ref2 = 23*5*7**1.5/(7+42*15)
    assert abs((res2-ref2)/ref2) < 1e-14


def test_MassAction():
    ma = MassAction([3.14], ['my_rate'])
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, ma, {'B': 1})
    arg1 = {'A': 11, 'B': 13, 'C': 17}
    arg2 = {'A': 11, 'B': 13, 'C': 17, 'my_rate': 2.72}
    res1 = r.rate_expr()(arg1, reaction=r)
    res2 = r.rate_expr()(arg2, reaction=r)
    ref1 = 3.14*13*11**2
    ref2 = 2.72*13*11**2
    assert abs(res1 - ref1) < 1e-14
    assert abs(res2 - ref2) < 1e-12
    rat1 = r.rate(arg1)
    rat2 = r.rate(arg2)
    for key, coeff in [('A', -2), ('B', -2), ('C', 1)]:
        assert abs(rat1[key] - ref1*coeff) < 2e-12
        assert abs(rat2[key] - ref2*coeff) < 2e-12


def test_Expr__from_callback():
    def rate_coeff(args, T, reaction, backend=math):
        return args[0]*backend.exp(args[1]/T)
    RateExpr = Expr.from_callback(rate_coeff, parameter_keys=('Tem',), nargs=2)
    k1 = RateExpr([2.1e10, -5132.2])
    rxn = Reaction({'H2': 2, 'O2': 1}, {'H2O': 2}, MassAction(k1))
    ma = rxn.rate_expr()
    res = ma({'Tem': 273.15, 'H2': 7, 'O2': 13}, reaction=rxn)
    ref = 7*7*13*2.1e10*math.exp(-5132.2/273.15)
    assert abs((res-ref)/ref) < 1e-14


def test_MassAction__subclass_from_callback():
    def rate_coeff(variables, all_args, reaction, backend):
        return all_args[0]*backend.exp(all_args[1]/variables['temperature'])
    CustomMassAction = MassAction.subclass_from_callback(
        rate_coeff, cls_attrs=dict(parameter_keys=('temperature',), nargs=2))
    k1 = CustomMassAction([2.1e10, -5132.2])
    rxn = Reaction({'H2': 2, 'O2': 1}, {'H2O': 2}, k1)
    cma = rxn.rate_expr()
    res = cma({'temperature': 273.15, 'H2': 7, 'O2': 13}, reaction=rxn)
    ref = 7*7*13*2.1e10*math.exp(-5132.2/273.15)
    assert abs((res-ref)/ref) < 1e-14


@requires(units_library)
def test_MassAction__subclass_from_callback__units():
    def rate_coeff(variables, all_args, backend, **kwargs):
        return all_args[0]*backend.exp(all_args[1]/variables['temperature'])
    CustomMassAction = MassAction.subclass_from_callback(
        rate_coeff, cls_attrs=dict(parameter_keys=('temperature',), nargs=2))
    k1 = CustomMassAction([2.1e10/u.molar**2/u.second, -5132.2*u.kelvin])
    rxn = Reaction({'H2': 2, 'O2': 1}, {'H2O': 2}, k1)
    variables = {
        'temperature': 491.67*u.rankine,
        'H2': 7000*u.mol/u.metre**3,
        'O2': 13*u.molar
    }
    cma = rxn.rate_expr()
    res = cma(variables, backend=Backend(), reaction=rxn)
    ref = 7*7*13*2.1e10*math.exp(-5132.2/273.15) * u.molar/u.second
    assert allclose(res, ref)


def test_ArrheniusMassAction():
    A, Ea_over_R = 1.2e11, 40e3/8.3145
    ama = MassAction(Arrhenius([A, Ea_over_R]))
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, ama, {'B': 1})
    T_ = 'temperature'

    def ref(v):
        return 1.2e11*math.exp(-Ea_over_R/v[T_])*v['B']*v['A']**2

    ma = r.rate_expr()
    for params in [(11., 13., 17., 311.2),
                   (12, 8, 5, 270)]:
        var = dict(zip(['A', 'B', 'C', T_], params))
        ref_var = ref(var)
        assert abs((ma(var, reaction=r) - ref_var)/ref_var) < 1e-14

    with pytest.raises(ValueError):
        Arrhenius([A, Ea_over_R, 1, A])

    # assert ama.as_mass_action({T_: 273.15}).args[0] == A*math.exp(-Ea_over_R/273.15)


@requires(units_library)
def test_specialfraction_rateexpr__units():
    _k, _kprime = 3.5 * u.s**-1 * u.molar**-0.5, 9.2
    rsys = _get_SpecialFraction_rsys(_k, _kprime)
    r = rsys.rxns[0]

    conc = {'H2': 2*u.molar, 'Br2': 3000*u.mol/u.m**3, 'HBr': 5*u.molar}
    ma = r.rate_expr()

    def _check(k, kprime, c):
        ref = k*c['H2']*c['Br2']**1.5/(c['Br2'] + kprime*c['HBr'])
        assert abs(ma(c, reaction=r) - ref) < 1e-15*u.molar/u.second

    _check(_k, _kprime, conc)
    alt_k, alt_kprime = 2*u.s**-1*u.molar**-0.5, 5
    _check(alt_k, alt_kprime, dict(list(conc.items()) + [
        ('k_HBr', alt_k), ('kprime_HBr', alt_kprime)]))


@requires(units_library, 'numpy')
@pytest.mark.parametrize('R_from_constants', [False, True])
def test_ArrheniusMassAction__units(R_from_constants):
    import numpy as np
    Ea = 40e3*u.J/u.mol
    R = default_constants.molar_gas_constant if R_from_constants else (8.3145*u.J/u.mol/u.K)
    A, Ea_over_R = 1.2e11/u.molar**2/u.second, Ea/R
    ref1 = A*np.exp(-to_unitless(Ea_over_R/(290*u.K)))
    arrh = Arrhenius([A, Ea_over_R])
    assert allclose(arrh({'temperature': 290*u.K}), ref1)
    ama = MassAction(arrh)
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, ama, {'B': 1})
    T_ = 'temperature'

    def ref(v):
        return 1.2e11/u.molar**2/u.second*math.exp(
            -Ea_over_R.simplified/v[T_])*v['B']*v['A']**2
    ma = r.rate_expr()

    for params in [(11.*u.molar, 13.*u.molar, 17.*u.molar, 311.2*u.kelvin),
                   (12*u.molar, 8*u.molar, 5*u.molar, 270*u.kelvin)]:
        var = dict(zip(['A', 'B', 'C', T_], params))
        ref_val = ref(var)
        assert abs((ma(var, reaction=r) - ref_val)/ref_val) < 1e-14


def test_Radiolytic():
    r = Radiolytic([2.1e-7])
    res = r({'doserate': 0.15, 'density': 0.998})
    assert abs(res - 0.15*0.998*2.1e-7) < 1e-15


@requires(parsing_library)
def test_Radiolytic__parsing():
    rxn = Reaction.from_string("-> H + OH; Radiolytic({'radiolytic_yield': 2.1e-7})", None)
    res = rxn.rate({'doserate': 0.15, 'density': 0.998})
    ref = 0.15*0.998*2.1e-7
    assert abs((res['H'] - ref)/ref) < 1e-15
    assert abs((res['OH'] - ref)/ref) < 1e-15
    gval, = rxn.rate_expr().g_values({}).values()
    assert abs(gval - 2.1e-7) < 1e-15


@requires(parsing_library, units_library)
def test_Radiolytic__parsing__units():
    rxn = Reaction.from_string("-> H + OH; Radiolytic({'radiolytic_yield': 2.1e-7*mol/J})", None)
    assert rxn.reac == {}
    assert rxn.prod == {'H': 1, 'OH': 1}
    res = rxn.rate({'doserate': 0.15*u.gray/u.s, 'density': 0.998*u.kg/u.dm3})
    ref = 0.15*0.998*2.1e-7*u.molar/u.second
    assert abs((res['H'] - ref)/ref) < 1e-15
    assert abs((res['OH'] - ref)/ref) < 1e-15


@requires(units_library)
def test_Radiolytic__units():

    def _check(r):
        res = r({'doserate': 0.15*u.gray/u.second,
                 'density': 0.998*u.kg/u.decimetre**3})
        ref = 0.15*0.998*2.1e-7*u.molar/u.second
        assert abs(to_unitless((res - ref)/ref)) < 1e-15

    _check(Radiolytic([2.1e-7*u.mol/u.joule]))
    _check(Radiolytic([2.0261921896167396*u.per100eV]))


@requires(units_library)
def test_Radioyltic__Reaction_html():
    rate = Radiolytic([2.1*u.per100eV])
    rxn = Reaction({}, {'H': 1}, rate)
    H = Substance.from_formula('H')
    html = rxn.html({'H': H}, with_param=True)
    assert html == ' &rarr; H&#59; %s' % str(rate)


def test_mk_Radiolytic():
    R1 = mk_Radiolytic()
    R2 = mk_Radiolytic()
    assert R1 is R2

    RABG = mk_Radiolytic('alpha', 'beta', 'gamma')
    rxn = Reaction({}, {'H': 2}, RABG([3, 5, 7], 'ya yb yg'.split()))
    rat = rxn.rate({'doserate_alpha': 11, 'doserate_beta': 13, 'doserate_gamma': 17, 'density': .7})
    assert abs(rat['H'] - .7*2*(3*11 + 5*13 + 7*17)) < 1e-13
    assert RABG.parameter_keys == ('density', 'doserate_alpha', 'doserate_beta', 'doserate_gamma')
    assert RABG.argument_names == tuple('radiolytic_yield_%s' % k for k in 'alpha beta gamma'.split())
    assert rxn.param.unique_keys == ('ya', 'yb', 'yg')
    rat2 = rxn.rate({'doserate_alpha': 11, 'doserate_beta': 13, 'doserate_gamma': 17, 'density': .7,
                     'ya': 23, 'yb': 29, 'yg': 31})
    assert abs(rat2['H'] - .7*2*(23*11 + 29*13 + 31*17)) < 1e-13


class Eyring4(Expr):
    nargs = 4

    def __call__(self, variables, reaction, backend=math):
        Sact_fact, Hact_exp, Sref_fact, Href_exp = self.all_args(variables, backend=backend)
        T = variables['temperature']
        return T * Sact_fact / Sref_fact * backend.exp((Href_exp-Hact_exp)/T)


def test_EyringMassAction():
    args = kB_h_times_exp_dS_R, dH_over_R, c0 = 1.2e11/273.15, 40e3/8, 1
    ama = MassAction(Eyring(args, ('Sfreq', 'Hact')))
    rxn1 = Reaction({'A': 2, 'B': 1}, {'C': 1}, ama, {'B': 1})
    T_ = 'temperature'

    def ref(v):
        return v.get('Sfreq', 1.2e11/273.15)*v[T_]*math.exp(-v.get('Hact', 40e3/8)/v[T_])*v['B']*v['A']**2

    for params in [(11., 13., 17., 311.2),
                   (12, 8, 5, 270)]:
        var = dict(zip(['A', 'B', 'C', T_], params))
        ref_val = ref(var)
        assert abs((ama(var, reaction=rxn1) - ref_val)/ref_val) < 1e-14

    with pytest.raises(ValueError):
        MassAction(Eyring([1, 1, 1, 1, 1]))

    # assert ama.as_mass_action({T_: 273.15}).args[0] == 1.2e11*math.exp(-40e3/8/273.15)

    ama2 = MassAction(Eyring4([1.2e11/273, 40e3/8, 1.2, 1e3], ('Sfreq', 'Hact', 'Sref', 'Href')))
    rxn2 = Reaction({'C': 1}, {'A': 2, 'B': 2}, ama2)
    var2 = {'C': 29, 'temperature': 273}

    def ref2(var):
        return var['C']*var.get('temperature', 273)*var.get('Sfreq', 1.2e11/273)/var.get('Sref', 1.2)*math.exp(
            (var.get('Href', 1e3) - var.get('Hact', 5e3))/var.get('temperature', 273))

    r2 = ref2(var2)
    assert abs((ama2(var2, reaction=rxn2) - r2)/r2) < 1e-14

    rsys = ReactionSystem([rxn1, rxn2])
    var3 = {'A': 11, 'B': 13, 'C': 17, 'temperature': 298, 'Sfreq': 1.2e11/298}
    rates = rsys.rates(var3)
    rf3 = ref(var3)
    rb3 = ref2(var3)
    ref_rates = {'A': 2*(rb3 - rf3), 'B': 2*(rb3 - rf3), 'C': rf3 - rb3}
    for k, v in ref_rates.items():
        assert abs((rates[k] - v)/v) < 1e-14
