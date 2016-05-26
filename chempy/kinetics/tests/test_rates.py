# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math

import pytest

from chempy import Reaction, ReactionSystem, Substance
from chempy.units import allclose, Backend, to_unitless, units_library, default_units as u
from chempy.util.parsing import parsing_library
from chempy.util.testing import requires
from ..rates import RateExpr, MassAction, ArrheniusMassAction, Radiolytic, mk_Radiolytic, EyringMassAction


class SpecialFraction(RateExpr):
    """
    Example from Atkins, De Paula, Physical Chemistry
    H2 + Br2 -> 2HBr
    ½ dHBr/dt = k[H2][Br2]**(3/2) / ([Br2] + k'[HBr])
    """
    def __call__(self, variables, backend=math):
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
    p = rsys.rxns[0].param

    conc = {'H2': 2, 'Br2': 3, 'HBr': 5}

    def _check(k, kprime, c):
        ref = k*c['H2']*c['Br2']**1.5/(c['Br2'] + kprime*c['HBr'])
        assert abs(p(c) - ref) < 1e-15

    _check(11, 13, conc)


def test_RateExpr__subclass_from_callback():
    SF = RateExpr.subclass_from_callback(
        lambda v, a, backend: a[0]*v['H2']*v['Br2']**(3/2) / (v['Br2'] + a[1]*v['HBr'])
    )
    ratex = SF([11, 13], ['k_HBr', 'kprime_HBr'])
    Reaction({'H2': 1, 'Br2': 1}, {'HBr': 2}, ratex)

    res1 = ratex({'H2': 5, 'Br2': 7, 'HBr': 15})
    ref1 = 11*5*7**1.5/(7+13*15)
    assert abs((res1-ref1)/ref1) < 1e-14

    res2 = ratex({'H2': 5, 'Br2': 7, 'HBr': 15, 'k_HBr': 23, 'kprime_HBr': 42})
    ref2 = 23*5*7**1.5/(7+42*15)
    assert abs((res2-ref2)/ref2) < 1e-14


def test_MassAction():
    ma = MassAction([3.14], ['my_rate'])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, ma, {'B': 1})
    assert abs(ma({'A': 11, 'B': 13, 'C': 17}) - 3.14*13*11**2) < 1e-14
    assert abs(ma({'A': 11, 'B': 13, 'C': 17, 'my_rate': 2.72}) - 2.72*13*11**2) < 1e-12
    ma2 = ma.as_mass_action({})
    assert ma == ma2
    assert abs(ma2({'A': 11, 'B': 13, 'C': 17}) - 3.14*13*11**2) < 1e-14


def test_MassAction__subclass_from_callback():
    def rate_coeff(variables, all_args, backend):
        return all_args[0]*backend.exp(all_args[1]/variables['temperature'])
    CustomMassAction = MassAction.subclass_from_callback(
        rate_coeff, cls_attrs=dict(parameter_keys=('temperature',), nargs=2))
    k1 = CustomMassAction([2.1e10, -5132.2], rxn=Reaction({'H2': 2, 'O2': 1}, {'H2O': 2}))
    res = k1({'temperature': 273.15, 'H2': 7, 'O2': 13})
    ref = 7*7*13*2.1e10*math.exp(-5132.2/273.15)
    assert abs((res-ref)/ref) < 1e-14


@requires(units_library)
def test_MassAction__subclass_from_callback__units():
    def rate_coeff(variables, all_args, backend):
        return all_args[0]*backend.exp(all_args[1]/variables['temperature'])
    CustomMassAction = MassAction.subclass_from_callback(
        rate_coeff, cls_attrs=dict(parameter_keys=('temperature',), nargs=2))
    k1 = CustomMassAction([2.1e10/u.molar**2/u.second, -5132.2*u.kelvin], rxn=Reaction({'H2': 2, 'O2': 1}, {'H2O': 2}))
    res = k1({'temperature': 491.67*u.rankine, 'H2': 7000*u.mol/u.metre**3, 'O2': 13*u.molar}, backend=Backend())
    ref = 7*7*13*2.1e10*math.exp(-5132.2/273.15) * u.molar/u.second
    assert allclose(res, ref)


def test_ArrheniusMassAction():
    A, Ea_over_R = 1.2e11, 40e3/8.3145
    ama = ArrheniusMassAction([A, Ea_over_R])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, ama, {'B': 1})
    T_ = 'temperature'

    def ref(v):
        return 1.2e11*math.exp(-Ea_over_R/v[T_])*v['B']*v['A']**2

    for params in [(11., 13., 17., 311.2),
                   (12, 8, 5, 270)]:
        var = dict(zip(['A', 'B', 'C', T_], params))
        r = ref(var)
        assert abs((ama(var) - r)/r) < 1e-14

    with pytest.raises(ValueError):
        ArrheniusMassAction([A, Ea_over_R, A])

    assert ama.as_mass_action({T_: 273.15}).args[0] == A*math.exp(-Ea_over_R/273.15)


@requires(units_library)
def test_specialfraction_rateexpr__units():
    _k, _kprime = 3.5 * u.s**-1 * u.molar**-0.5, 9.2
    rsys = _get_SpecialFraction_rsys(_k, _kprime)
    p = rsys.rxns[0].param

    conc = {'H2': 2*u.molar, 'Br2': 3000*u.mol/u.m**3, 'HBr': 5*u.molar}

    def _check(k, kprime, c):
        ref = k*c['H2']*c['Br2']**1.5/(c['Br2'] + kprime*c['HBr'])
        assert abs(p(c) - ref) < 1e-15*u.molar/u.second

    _check(_k, _kprime, conc)
    alt_k, alt_kprime = 2*u.s**-1*u.molar**-0.5, 5
    _check(alt_k, alt_kprime, dict(list(conc.items()) + [
        ('k_HBr', alt_k), ('kprime_HBr', alt_kprime)]))


@requires(units_library)
def test_ArrheniusMassAction__units():
    A, Ea_over_R = 1.2e11/u.molar**2/u.second, 40e3/8.3145*u.kelvin
    ama = ArrheniusMassAction([A, Ea_over_R])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, ama, {'B': 1})
    T_ = 'temperature'

    def ref(v):
        return 1.2e11/u.molar**2/u.second*math.exp(
            -Ea_over_R/v[T_])*v['B']*v['A']**2

    for params in [(11.*u.molar, 13.*u.molar, 17.*u.molar, 311.2*u.kelvin),
                   (12*u.molar, 8*u.molar, 5*u.molar, 270*u.kelvin)]:
        var = dict(zip(['A', 'B', 'C', T_], params))
        r = ref(var)
        assert abs((ama(var) - r)/r) < 1e-14


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


def test_EyringMassAction():
    args = kB_h_times_exp_dS_R, dH_over_R = 1.2e11/273.15, 40e3/8.3145
    ama = EyringMassAction(args)
    Reaction({'A': 2, 'B': 1}, {'C': 1}, ama, {'B': 1})
    T_ = 'temperature'

    def ref(v):
        return 1.2e11/273.15*v[T_]*math.exp(-40e3/8.3145/v[T_])*v['B']*v['A']**2

    for params in [(11., 13., 17., 311.2),
                   (12, 8, 5, 270)]:
        var = dict(zip(['A', 'B', 'C', T_], params))
        r = ref(var)
        assert abs((ama(var) - r)/r) < 1e-14

    with pytest.raises(ValueError):
        EyringMassAction([1, 1, 1])

    assert ama.as_mass_action({T_: 273.15}).args[0] == 1.2e11*math.exp(-40e3/8.3145/273.15)
