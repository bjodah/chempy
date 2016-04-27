# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math

import pytest

from chempy import Reaction, ReactionSystem, Substance
from chempy.units import Backend, to_unitless, units_library, default_units as u
from chempy.util.testing import requires
from ..rates import (
    RateExpr, MassAction, ArrheniusMassAction, Radiolytic, TPolyMassAction,
    RTPolyMassAction, Log10TPolyMassAction, TPolyInLog10MassAction, TPolyRadiolytic,
    TPoly, RTPoly, PiecewiseTPolyMassAction, Log10PiecewiseRTPolyMassAction
)


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


def test_RateExpr__from_callback():
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


@requires(units_library)
def test_Radiolytic__units():

    def _check(r):
        res = r({'doserate': 0.15*u.gray/u.second,
                 'density': 0.998*u.kg/u.decimetre**3})
        ref = 0.15*0.998*2.1e-7*u.molar/u.second
        assert abs(to_unitless((res - ref)/ref)) < 1e-15

    _check(Radiolytic([2.1e-7*u.mol/u.joule]))
    _check(Radiolytic([2.0261921896167396*u.per100eV]))


def test_TPolyMassAction():
    r = TPolyMassAction([273.15, 7, .2, .03, .004])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, r, {'B': 1})
    res = r({'A': 11, 'B': 13, 'temperature': 298.15})
    ref = 7 + .2*25 + 0.03 * 25**2 + 0.004 * 25**3
    assert abs(res - ref*13*11**2) < 1e-15


def test_TPolyMassAction__2():
    rate = TPolyMassAction([100, 2, 5, 7])
    assert rate.args == [100, 2, 5, 7]
    Reaction({'A': 2, 'B': 1}, {'P': 1}, rate)
    res = rate({'A': 11, 'B': 13, 'temperature': 273.15})
    x = 273.15-100
    ref = 11*11*13*(2 + 5*x + 7*x**2)
    assert abs((res - ref)/ref) < 1e-14


@requires(units_library)
def test_TPolyMassAction__units():
    Mps = u.molar/u.second
    kunit = 1/u.molar**2/u.second
    r = TPolyMassAction([273.15*u.K, 7*kunit, .2*kunit/u.K, .03*kunit/u.K**2, .004*kunit/u.K**3])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, r, {'B': 1})
    res = r({'A': 11*u.molar, 'B': 13*u.molar, 'temperature': 298.15*u.K})
    ref = 7 + .2*25 + 0.03 * 25**2 + 0.004 * 25**3
    assert abs(res - ref*13*11**2*Mps) < 1e-15


def test_RTPolyMassAction():
    r = RTPolyMassAction([273.15, 7, .2, .03, .004])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, r, {'B': 1})
    res = r({'A': 11, 'B': 13, 'temperature': 298.15})
    ref = 7 + .2/25 + 0.03 / 25**2 + 0.004 / 25**3
    assert abs(res - ref*13*11**2) < 1e-15


@requires(units_library)
def test_RTPolyMassAction__units():
    Mps = u.molar/u.second
    kunit = 1/u.molar**2/u.second
    r = RTPolyMassAction([273.15*u.K, 7*kunit, .2*kunit*u.K, .03*kunit*u.K**2, .004*kunit*u.K**3])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, r, {'B': 1})
    res = r({'A': 11*u.molar, 'B': 13*u.molar, 'temperature': 298.15*u.K})
    ref = 7 + .2/25 + 0.03 / 25**2 + 0.004 / 25**3
    assert abs(res - ref*13*11**2*Mps) < 1e-15


def test_Log10TPolyMassAction():
    r = Log10TPolyMassAction([1, 273.15, .7, .02, .003, .0004])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, r, {'B': 1})
    res = r({'A': 11, 'B': 13, 'temperature': 298.15})
    ref = 10**(.7 + .02*25 + 0.003 * 25**2 + 0.0004 * 25**3)
    assert abs(res - ref*13*11**2) < 1e-15


@requires(units_library)
def test_Log10TPolyMassAction__units():
    Mps = u.molar/u.second
    kunit = 1/u.molar**2/u.second
    r = Log10TPolyMassAction([kunit, 273.15*u.K, .7, .02/u.K, .003/u.K**2, .0004/u.K**3])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, r, {'B': 1})
    res = r({'A': 11*u.molar, 'B': 13*u.molar, 'temperature': 298.15*u.K})
    ref = 10**(.7 + .02*25 + 0.003 * 25**2 + 0.0004 * 25**3)
    assert abs(res - ref*13*11**2*Mps) < 1e-15


def test_TPolyInLog10MassAction():
    r = TPolyInLog10MassAction([1, 2, 0.3, .2, .03, .004])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, r, {'B': 1})
    res = r({'A': 11, 'B': 13, 'temperature': 298.15})
    _T = math.log10(298.15) - 2
    ref = .3 + .2*_T + 0.03 * _T**2 + 0.004 * _T**3
    assert abs(res - ref*13*11**2) < 1e-15


@requires(units_library)
def test_TPolyInLog10MassAction__units():
    Mps = u.molar/u.second
    kunit = 1/u.molar**2/u.second
    r = TPolyInLog10MassAction([u.K, 2, 0.3*kunit, .2*kunit, .03*kunit, .004*kunit])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, r, {'B': 1})
    res = r({'A': 11*u.molar, 'B': 13*u.molar, 'temperature': 298.15*u.K}, backend=Backend())
    _T = math.log10(298.15) - 2
    ref = .3 + .2*_T + 0.03 * _T**2 + 0.004 * _T**3
    assert abs(res - ref*13*11**2*Mps) < 1e-15


def test_TPolyRadiolytic():
    r = TPolyRadiolytic([273.15, 1.85e-7, 1e-9])
    res = r({'doserate': 0.15, 'density': 0.998, 'temperature': 298.15})
    assert abs(res - 0.15*0.998*2.1e-7) < 1e-15


@requires(units_library)
def test_TPolyRadiolytic__units():

    def _check(r):
        res = r({'doserate': 0.15*u.gray/u.second,
                 'density': 0.998*u.kg/u.decimetre**3,
                 'temperature': 298.15*u.K})
        ref = 0.15*0.998*2.1e-7*u.molar/u.second
        assert abs(to_unitless((res - ref)/ref)) < 1e-15

    _check(TPolyRadiolytic([273.15*u.K, 1.85e-7*u.mol/u.joule, 1e-9*u.mol/u.joule/u.K]))


@requires(units_library)
def test_Radioyltic__Reaction_html():
    rate = Radiolytic([2.1*u.per100eV])
    rxn = Reaction({}, {'H': 1}, rate)
    H = Substance.from_formula('H')
    html = rxn.html({'H': H}, with_param=True)
    assert html == ' &rarr; H&#59; %s' % str(rate)


def test_PiecewiseTPolyMassAction():
    tp1 = TPoly([0, 10, 0.1])
    tp2 = TPoly([273.15, 37.315, -0.1])
    pwp = PiecewiseTPolyMassAction.from_polynomials([(0, 273.15), (273.15, 373.15)], [tp1, tp2])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, pwp, {'B': 1})
    res1 = pwp({'A': 11, 'B': 13, 'temperature': 198.15})
    ref1 = 11*11*13*29.815
    assert abs((res1-ref1)/ref1) < 1e-14
    res2 = pwp({'A': 11, 'B': 13, 'temperature': 298.15})
    ref2 = 11*11*13*(37.315 - 25*0.1)
    assert abs((res2-ref2)/ref2) < 1e-14
    with pytest.raises(ValueError):
        pwp({'A': 11, 'B': 13, 'temperature': 398.15})


@requires(units_library)
def test_PiecewiseTPolyMassAction__units():
    tp1 = TPoly([0*u.kelvin, 10/u.molar**2/u.second, 0.1/u.molar**2/u.second/u.kelvin])
    tp2 = TPoly([273.15*u.kelvin, 37.315/u.molar**2/u.second, -0.1/u.molar**2/u.second/u.kelvin])
    pwp = PiecewiseTPolyMassAction.from_polynomials([
        (0*u.kelvin, 273.15*u.kelvin),
        (273.15*u.kelvin, 373.15*u.kelvin)
    ], [tp1, tp2])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, pwp, {'B': 1})
    res1 = pwp({'A': 11*u.molar, 'B': 13*u.molar, 'temperature': 198.15*u.kelvin})
    ref1 = 11*11*13*29.815 * u.molar/u.second
    assert abs((res1-ref1)/ref1) < 1e-14


@requires('sympy')
def test_PiecewiseTPolyMassAction__sympy():
    import sympy as sp
    tp1 = TPoly([0, 10, 0.1])
    tp2 = TPoly([273.15, 37.315, -0.1])
    pwp = PiecewiseTPolyMassAction.from_polynomials([(0, 273.15), (273.15, 373.15)], [tp1, tp2])
    T = sp.Symbol('T')
    Reaction({'A': 2, 'B': 1}, {'C': 1}, pwp, {'B': 1})
    res1 = pwp({'A': 11, 'B': 13, 'temperature': T}, backend=sp)
    ref1 = 11**2 * 13 * sp.Piecewise(
        (10+0.1*T, sp.And(0 <= T, T <= 273.15)),
        (37.315 - 0.1*(T-273.15), sp.And(273.15 <= T, T <= 373.15))
    )
    assert res1 == ref1


def test_Log10PiecewiseRTPolyMassAction():
    p1 = RTPoly([0, 12.281, -3.768e2, -6.673e4, -1.075e7])
    p2 = RTPoly([0, -47.532, 4.92, -1.036, 0.0])
    bounds = [(293.15, 423.15), (423.15, 623.15)]
    k_unit = 1
    ratex = Log10PiecewiseRTPolyMassAction.from_polynomials(bounds, [p1, p2], [k_unit])
    Reaction({'e-(aq)': 2}, {'H2': 1, 'OH-': 2}, ratex, {'H2O': 2})
    res = ratex({'e-(aq)': 1e-13, 'temperature': 293.15})
    ref = 6.20e9*1e-26
    assert abs((res-ref)/ref) < 6e-3
