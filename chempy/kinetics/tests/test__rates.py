# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math

import pytest

from chempy import Reaction
from chempy.units import allclose, Backend, to_unitless, units_library, default_units as u
from chempy.util._expr import Log10, Constant
from chempy.util.testing import requires
from ..rates import MassAction, Radiolytic
from .._rates import TPoly, RTPoly, ShiftedTPoly, ShiftedRTPoly, ShiftedLog10TPoly, TPiecewise


def test_TPolyMassAction():
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, inact_reac={'B': 1})
    p = MassAction(ShiftedTPoly([273.15, 7, .2, .03, .004]))
    res = p({'A': 11, 'B': 13, 'temperature': 298.15}, reaction=r)
    ref = 7 + .2*25 + 0.03 * 25**2 + 0.004 * 25**3
    assert abs(res - ref*13*11**2) < 1e-15


def test_ShiftedTPoly_MassAction():
    rate = MassAction(ShiftedTPoly([100, 2, 5, 7]))
    assert rate.args[0].args == [100, 2, 5, 7]
    r = Reaction({'A': 2, 'B': 1}, {'P': 1}, rate)
    res = r.rate_expr()({'A': 11, 'B': 13, 'temperature': 273.15}, reaction=r)
    x = 273.15-100
    ref = 11*11*13*(2 + 5*x + 7*x**2)
    assert abs((res - ref)/ref) < 1e-14


@requires(units_library)
def test_ShiftedTPoly__units():
    stp1 = ShiftedTPoly([273.15*u.kelvin, 5, 7/u.kelvin])
    allclose(stp1({'temperature': 274.15*u.kelvin}), 5+7)
    allclose(stp1({'temperature': 273.15*u.kelvin}), 5)
    stp2 = ShiftedTPoly([273.15*u.kelvin, 5*u.m, 7*u.m/u.kelvin, 13*u.m*u.kelvin**-2])
    allclose(stp2({'temperature': 274.15*u.kelvin}), (5+7+13)*u.m)
    allclose(stp2({'temperature': 273.15*u.kelvin}), 5*u.m)


@requires(units_library)
def test_TPolyMassAction__units():
    Mps = u.molar/u.second
    kunit = 1/u.molar**2/u.second
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, None, {'B': 1})
    p = MassAction(ShiftedTPoly([273.15*u.K, 7*kunit, .2*kunit/u.K, .03*kunit/u.K**2, .004*kunit/u.K**3]))
    res = p({'A': 11*u.molar, 'B': 13*u.molar, 'temperature': 298.15*u.K}, reaction=r)
    ref = 7 + .2*25 + 0.03 * 25**2 + 0.004 * 25**3
    assert abs(res - ref*13*11**2*Mps) < 1e-15


def test_RTPolyMassAction():
    p = MassAction(ShiftedRTPoly([273.15, 7, .2, .03, .004]))
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, inact_reac={'B': 1})
    res = p({'A': 11, 'B': 13, 'temperature': 298.15}, reaction=r)
    ref = 7 + .2/25 + 0.03 / 25**2 + 0.004 / 25**3
    assert abs(res - ref*13*11**2) < 1e-15


@requires(units_library)
def test_RTPolyMassAction__units():
    Mps = u.molar/u.second
    kunit = 1/u.molar**2/u.second
    p = MassAction(ShiftedRTPoly([273.15*u.K, 7*kunit, .2*kunit*u.K, .03*kunit*u.K**2, .004*kunit*u.K**3]))
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, None, {'B': 1})
    res = p({'A': 11*u.molar, 'B': 13*u.molar, 'temperature': 298.15*u.K}, reaction=r)
    ref = 7 + .2/25 + 0.03 / 25**2 + 0.004 / 25**3
    assert abs(res - ref*13*11**2*Mps) < 1e-15


def test_Log10TPolyMassAction():
    p = MassAction(10**ShiftedTPoly([273.15, .7, .02, .003, .0004]))
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, p, {'B': 1})
    res = p({'A': 11, 'B': 13, 'temperature': 298.15}, reaction=r)
    ref = 10**(.7 + .02*25 + 0.003 * 25**2 + 0.0004 * 25**3)
    assert abs(res - ref*13*11**2) < 1e-15


@requires(units_library)
def test_Log10TPolyMassAction__units():
    Mps = u.molar/u.second
    kunit = 1/u.molar**2/u.second
    p = MassAction(Constant(kunit)*10**ShiftedTPoly([273.15*u.K, .7, .02/u.K, .003/u.K**2, .0004/u.K**3]))
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, p, {'B': 1})
    res = p({'A': 11*u.molar, 'B': 13*u.molar, 'temperature': 298.15*u.K}, reaction=r)
    ref = 10**(.7 + .02*25 + 0.003 * 25**2 + 0.0004 * 25**3)
    assert abs(res - ref*13*11**2*Mps) < 1e-15


def test_TPolyInLog10MassAction():
    p = MassAction(ShiftedLog10TPoly([2, 0.3, .2, .03, .004]))
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, p, {'B': 1})
    lgT = Log10('temperature')
    lgTref = Log10('Tref')
    res = p({'A': 11, 'B': 13, 'temperature': 298.15,
             'log10_temperature': lgT, 'log10_Tref': lgTref}, reaction=r)
    _T = math.log10(298.15) - 2
    ref = .3 + .2*_T + 0.03 * _T**2 + 0.004 * _T**3
    assert abs(res - ref*13*11**2) < 1e-15


@requires(units_library)
def test_TPolyInLog10MassAction__units():
    Mps = u.molar/u.second
    kunit = 1/u.molar**2/u.second
    p = MassAction(Constant(kunit)*ShiftedLog10TPoly([2, 0.3, .2, .03, .004]))
    lgT = Log10('temperature'/Constant(u.K))
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, p, {'B': 1})
    res = p({'A': 11*u.molar, 'B': 13*u.molar, 'temperature': 298.15*u.K,
             'log10_temperature': lgT}, backend=Backend(), reaction=r)
    _T = math.log10(298.15) - 2
    ref = (.3 + .2*_T + 0.03 * _T**2 + 0.004 * _T**3) * 13 * 11**2 * Mps
    assert abs(res - ref) < 1e-15*Mps


def test_TPolyRadiolytic():
    tpr = Radiolytic(ShiftedTPoly([273.15, 1.85e-7, 1e-9]))
    res = tpr({'doserate': 0.15, 'density': 0.998, 'temperature': 298.15})
    assert abs(res - 0.15*0.998*2.1e-7) < 1e-15


@requires(units_library)
def test_TPolyRadiolytic__units():

    def _check(r):
        res = r({'doserate': 0.15*u.gray/u.second,
                 'density': 0.998*u.kg/u.decimetre**3,
                 'temperature': 298.15*u.K})
        ref = 0.15*0.998*2.1e-7*u.molar/u.second
        assert abs(to_unitless((res - ref)/ref)) < 1e-15

    _check(Radiolytic(ShiftedTPoly([273.15*u.K, 1.85e-7*u.mol/u.joule, 1e-9*u.mol/u.joule/u.K])))


def test_TPiecewisePolyMassAction():
    tp1 = TPoly([10, 0.1])
    tp2 = ShiftedTPoly([273.15, 37.315, -0.1])
    pwp = MassAction(TPiecewise([0, tp1, 273.15, tp2, 373.15]))
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, inact_reac={'B': 1})
    res1 = pwp({'A': 11, 'B': 13, 'temperature': 198.15}, reaction=r)
    ref1 = 11*11*13*29.815
    assert abs((res1-ref1)/ref1) < 1e-14
    res2 = pwp({'A': 11, 'B': 13, 'temperature': 298.15}, reaction=r)
    ref2 = 11*11*13*(37.315 - 25*0.1)
    assert abs((res2-ref2)/ref2) < 1e-14
    with pytest.raises(ValueError):
        pwp({'A': 11, 'B': 13, 'temperature': 398.15}, reaction=r)


@requires(units_library)
def test_PiecewiseTPolyMassAction__units():
    tp1 = TPoly([10/u.molar**2/u.second, 0.1/u.molar**2/u.second/u.kelvin])
    tp2 = ShiftedTPoly([273.15*u.kelvin, 37.315/u.molar**2/u.second, -0.1/u.molar**2/u.second/u.kelvin])
    pwp = MassAction(TPiecewise([0*u.kelvin, tp1, 273.15*u.kelvin, tp2, float('inf')*u.kelvin]))
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, pwp, {'B': 1})
    res1 = pwp({'A': 11*u.molar, 'B': 13*u.molar, 'temperature': 198.15*u.kelvin}, reaction=r)
    ref1 = 11*11*13*29.815 * u.molar/u.second
    assert abs((res1-ref1)/ref1) < 1e-14


@requires('sympy')
def test_PiecewiseTPolyMassAction__sympy():
    import sympy as sp
    tp1 = TPoly([10, 0.1])
    tp2 = ShiftedTPoly([273.15, 37.315, -0.1])
    pwp = MassAction(TPiecewise([0, tp1, 273.15, tp2, 373.15]))
    T = sp.Symbol('T')
    r = Reaction({'A': 2, 'B': 1}, {'C': 1}, inact_reac={'B': 1})
    res1 = pwp({'A': 11, 'B': 13, 'temperature': T}, backend=sp, reaction=r)
    ref1 = 11**2 * 13 * sp.Piecewise(
        (10+0.1*T, sp.And(0 <= T, T <= 273.15)),
        (37.315 - 0.1*(T-273.15), sp.And(273.15 <= T, T <= 373.15)),
        (sp.Symbol('NAN'), True)
    )
    assert res1 == ref1


def test_Log10PiecewiseRTPolyMassAction():
    p1 = RTPoly([12.281, -3.768e2, -6.673e4, -1.075e7])
    p2 = RTPoly([-47.532, 4.92, -1.036, 0.0])
    ratex = MassAction(10**TPiecewise([293.15, p1, 423.15, p2, 623.15]))
    r = Reaction({'e-(aq)': 2}, {'H2': 1, 'OH-': 2}, ratex, {'H2O': 2})
    res = ratex({'e-(aq)': 1e-13, 'temperature': 293.15}, reaction=r)
    ref = 6.20e9*1e-26
    assert abs((res-ref)/ref) < 6e-3


@requires(units_library)
def test_TPiecewise():
    expr0 = ShiftedTPoly([273.15*u.K, 10/u.molar/u.s, 0.1/u.molar/u.s/u.K])
    expr1 = ShiftedTPoly([298.15*u.K, 12.5/u.molar/u.s, 0/u.molar/u.s/u.K, 2/u.molar/u.s/u.K**2])
    pwma = MassAction(TPiecewise([273.15*u.K, expr0, 298.15*u.K, expr1, 373.15*u.K]))
    r = Reaction({'e-(aq)': 2}, {'H2': 1, 'OH-': 2}, inact_reac={'H2O': 2})
    res0 = pwma({'temperature': 293.15*u.K, 'e-(aq)': 1e-13*u.molar}, reaction=r)
    ref0 = 12*1e-26 * u.molar/u.s
    assert allclose(res0, ref0)
    assert not allclose(res0, 2*ref0)

    res1 = pwma({'temperature': 300.15*u.K, 'e-(aq)': 2e-13*u.molar}, reaction=r)
    ref1 = 20.5*4e-26 * u.molar/u.s
    assert allclose(res1, ref1)
    assert not allclose(res1, ref1/2)


@requires(units_library)
def test_MassAction__rate_coeff():
    perMolar_perSecond = u.perMolar_perSecond

    p1 = MassAction(Constant(perMolar_perSecond) * 10**RTPoly([1, 2*u.kelvin, 3*u.kelvin**2]))
    rcoeff1 = p1.rate_coeff({'temperature': 283.15*u.K})
    ref1 = 10**(1 + 2/283.15 + 3/283.15**2) * perMolar_perSecond
    assert allclose(rcoeff1, ref1)
    rxn1 = Reaction({'A', 'B'}, {'C'}, p1)
    rat1 = rxn1.rate({'A': 2, 'B': 3, 'temperature': 283.15*u.K})
    assert allclose(rat1['A'], -2*3*ref1)
    assert allclose(rat1['B'], -2*3*ref1)
    assert allclose(rat1['C'], 2*3*ref1)
