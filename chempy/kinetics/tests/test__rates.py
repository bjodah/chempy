# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math

import pytest

from chempy import Reaction
from chempy.units import allclose, Backend, to_unitless, units_library, default_units as u
from chempy.util.testing import requires
from .._rates import (
    TPolyMassAction, RTPolyMassAction, Log10TPolyMassAction, TPolyInLog10MassAction, TPolyRadiolytic,
    TPoly, RTPoly, PiecewiseTPolyMassAction, Log10PiecewiseRTPolyMassAction, TPiecewise
)


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


@requires(units_library)
def test_TPiecewise():
    expr0 = TPolyMassAction([273.15*u.K, 10/u.molar/u.s, 0.1/u.molar/u.s/u.K])
    expr1 = TPolyMassAction([298.15*u.K, 12.5/u.molar/u.s, 0/u.molar/u.s/u.K, 2/u.molar/u.s/u.K**2])
    pw = TPiecewise([273.15*u.K, 298.15*u.K, expr0, 298.15*u.K, 373.15*u.K, expr1])
    Reaction({'e-(aq)': 2}, {'H2': 1, 'OH-': 2}, pw, {'H2O': 2})
    res0 = pw({'temperature': 293.15*u.K, 'e-(aq)': 1e-13*u.molar})
    ref0 = 12*1e-26 * u.molar/u.s
    assert allclose(res0, ref0)
    assert not allclose(res0, 2*ref0)

    res1 = pw({'temperature': 300.15*u.K, 'e-(aq)': 2e-13*u.molar})
    ref1 = 20.5*4e-26 * u.molar/u.s
    assert allclose(res1, ref1)
    assert not allclose(res1, ref1/2)
