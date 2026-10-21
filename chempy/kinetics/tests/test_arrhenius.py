# -*- coding: utf-8 -*-


import math

import pytest

from chempy.chemistry import Reaction
from chempy.util.testing import requires
from chempy.units import units_library, allclose, default_units as u

from ..arrhenius import arrhenius_equation, ArrheniusParam, ArrheniusParamWithUnits


def test_arrhenius_equation():
    assert abs(arrhenius_equation(3, 831.4472, 100) - 3 / 2.7182818) < 1e-7


_A1, _Ea1, _T1, _k1 = 1e10, 42e3, 273.15, 1e10 * math.exp(-42e3 / (8.3145 * 273.15))


def test_ArrheniusParam():
    k = ArrheniusParam(_A1, _Ea1)(_T1)
    assert abs((k - _k1) / _k1) < 1e-4


@requires("numpy")
def test_ArrheniusParam__from_rateconst_at_T():
    ap = ArrheniusParam.from_rateconst_at_T(_Ea1, (_T1, _k1))
    assert abs((ap.A - _A1) / _A1) < 1e-4


def _get_ref2_units():
    A__s = 1e10
    act_J__mol = 42e3
    freezing_K = 273.15

    class ValueHolder:
        A = A__s / u.s
        Ea = act_J__mol * u.J / u.mol
        T = freezing_K * u.K
        k = A__s / u.s * math.exp(-act_J__mol / (8.3145 * freezing_K))

    return ValueHolder()


@requires(units_library)
def test_ArrheniusParamWithUnits():
    _2 = _get_ref2_units()
    ap = ArrheniusParamWithUnits(_2.A, _2.Ea)

    k = ap(_2.T)
    assert abs((k - _2.k) / _2.k) < 1e-4

    r = Reaction({"H2O2": 1}, {"OH": 2}, ap)
    ratc = r.rate_expr().rate_coeff({"temperature": _2.T})
    assert allclose(ratc, _2.k, rtol=1e-4)


@requires(units_library)
def test_ArrheniusParamWithUnits__from_rateconst_at_T():
    _2 = _get_ref2_units()
    apu = ArrheniusParamWithUnits.from_rateconst_at_T(_2.Ea, (_2.T, _2.k))
    assert allclose(apu(_2.T), _2.k)


@requires(units_library)
def test_arrhenius_equation__units_constants_None():
    _2 = _get_ref2_units()
    k = arrhenius_equation(_2.A, _2.Ea, _2.T, constants=None, units=u)
    assert abs((k - _2.k) / _2.k) < 1e-4

    ap = ArrheniusParam.from_rateconst_at_T(_2.Ea, (_2.T, _2.k), units=u)
    assert abs((ap.A - _2.A) / _2.A) < 1e-4


@requires(units_library)
def test_arrhenius_equation__legacy_unit_namespace():
    from types import SimpleNamespace

    legacy_units = SimpleNamespace(Joule=u.joule, Kelvin=u.kelvin, mol=u.mol)
    ref = _get_ref2_units()
    result = arrhenius_equation(ref.A, ref.Ea, ref.T, units=legacy_units)
    assert allclose(result, ref.k, rtol=1e-4)


@requires(units_library)
@pytest.mark.parametrize("backend", ["math", "numpy"])
def test_ArrheniusParam__from_rateconst_at_T_scaled_units(backend):
    ref = _get_ref2_units()
    Ea = ref.Ea.rescale(u.kilojoule / u.mol)
    T = ref.T.rescale(u.mK)
    ap = ArrheniusParam.from_rateconst_at_T(
        Ea, (T, ref.k), units=u, backend=__import__(backend)
    )
    assert allclose(ap.A, ref.A, rtol=1e-4)
    assert allclose(ap(T, units=u), ref.k)
