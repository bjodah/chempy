# -*- coding: utf-8 -*-


import math

import pytest

from ..eyring import eyring_equation, EyringParam, EyringParamWithUnits
from chempy.util.testing import requires
from chempy.units import allclose, units_library, default_units as u, default_constants

_kB_over_h = 2.083664399411865234375e10
_R = 8.314472


def test_eyring_equation():
    dH, dS = 40e3, 1e2
    T = 123.45
    ref = _kB_over_h * T * math.exp(dS / _R) * math.exp(-dH / _R / T)
    res = eyring_equation(dH, dS, T)
    assert abs((res - ref) / ref) < 1e-10


def test_EyringParam():
    T = 273.15
    k = EyringParam(40e3, 1e2)(T)
    ref = _kB_over_h * T * math.exp(1e2 / _R) * math.exp(-40e3 / _R / T)
    assert abs((k - ref) / ref) < 1e-7


@requires(units_library)
def test_EyringParamWithUnits():
    mol = u.mol
    J = u.joule
    K = u.kelvin
    T = 273.15 * K
    k = EyringParamWithUnits(42e3 * J / mol, 123 * J / mol / K)(T)
    ref = _kB_over_h * 273.15 * math.exp(123 / _R) * math.exp(-42e3 / _R / 273.15)
    ref /= u.second
    assert allclose(k, ref)


@requires(units_library)
def test_eyring_equation__units_constants_None():
    mol = u.mol
    J = u.joule
    K = u.kelvin
    T = 273.15 * K
    dH = 42e3 * J / mol
    dS = 123 * J / mol / K
    k = eyring_equation(dH, dS, T, constants=None, units=u)
    ref = _kB_over_h * 273.15 * math.exp(123 / _R) * math.exp(-42e3 / _R / 273.15)
    ref /= u.second
    assert allclose(k, ref)

    ep = EyringParam(dH, dS)
    assert allclose(ep(T, units=u), ref)
    assert allclose(ep.dH_over_R(units=u), dH / (_R * J / mol / K))


@requires(units_library)
@pytest.mark.parametrize("constants", [None, default_constants])
@pytest.mark.parametrize("backend", ["math", "numpy"])
def test_eyring_equation__scaled_units(constants, backend):
    backend = __import__(backend)
    dH = 42 * u.kilojoule / u.mol
    dS = 0.123 * u.kilojoule / u.mol / u.K
    T = 273150 * u.mK
    kwargs = {"units": u, "constants": constants}
    reference = eyring_equation(
        42e3 * u.J / u.mol, 123 * u.J / u.mol / u.K, 273.15 * u.K, **kwargs
    )
    assert allclose(eyring_equation(dH, dS, T, backend=backend, **kwargs), reference)

    ep = EyringParam(dH, dS)
    prefactor = ep.kB_h_times_exp_dS_R(backend=backend, **kwargs)
    expected = EyringParam(42e3 * u.J / u.mol, 123 * u.J / u.mol / u.K)
    assert allclose(prefactor, expected.kB_h_times_exp_dS_R(**kwargs))
