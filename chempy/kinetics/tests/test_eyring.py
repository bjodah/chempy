# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function)

import math

from ..eyring import (
    eyring_equation, EyringParam, EyringParamWithUnits
)
from chempy.util.testing import requires
from chempy.units import allclose, units_library, default_units as u


_kB_over_h = 2.083664399411865234375e10
_R = 8.314472


def test_eyring_equation():
    dH, dS = 40e3, 1e2
    T = 123.45
    ref = _kB_over_h * T * math.exp(dS/_R) * math.exp(-dH/_R/T)
    res = eyring_equation(dH, dS, T)
    assert abs((res - ref)/ref) < 1e-10


def test_EyringParam():
    T = 273.15
    k = EyringParam(40e3, 1e2)(T)
    ref = _kB_over_h * T * math.exp(1e2/_R) * math.exp(-40e3/_R/T)
    assert abs((k - ref)/ref) < 1e-7


@requires(units_library)
def test_EyringParamWithUnits():
    mol = u.mol
    J = u.joule
    K = u.kelvin
    T = 273.15*K
    k = EyringParamWithUnits(42e3 * J/mol, 123*J/mol/K)(T)
    ref = _kB_over_h * 273.15 * math.exp(123/_R) * math.exp(-42e3/_R/273.15)
    ref /= u.second
    assert allclose(k, ref)
