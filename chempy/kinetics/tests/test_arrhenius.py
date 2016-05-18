# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function)

import math

from ..arrhenius import (
    arrhenius_equation, ArrheniusParam, ArrheniusParamWithUnits
)
from chempy.util.testing import requires
from chempy.units import default_units, units_library


def test_arrhenius_equation():
    assert abs(arrhenius_equation(3, 831.4472, 100) - 3/2.7182818) < 1e-7


def test_ArrheniusParam():
    k = ArrheniusParam(1e10, 42e3)(273.15)
    ref = 1e10 * math.exp(-42e3/(8.3145*273.15))
    assert abs((k - ref)/ref) < 1e-4


@requires(units_library)
def test_ArrheniusParamWithUnits():
    s = default_units.second
    mol = default_units.mol
    J = default_units.joule
    K = default_units.kelvin
    k = ArrheniusParamWithUnits(1e10/s, 42e3 * J/mol)(273.15*K)
    ref = 1e10/s * math.exp(-42e3/(8.3145*273.15))
    assert abs((k - ref)/ref) < 1e-4
