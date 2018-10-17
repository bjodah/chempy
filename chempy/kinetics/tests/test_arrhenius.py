# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function)

import math

from chempy.chemistry import Reaction
from chempy.util.testing import requires
from chempy.units import units_library, allclose, default_units as u

from ..arrhenius import (
    arrhenius_equation, ArrheniusParam, ArrheniusParamWithUnits
)


def test_arrhenius_equation():
    assert abs(arrhenius_equation(3, 831.4472, 100) - 3/2.7182818) < 1e-7


def test_ArrheniusParam():
    k = ArrheniusParam(1e10, 42e3)(273.15)
    ref = 1e10 * math.exp(-42e3/(8.3145*273.15))
    assert abs((k - ref)/ref) < 1e-4


@requires(units_library)
def test_ArrheniusParamWithUnits():
    s = u.second
    mol = u.mol
    J = u.joule
    K = u.kelvin
    act_J__mol = 42e3
    ap = ArrheniusParamWithUnits(1e10/s, act_J__mol * J/mol)
    freezing_K = 273.15
    k = ap(freezing_K*K)
    ref = 1e10/s * math.exp(-act_J__mol/(8.3145*freezing_K))
    assert abs((k - ref)/ref) < 1e-4

    r = Reaction({'H2O2': 1}, {'OH': 2}, ap)
    ratc = r.rate_expr().rate_coeff({'temperature': freezing_K*u.K})
    assert allclose(ratc, ref, rtol=1e-4)
