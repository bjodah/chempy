# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from .._solution import Solution, QuantityDict
from ..util.testing import requires
from ..units import units_library, to_unitless, default_units as u


@requires(units_library)
def test_QuantityDict():
    # QuantityDict * scalar_quantity
    c = QuantityDict(u.molar, {})
    c['H2O'] = 55.4 * u.molar
    V = .4*u.dm3
    n = c*V
    assert isinstance(n, QuantityDict)

    # For the following to work major changes to quantities needs to be made:
    # n = V*c
    # assert isinstance(n, QuantityDict)

    assert n.isclose({'H2O': 55.4*.4*u.mol})
    assert abs(to_unitless(n['H2O'], u.mol) - 55.4*.4) < 1e-14


@requires(units_library)
def test_Solution():
    s1 = Solution(0.1*u.dm3, {'CH3OH': 0.1 * u.molar})
    s2 = Solution(0.3*u.dm3, {'CH3OH': 0.4 * u.molar, 'Na+': 2e-3*u.molar,
                              'Cl-': 2e-3*u.molar})
    s3 = s1 + s2
    assert abs(to_unitless(s3.volume - 4e-4 * u.m**3, u.dm3)) < 1e-15

    assert s3.concentrations.isclose({'CH3OH': 0.325*u.molar, 'Na+': 1.5e-3*u.molar,
                                      'Cl-': 1.5e-3*u.molar})

    s4 = s3.dissolve({'CH3OH': 1*u.gram})
    mw = (12.011 + 4*1.008 + 15.999)
    assert abs(s4.concentrations['CH3OH'] - (0.325 + 1/mw/.4)*u.molar) < 1e-7
