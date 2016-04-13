# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


from ..henry import Henry, HenryWithUnits
from ..units import units_library, allclose, default_units as u
from ..util.testing import requires


def test_Henry():
    kH_O2 = Henry(1.2e-3, 1800, ref='carpenter_1966')
    assert abs(kH_O2(298.15) - 1.2e-3) < 1e-4
    assert abs(kH_O2.get_c_at_T_and_P(290, 1) - 0.001421892) < 1e-8
    assert abs(kH_O2.get_P_at_T_and_c(310, 1e-3) - 1.05) < 1e-3


@requires(units_library, 'numpy')
def test_Henry__with_units():
    import numpy as np

    kH_H2 = HenryWithUnits(7.8e-4*u.molar/u.atm, 640*u.K, ref='dean_1992')

    assert allclose(kH_H2.get_kH_at_T(
        300*u.K), 7.697430323e-4*u.molar/u.atm)
    kH = kH_H2.get_c_at_T_and_P(
        np.linspace(297.5, 298.65, 3)*u.K, .1*u.bar)
    assert allclose(kH, 7.7e-5*u.molar, rtol=1e-5, atol=1e-6*u.molar)
    kH = kH_H2.get_P_at_T_and_c(
        298.15*u.K, np.linspace(2e-3, 2.1e-3, 3)*u.molar)
    assert allclose(kH, 2.65*u.atm, rtol=1e-5, atol=0.2*u.bar)


@requires(units_library)
def test_HenryWithUnits():
    kH_H2 = HenryWithUnits(7.8e-4*u.molar/u.atm, 640*u.K, ref='dean_1992')
    Hcp = kH_H2(300*u.K)
    assert allclose(Hcp, 7.697430323e-4*u.molar/u.atm)
