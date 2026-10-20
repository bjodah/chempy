# -*- coding: utf-8 -*-

import pytest

from ..einstein_smoluchowski import electrical_mobility_from_D
from ..units import default_units, allclose, default_constants, to_unitless
from ..util.testing import requires
from ..units import units_library


def test_electrical_mobility_from_D():
    D = 3
    z = -2
    T = 100
    mu = electrical_mobility_from_D(D, z, T)
    e = 1.60217657e-19
    kB = 1.3806488e-23
    ref = z * e * D / (kB * T)
    assert allclose(mu, ref, rtol=1e-5)

    mu2 = electrical_mobility_from_D(3, -2, 100)
    assert allclose(mu2, -2 * 1.60217657e-19 * 3 / 1.3806488e-23 / 100, rtol=1e-5)


@requires(units_library)
@pytest.mark.parametrize("constants", [None, default_constants])
def test_electrical_mobility_from_D__with_units(constants):
    metre = default_units.metre
    second = default_units.second
    kelvin = default_units.kelvin
    coulomb = default_units.coulomb
    kilogram = default_units.kilogram

    D = 3 * metre**2 / second
    z = -2
    T = 100 * kelvin
    mu = electrical_mobility_from_D(D, z, T, constants, default_units)
    e = 1.60217657e-19 * coulomb
    kB = 1.3806488e-23 * metre**2 * kilogram * second**-2 / kelvin
    ref = z * e * D / (kB * T)
    assert allclose(mu, ref, rtol=1e-5)
    assert allclose(
        to_unitless(mu, metre**2 / default_units.volt / second),
        electrical_mobility_from_D(3, z, 100),
        rtol=1e-5,
    )


@requires(units_library)
def test_electrical_mobility_from_D__issue_247():
    u = default_units
    D, T = 2.03e-9 * u.m**2 / u.s, 298.15 * u.K
    mu = electrical_mobility_from_D(D, 1, T, units=u)
    assert abs(to_unitless(mu, u.m**2 / u.V / u.s) - 7.9011e-8) < 1e-11
