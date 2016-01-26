from ..einstein_smoluchowski import electrical_mobility_from_D
from ..units import default_units, allclose, default_constants


def test_electrical_mobility_from_D():
    metre = default_units.metre
    second = default_units.second
    kelvin = default_units.kelvin
    coulomb = default_units.coulomb
    kilogram = default_units.kilogram

    D = 3*metre/second**2
    z = -2
    T = 100*kelvin
    mu = electrical_mobility_from_D(D, z, T, default_constants, default_units)
    e = 1.60217657e-19 * coulomb
    kB = 1.3806488e-23 * metre**2 * kilogram * second**-2 / kelvin
    ref = z*e*D/(kB*T)
    assert allclose(mu, ref, rtol=1e-5)
