import warnings

from chempy.units import allclose

from ..water_permittivity_bradley_pitzer_1979 import water_permittivity
from chempy.util.testing import requires
from chempy.units import linspace, units_library, default_units as u


def test_water_permittivity():
    warnings.filterwarnings("error")
    abs(water_permittivity(273.15 + 0) - 80) < 1.0
    abs(water_permittivity(273.15 + 20) - 80.1) < 0.2
    abs(water_permittivity(273.15 + 100) - 55.3) < 0.5

    try:
        water_permittivity(1)
    except UserWarning:
        pass  # good: warning raised
    else:
        raise
    warnings.resetwarnings()


@requires(units_library)
def test_water_permittivity__units():
    assert allclose(water_permittivity(
        298.15*u.K, 1*u.bar,
        units=u), 78.38436874203077)
    assert allclose(water_permittivity(
        linspace(297.5, 298.65)*u.K, 1*u.bar,
        units=u), 78, rtol=1e-2, atol=1e-2)
