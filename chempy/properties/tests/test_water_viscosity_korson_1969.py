import warnings

import pytest

from ..water_viscosity_korson_1969 import water_viscosity
from chempy.util.testing import requires
from chempy.units import units_library, default_units as u, allclose, to_unitless


def test_water_viscosity():
    warnings.filterwarnings("error")  # Table II (p. 38):
    assert abs(water_viscosity(273.15 + 0) - 1.7916) < 5e-4
    assert abs(water_viscosity(273.15 + 5) - 1.5192) < 5e-4
    assert abs(water_viscosity(273.15 + 10) - 1.3069) < 5e-4
    assert abs(water_viscosity(273.15 + 15) - 1.1382) < 5e-4
    assert abs(water_viscosity(273.15 + 20) - 1.0020) < 5e-4
    assert abs(water_viscosity(273.15 + 25) - 0.8903) < 5e-4
    assert abs(water_viscosity(273.15 + 30) - 0.7975) < 5e-4
    assert abs(water_viscosity(273.15 + 35) - 0.7195) < 5e-4
    assert abs(water_viscosity(273.15 + 40) - 0.6532) < 5e-4
    assert abs(water_viscosity(273.15 + 45) - 0.5963) < 5e-4
    assert abs(water_viscosity(273.15 + 50) - 0.5471) < 5e-4
    assert abs(water_viscosity(273.15 + 55) - 0.5042) < 5e-4
    assert abs(water_viscosity(273.15 + 60) - 0.4666) < 5e-4
    assert abs(water_viscosity(273.15 + 65) - 0.4334) < 5e-4
    assert abs(water_viscosity(273.15 + 70) - 0.4039) < 5e-4
    assert abs(water_viscosity(273.15 + 75) - 0.3775) < 5e-4
    assert abs(water_viscosity(273.15 + 80) - 0.3538) < 5e-4
    assert abs(water_viscosity(273.15 + 85) - 0.3323) < 5e-4
    assert abs(water_viscosity(273.15 + 90) - 0.3128) < 5e-4
    assert abs(water_viscosity(273.15 + 95) - 0.2949) < 6e-4
    assert abs(water_viscosity(273.15 + 100) - 0.2783) < 2e-3
    warnings.resetwarnings()


@requires(units_library)
@pytest.mark.parametrize("temperature_unit", ["K", "mK"])
def test_water_viscosity__units(temperature_unit):
    import numpy as np

    T = (298.15 * u.K).rescale(getattr(u, temperature_unit))
    eta = water_viscosity(T, units=u)
    assert allclose(eta, 0.8903 * u.centipoise, rtol=1e-3)
    assert abs(to_unitless(eta, u.centipoise) - water_viscosity(298.15)) < 1e-12

    T_arr = (np.array([273.15, 293.15, 298.15]) * u.K).rescale(
        getattr(u, temperature_unit)
    )
    eta_arr = water_viscosity(T_arr, units=u)
    assert allclose(
        eta_arr,
        np.array([1.7916, 1.0020, 0.8903]) * u.centipoise,
        rtol=1e-3,
    )


@requires(units_library)
def test_water_viscosity__defaults_and_reference_viscosity():
    assert allclose(water_viscosity(units=u), water_viscosity() * u.centipoise)
    eta20 = 2.004e-3 * u.Pa * u.s
    assert allclose(water_viscosity(293.15 * u.K, eta20=eta20, units=u), eta20)
    assert allclose(
        water_viscosity(298.15 * u.K, eta20=eta20, units=u),
        2 * water_viscosity(298.15) * u.centipoise,
    )


@requires(units_library)
def test_water_viscosity__invalid_units():
    with pytest.raises(ValueError):
        water_viscosity(298.15 * u.s, units=u)


@requires(units_library)
@pytest.mark.parametrize("T", [272.15, 374.15])
def test_water_viscosity__units_temperature_warning(T):
    with pytest.warns(UserWarning, match="outside range"):
        water_viscosity(T * u.K, units=u)
