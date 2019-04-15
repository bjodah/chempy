import warnings
from chempy.util.testing import requires
from chempy.units import units_library
from ..water_diffusivity_holz_2000 import water_self_diffusion_coefficient as w_sd


def test_water_self_diffusion_coefficient():
    warnings.filterwarnings("error")
    assert abs(w_sd(273.15 + 0.) - 1.099e-9) < 0.027e-9
    assert abs(w_sd(273.15 + 4.) - 1.261e-9) < 0.011e-9
    assert abs(w_sd(273.15 + 10) - 1.525e-9) < 0.007e-9
    assert abs(w_sd(273.15 + 15) - 1.765e-9) < 0.006e-9
    assert abs(w_sd(273.15 + 20) - 2.023e-9) < 0.001e-9
    assert abs(w_sd(273.15 + 25) - 2.299e-9) < 0.001e-9
    assert abs(w_sd(273.15 + 30) - 2.594e-9) < 0.001e-9
    assert abs(w_sd(273.15 + 35) - 2.907e-9) < 0.004e-9

    try:
        w_sd(1)
    except UserWarning:
        pass  # good warning raised
    else:
        raise
    warnings.resetwarnings()


@requires(units_library)
def test_water_self_diffusion_coefficient__units():
    from chempy.units import allclose, linspace, default_units as u
    unit = u.m**2/u.s
    assert allclose(1e9*w_sd(298.15*u.K, units=u),
                    2.299*unit, rtol=1e-3, atol=1e-8*unit)
    assert allclose(1e9*w_sd(linspace(297, 299)*u.K, units=u),
                    2.299*u.m**2/u.s, rtol=5e-2, atol=1e-2*unit)
