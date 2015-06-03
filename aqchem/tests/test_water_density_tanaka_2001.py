import warnings

from ..water_density_tanaka_2001 import water_density


def test_water_density():
    warnings.filterwarnings("error")
    assert abs(water_density(273.15 + 0) - 999.8395) < 0.004
    assert abs(water_density(273.15 + 4) - 999.9720) < 0.003
    assert abs(water_density(273.15 + 10) - 999.7026) < 0.0003
    assert abs(water_density(273.15 + 15) - 999.1026) < 0.0001
    assert abs(water_density(273.15 + 20) - 998.2071) < 0.0005
    assert abs(water_density(273.15 + 22) - 997.7735) < 0.0007
    assert abs(water_density(273.15 + 25) - 997.0479) < 0.0009
    assert abs(water_density(273.15 + 30) - 995.6502) < 0.0016
    assert abs(water_density(273.15 + 40) - 992.2) < 0.02

    try:
        water_density(1)
    except UserWarning:
        pass  # good warning raised
    else:
        raise
    warnings.resetwarnings()

    try:
        import quantities as pq
        import numpy as np
        assert np.allclose(water_density(298.15*pq.K, units=pq),
                           997.047021671824*pq.kg/pq.m**3)
        assert np.allclose(water_density(np.linspace(297, 299)*pq.K, units=pq),
                           997*pq.kg/pq.m**3, rtol=1e-3, atol=1e-3)
    except ImportError:
        pass
