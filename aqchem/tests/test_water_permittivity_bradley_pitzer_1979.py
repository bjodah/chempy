import warnings

from ..water_permittivity_bradley_pitzer_1979 import water_permittivity


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

    try:
        import quantities as pq
        import numpy as np
        assert np.allclose(water_permittivity(
            298.15*pq.K, 1*pq.bar,
            units=pq), 78.38436874203077)
    except ImportError:
        pass
