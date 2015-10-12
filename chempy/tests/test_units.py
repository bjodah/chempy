import numpy as np
from ..units import allclose, default_units

second = default_units.second
hour = default_units.hour


def test_allclose():
    a = np.linspace(2, 3)*second
    b = np.linspace(2/3600., 3/3600.)*hour
    assert allclose(a, b)
