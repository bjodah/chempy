# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

try:
    import numpy
    import scipy
except ImportError:
    numpy = None
    scipy = None
else:
    from ..ethanol_water_density_lange import (
        aqueous_ethanol_density
    )

from chempy.util.testing import requires


@requires('numpy', 'scipy')
def test_aqueous_ethanol_density():
    rho = aqueous_ethanol_density(.01, 283.15)
    assert abs(997.85 - rho) < .01
