# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from ..nernst import nernst_potential
from chempy.util.testing import requires
from chempy.units import default_units, default_constants, units_library, allclose


def test_nernst_potential():
    """
    Test cases obtained from textbook examples of Nernst potential in cellular
    membranes. 310K = 37C, typical mammalian cell environment temperature.
    """
    # Sodium in cells
    assert abs(1000 * nernst_potential(145, 15, 1, 310) - 60.605) < 1e-4
    # Potassium in cells
    assert abs(1000 * nernst_potential(4, 150, 1, 310) - (-96.8196)) < 1e-4
    # Calcium in cells
    assert abs(1000 * nernst_potential(2, 7e-5, 2, 310) - 137.0436) < 1e-4
    # Chloride in cells
    assert abs(1000 * nernst_potential(110, 10, -1, 310) - (-64.0567)) < 1e-4


@requires(units_library)
def test_nernst_potential__units():
    J = default_units.joule
    K = default_units.kelvin
    coulomb = default_units.coulomb
    v = nernst_potential(145, 15, 1, 310*K, default_constants)
    assert allclose(1000 * v, 60.605*J/coulomb, rtol=1e-4)
