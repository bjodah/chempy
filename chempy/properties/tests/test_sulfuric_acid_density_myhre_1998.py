# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from ..sulfuric_acid_density_myhre_1998 import (
    sulfuric_acid_density, density_from_concentration
)


def test_sulfuric_acid_density():
    rho = sulfuric_acid_density(.1, 298)
    assert abs(1063.8 - rho) < .1


def test_density_from_concentration():
    rho = density_from_concentration(1000)
    assert abs(1058.5 - rho) < .1


def test_density_from_concentration__units():
    from chempy.units import default_units as units
    rho = density_from_concentration(0.4*units.molar, units=units)
    assert abs(1.024*units.kg/units.decimetre**3/rho - 1) < 1e-3
