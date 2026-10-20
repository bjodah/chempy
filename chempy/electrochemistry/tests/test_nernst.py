# -*- coding: utf-8 -*-

import math

import pytest

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
    v = nernst_potential(145, 15, 1, 310 * K, default_constants)
    assert allclose(1000 * v, 60.605 * J / coulomb, rtol=1e-4)

    v2 = nernst_potential(145, 15, 1, 310 * K, units=default_units)
    assert allclose(1000 * v2, 60.605 * J / coulomb, rtol=1e-4)


@requires(units_library)
@pytest.mark.parametrize("constants", [None, default_constants])
@pytest.mark.parametrize("backend", ["math", "numpy"])
def test_nernst_potential__mixed_concentration_units(constants, backend):
    u = default_units
    kwargs = {"units": u} if constants is None else {"constants": constants}
    result = nernst_potential(
        145 * u.millimolar,
        0.015 * u.molar,
        1,
        310 * u.K,
        backend=__import__(backend),
        **kwargs
    )
    assert allclose(result, 60.605 * u.mV, rtol=1e-4)


@requires(units_library)
def test_nernst_potential__incompatible_concentrations():
    u = default_units
    with pytest.raises(ValueError):
        nernst_potential(145 * u.m, 15 * u.s, 1, 310 * u.K, units=u)


@requires("sympy")
def test_nernst_potential__symbolic():
    import sympy as sp

    outside, inside, temperature = sp.symbols("outside inside temperature")
    result = nernst_potential(outside, inside, 1, temperature, backend=sp)
    reference = 8.3144598 / 96485.33289 * temperature * sp.log(outside / inside)
    assert sp.simplify(result - reference) == 0
    assert math.isfinite(
        float(result.subs({outside: 145, inside: 15, temperature: 310}))
    )
