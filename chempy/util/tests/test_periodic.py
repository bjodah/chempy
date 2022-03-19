# -*- coding: utf-8 -*-

import pytest

from ..periodic import (
    atomic_number,
    mass_from_composition,
    relative_atomic_masses,
    groups,
    symbols,
)
from ..testing import requires
from ..parsing import formula_to_composition, parsing_library


def test_atomic_number():
    assert atomic_number("U") == 92
    assert atomic_number("u") == 92
    assert atomic_number("carbon") == 6
    assert atomic_number("oganesson") == 118
    with pytest.raises(ValueError):
        atomic_number("unobtainium")


def test_mass_from_composition():
    mass = mass_from_composition({11: 1, 9: 1})
    assert abs(41.988172443 - mass) < 1e-7


def test_relative_atomic_masses():
    assert relative_atomic_masses[0] == 1.008


def test_groups():
    assert groups[1] == (1, 3, 11, 19, 37, 55, 87)
    assert groups[2] == (4, 12, 20, 38, 56, 88)
    assert groups[13] == (5, 13, 31, 49, 81, 113)
    assert groups[14] == (6, 14, 32, 50, 82, 114)
    assert groups[15] == (7, 15, 33, 51, 83, 115)
    assert groups[16] == (8, 16, 34, 52, 84, 116)
    assert groups[17] == (9, 17, 35, 53, 85, 117)
    assert groups[18] == (2, 10, 18, 36, 54, 86, 118)


@requires(parsing_library)
def test_mass_from_composition__formula():
    mass = mass_from_composition(formula_to_composition("NaF"))
    assert abs(41.988172443 - mass) < 1e-7

    Fminus = mass_from_composition(formula_to_composition("F-"))
    assert abs(Fminus - 18.998403163 - 5.489e-4) < 1e-7


def test_molar_masses_of_the_elements_gh172():
    from chempy import Substance

    previous_mass = 0
    for symbol in symbols:
        this_mass = Substance.from_formula(symbol).mass
        # Note that any average atomic masses of naturally occurring isotopes loose their
        # meaning for trans-uranium elements. The numbers are almost to be considered arbitrary
        # and the user needs to know what isotope they have at hand (and use isotopic mass).
        if symbol in ("K", "Ni", "I", "Pa", "Np", "Am"):
            assert this_mass < previous_mass
        elif symbol in ("Bk", "Og"):
            assert this_mass == previous_mass
        else:
            assert this_mass > previous_mass
        previous_mass = this_mass

    assert Substance.from_formula("Hs").mass == 271
