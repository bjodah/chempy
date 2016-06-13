# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import pytest

from ..periodic import atomic_number, mass_from_composition, relative_atomic_masses, groups
from ..testing import requires
from ..parsing import formula_to_composition, parsing_library


def test_atomic_number():
    assert atomic_number('U') == 92
    assert atomic_number('carbon') == 6
    assert atomic_number('ununpentium') == 115
    with pytest.raises(ValueError):
        atomic_number('unobtainium')


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
    mass = mass_from_composition(formula_to_composition('NaF'))
    assert abs(41.988172443 - mass) < 1e-7

    Fminus = mass_from_composition(formula_to_composition('F/-'))
    assert abs(Fminus - 18.998403163 - 5.489e-4) < 1e-7
