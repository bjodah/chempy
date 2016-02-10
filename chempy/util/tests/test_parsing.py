# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import pytest

from ..parsing import (
    to_composition, relative_atomic_masses, mass_from_composition,
    to_latex
)


def test_to_composition():
    assert to_composition('H2O') == {1: 2, 8: 1}
    assert to_composition('Fe/3+') == {0: 3, 26: 1}
    assert to_composition('Na/+') == {0: 1, 11: 1}
    assert to_composition('Na+') == {0: 1, 11: 1}
    assert to_composition('Cl/-') == {0: -1, 17: 1}
    assert to_composition('Cl-') == {0: -1, 17: 1}
    assert to_composition('NaCl') == {11: 1, 17: 1}
    assert to_composition('NaCl(s)') == {11: 1, 17: 1}
    assert to_composition('Fe(SCN)2/+') == {0: 1, 6: 2, 7: 2, 16: 2, 26: 1}
    assert to_composition('Fe(SCN)2+') == {0: 1, 6: 2, 7: 2, 16: 2, 26: 1}
    assert to_composition('((H2O)2OH)12') == {1: 60, 8: 36}

    # Special case: solvated electron:
    assert to_composition('e-') == {0: -1}
    assert to_composition('e-(aq)') == {0: -1}

    # prefixes and suffixes
    assert to_composition('.NO2(g)') == {7: 1, 8: 2}
    assert to_composition('.NH2') == {1: 2, 7: 1}
    assert to_composition('ONOOH') == {1: 1, 7: 1, 8: 3}
    assert to_composition('.ONOO') == {7: 1, 8: 3}
    assert to_composition('.NO3/2-') == {0: -2, 7: 1, 8: 3}

    with pytest.raises(ValueError):
        to_composition('F-F')


def test_relative_atomic_masses():
    assert relative_atomic_masses[0] == 1.008


def test_mass_from_composition():
    mass = mass_from_composition(to_composition('NaF'))
    assert abs(41.988172443 - mass) < 1e-7

    Fminus = mass_from_composition(to_composition('F/-'))
    assert abs(Fminus - 18.998403163 - 5.489e-4) < 1e-7


def test_to_latex():
    assert to_latex('H2O') == 'H_{2}O'
    assert to_latex('C6H6/+') == 'C_{6}H_{6}^{+}'
    assert to_latex('Fe(CN)6/3-') == 'Fe(CN)_{6}^{3-}'
    assert to_latex('C18H38/2+') == 'C_{18}H_{38}^{2+}'
    assert to_latex('C18H38/+2') == 'C_{18}H_{38}^{2+}'
    assert to_latex('((H2O)2OH)12') == '((H_{2}O)_{2}OH)_{12}'
    assert to_latex('NaCl') == 'NaCl'
    assert to_latex('NaCl(s)') == 'NaCl(s)'
    assert to_latex('e-(aq)') == 'e^{-}(aq)'
    assert to_latex('.NO2(g)') == r'^\bullet NO_{2}(g)'
    assert to_latex('.NH2') == r'^\bullet NH_{2}'
    assert to_latex('ONOOH') == 'ONOOH'
    assert to_latex('.ONOO') == r'^\bullet ONOO'
    assert to_latex('.NO3/2-') == r'^\bullet NO_{3}^{2-}'
