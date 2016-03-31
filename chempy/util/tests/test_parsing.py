# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import pytest

from ..parsing import (
    to_composition, relative_atomic_masses, mass_from_composition,
    to_latex, to_reaction, atomic_number
)


def test_atomic_number():
    assert atomic_number('U') == 92
    assert atomic_number('carbon') == 6
    assert atomic_number('ununpentium') == 115
    with pytest.raises(ValueError):
        atomic_number('unobtainium')


def test_to_composition():
    assert to_composition('H2O') == {1: 2, 8: 1}
    assert to_composition('Fe/3+') == {0: 3, 26: 1}
    assert to_composition('Fe+3') == {0: 3, 26: 1}
    assert to_composition('Na/+') == {0: 1, 11: 1}
    assert to_composition('Na+1') == {0: 1, 11: 1}
    assert to_composition('Na+') == {0: 1, 11: 1}
    assert to_composition('Cl/-') == {0: -1, 17: 1}
    assert to_composition('Cl-') == {0: -1, 17: 1}
    assert to_composition('NaCl') == {11: 1, 17: 1}
    assert to_composition('NaCl(s)') == {11: 1, 17: 1}
    assert to_composition('Fe(SCN)2/+') == {0: 1, 6: 2, 7: 2, 16: 2, 26: 1}
    assert to_composition('Fe(SCN)2+') == {0: 1, 6: 2, 7: 2, 16: 2, 26: 1}
    assert to_composition('Fe(SCN)2+1') == {0: 1, 6: 2, 7: 2, 16: 2, 26: 1}
    assert to_composition('((H2O)2OH)12') == {1: 60, 8: 36}

    # Special case: solvated electron:
    assert to_composition('e-') == {0: -1}
    assert to_composition('e-1') == {0: -1}
    assert to_composition('e-(aq)') == {0: -1}

    # prefixes and suffixes
    assert to_composition('.NO2(g)') == {7: 1, 8: 2}
    assert to_composition('.NH2') == {1: 2, 7: 1}
    assert to_composition('ONOOH') == {1: 1, 7: 1, 8: 3}
    assert to_composition('.ONOO') == {7: 1, 8: 3}
    assert to_composition('.NO3/2-') == {0: -2, 7: 1, 8: 3}
    assert to_composition('.NO3-2') == {0: -2, 7: 1, 8: 3}

    with pytest.raises(ValueError):
        to_composition('F-F')

    assert to_composition('alpha-FeOOH(s)') == {1: 1, 8: 2, 26: 1}
    assert to_composition('epsilon-Zn(OH)2(s)') == {1: 2, 8: 2, 30: 1}


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
    assert to_latex('Fe(CN)6-3') == 'Fe(CN)_{6}^{3-}'
    assert to_latex('C18H38/2+') == 'C_{18}H_{38}^{2+}'
    assert to_latex('C18H38/+2') == 'C_{18}H_{38}^{2+}'
    assert to_latex('C18H38+2') == 'C_{18}H_{38}^{2+}'
    assert to_latex('((H2O)2OH)12') == '((H_{2}O)_{2}OH)_{12}'
    assert to_latex('NaCl') == 'NaCl'
    assert to_latex('NaCl(s)') == 'NaCl(s)'
    assert to_latex('e-(aq)') == 'e^{-}(aq)'
    assert to_latex('.NO2(g)') == r'^\bullet NO_{2}(g)'
    assert to_latex('.NH2') == r'^\bullet NH_{2}'
    assert to_latex('ONOOH') == 'ONOOH'
    assert to_latex('.ONOO') == r'^\bullet ONOO'
    assert to_latex('.NO3/2-') == r'^\bullet NO_{3}^{2-}'
    assert to_latex('.NO3-2') == r'^\bullet NO_{3}^{2-}'
    assert to_latex('alpha-FeOOH(s)') == r'\alpha-FeOOH(s)'
    assert to_latex('epsilon-Zn(OH)2(s)') == r'\varepsilon-Zn(OH)_{2}(s)'


def test_to_reaction():
    from chempy.chemistry import Reaction, Equilibrium
    rxn = to_reaction(
        "H+ + OH- -> H2O; 1.4e11; ref={'doi': '10.1039/FT9908601539'}",
        'H+ OH- H2O'.split(), '->', Reaction)
    assert rxn.__class__ == Reaction

    assert rxn.reac['H+'] == 1
    assert rxn.reac['OH-'] == 1
    assert rxn.prod['H2O'] == 1
    assert rxn.param == 1.4e11
    assert rxn.ref['doi'].startswith('10.')

    eq = to_reaction("H+ + OH- = H2O; 1e-14; ref='rt, [H2O] == 1 M'",
                     'H+ OH- H2O'.split(), '=', Equilibrium)
    assert eq.__class__ == Equilibrium

    assert eq.reac['H+'] == 1
    assert eq.reac['OH-'] == 1
    assert eq.prod['H2O'] == 1
    assert eq.ref.startswith('rt')

    for s in ['2 e-(aq) + (2 H2O) -> H2 + 2 OH- ; 1e6 ; ',
              '2 * e-(aq) + (2 H2O) -> 1 * H2 + 2 * OH- ; 1e6 ; ']:
        rxn2 = to_reaction(s, 'e-(aq) H2 OH- H2O'.split(), '->', Reaction)
        assert rxn2.__class__ == Reaction
        assert rxn2.reac['e-(aq)'] == 2
        assert rxn2.inact_reac['H2O'] == 2
        assert rxn2.prod['H2'] == 1
        assert rxn2.prod['OH-'] == 2
        assert rxn2.param == 1e6
