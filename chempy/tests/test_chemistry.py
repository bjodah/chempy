# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


from ..chemistry import Substance, Solute, Reaction


def test_Substance():
    s = Substance('Hp', formula='H{+}')
    assert s.composition == {0: -1, 1: 1}


def test_Solute():
    s = Solute('Hp', formula='H{+}', solid=True)
    assert abs(s.mass - 1.00794 + 5.5e-4) < 2e-5


def test_Reaction():
    substances = Hp, OHm, H2O = (
        Solute('H+', charge=1, mass=1),
        Solute('OH-', charge=-1, mass=17),
        Solute('H2O', charge=0, mass=18),
    )
    r1 = Reaction({Hp: 1, OHm: 1}, {H2O: 1})
    assert r1.mass_balance_violation(substances) == 0
    assert r1.charge_neutrality_violation(substances) == 0

    r2 = Reaction({Hp: 1, OHm: 1}, {H2O: 2})
    assert r2.mass_balance_violation(substances) != 0
    assert r2.charge_neutrality_violation(substances) == 0

    r3 = Reaction({Hp: 2, OHm: 1}, {H2O: 2})
    assert r3.mass_balance_violation(substances) != 0
    assert r3.charge_neutrality_violation(substances) != 0
