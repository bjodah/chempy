# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


from ..chemistry import Substance, Solute, Reaction


def test_Substance():
    s = Substance('Hp', formula='H{+}')
    assert s.composition == {0: 1, 1: 1}
    assert s.charge == 1
    assert abs(s.mass - 1.008) < 1e-3


def test_Solute():
    s = Solute('Hp', formula='H{+}', solid=True)
    assert abs(s.mass - 1.00794 + 5.5e-4) < 2e-5


def test_Reaction():
    substances = s_Hp, s_OHm, s_H2O = (
        Solute('H+', composition={0: 1, 1: 1}),
        Solute('OH-', composition={0: -1, 1: 1, 8: 1}),
        Solute('H2O', composition={0: 0, 1: 2, 8: 1}),
    )
    substance_names = Hp, OHm, H2O = [s.name for s in substances]
    substance_dict = {n: s for n, s in zip(substance_names, substances)}
    r1 = Reaction({Hp: 1, OHm: 1}, {H2O: 1})
    assert sum(r1.composition_violation(substance_dict)) == 0
    assert r1.charge_neutrality_violation(substance_dict) == 0

    r2 = Reaction({Hp: 1, OHm: 1}, {H2O: 2})
    assert sum(r2.composition_violation(substance_dict)) != 0
    assert r2.charge_neutrality_violation(substance_dict) == 0

    r3 = Reaction({Hp: 2, OHm: 1}, {H2O: 2})
    assert sum(r3.composition_violation(substance_dict)) != 0
    assert r3.charge_neutrality_violation(substance_dict) != 0
