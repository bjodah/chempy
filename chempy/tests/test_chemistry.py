# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


from ..chemistry import Substance, Solute


def test_Substance():
    s = Substance('Hp', formula='H{+}')
    assert s.composition == {0: -1, 1: 1}


def test_Solute():
    s = Solute('Hp', formula='H{+}', solid=True)
    assert abs(s.mass - 1.00794 + 5.5e-4) < 2e-5
