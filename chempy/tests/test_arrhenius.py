# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function)

from ..arrhenius import arrhenius_equation


def test_arrhenius_equation():
    assert abs(arrhenius_equation(3, 831.4472, 100) - 3/2.7182818) < 1e-7
