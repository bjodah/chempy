# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import sympy

from .._util import prodpow


def test_prodpow():
    a, b = sympy.symbols('a b')
    exprs = prodpow([a, b], [[0, 1], [1, 2]])
    assert exprs[0] == b
    assert exprs[1] == a*b**2
