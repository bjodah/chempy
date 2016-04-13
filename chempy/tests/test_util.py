# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


from .._util import prodpow
from ..util.testing import requires


def test_prodpow():
    result = prodpow([11, 13], [[0, 1], [1, 2]])
    assert result[0] == 13
    assert result[1] == 11*13*13


@requires('sympy')
def test_prodpow__symbols():
    import sympy
    a, b = sympy.symbols('a b')
    exprs = prodpow([a, b], [[0, 1], [1, 2]])
    assert exprs[0] == b
    assert exprs[1] == a*b**2
