# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)
import re


try:
    import numpy as np
except ImportError:
    np = None


# See https://github.com/python-quantities/python-quantities/pull/112
def format_units_html(udict, font='%s', mult=r'&sdot;', paren=False):
    '''
    Replace the units string provided with an equivalent html string.

    Exponentiation (m**2) will be replaced with superscripts (m<sup>2</sup>})

    No formating is done, change `font` argument to e.g.:
    '<span style="color: #0000a0">%s</span>' to have text be colored blue.

    Multiplication (*) are replaced with the symbol specified by the mult
    argument. By default this is the latex &sdot; symbol.  Other useful options
    may be '' or '*'.

    If paren=True, encapsulate the string in '(' and ')'

    '''
    from quantities.markup import format_units
    res = format_units(udict)
    if res.startswith('(') and res.endswith(')'):
        # Compound Unit
        compound = True
    else:
        # Not a compound unit
        compound = False
    # Replace exponentiation (**exp) with ^{exp}
    res = re.sub(r'\*{2,2}(?P<exp>\d+)', r'<sup>\g<exp></sup>', res)
    # Remove multiplication signs
    res = re.sub(r'\*', mult, res)
    if paren and not compound:
        res = '(%s)' % res
    res = font % res
    return res


def _patch_pow0_py35(pq):
    try:
        pq.metre**0
    except Exception:
        pq.quantity.Quantity.__pow__ = pq.quantity.check_uniform(lambda self, other: np.power(self, other))


def _patch_quantities(pq):
    # See https://github.com/python-quantities/python-quantities/pull/112
    if not hasattr(pq.dimensionality.Dimensionality, 'html'):
        pq.dimensionality.Dimensionality.html = property(lambda self: format_units_html(self))

    # See https://github.com/python-quantities/python-quantities/pull/116
    a = pq.UncertainQuantity([1, 2], pq.m, [.1, .2])
    if (-a).uncertainty[0] != a.uncertainty[0]:
        pq.UncertainQuantity.__neg__ = lambda self: self*-1
    a = pq.UncertainQuantity([1, 2], pq.m, [.1, .2])
    assert (-a).uncertainty[0] == (a*-1).uncertainty[0]

    # See https://github.com/python-quantities/python-quantities/pull/126
    _patch_pow0_py35(pq)
    assert (3*pq.m)**0 == 1*pq.dimensionless
