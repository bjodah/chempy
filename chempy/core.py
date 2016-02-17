# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function


def ionic_strength(molalities, charges, b0=1):
    """ Calculates the ionic strength

    Parameters
    ----------
    molalities: float
        optionally with unit (amount / mass)
    charges: iterable of integers
        charge of respective ion
    b0: float
        reference molality, optionally with unit (amount / mass)
        by IUPAC defines it as 1 mol/kg. (default: 1)

    Examples
    --------
    >>> ionic_strength([1e-3, 3e-3], [3, -1]) == .5 * (9 + 3) * 1e-3
    True

    """
    tot = 0
    if len(molalities) != len(charges):
        raise ValueError("molalities and charges of different lengths")
    for b, z in zip(molalities, charges):
        tot += (b/b0)*z**2
    return tot/2


class _ActivityProductBase(object):
    """ Baseclass for activity products """

    def __init__(self, stoich, *args):
        self.stoich = stoich
        self.args = args

    def __call__(self, c):
        pass
