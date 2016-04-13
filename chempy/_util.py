# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

try:
    import numpy as np

    def prodpow(bases, exponents):
        """
        Examples
        --------
        >>> prodpow([2, 3], np.array([[0, 1], [1, 2]]))
        array([ 3, 18])

        """
        exponents = np.asarray(exponents)
        return np.multiply.reduce(bases**exponents, axis=-1)

except ImportError:

    def prodpow(bases, exponents):
        """
        Examples
        --------
        >>> prodpow([2, 3], [[0, 1], [1, 2]])
        [3, 8]

        """
        result = []
        for row in exponents:
            res = 1
            for b, e in zip(bases, row):
                res *= b**e
            result.append(res)
        return result


def intdiv(p, q):
    """ Integer divsions which rounds toward zero

    Examples
    --------
    >>> intdiv(3, 2)
    1
    >>> intdiv(-3, 2)
    -1
    >>> -3 // 2
    -2

    """
    r = p // q
    if r < 0 and q*r != p:
        r += 1
    return r
