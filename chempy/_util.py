# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from functools import reduce
from operator import mul, add

try:
    import numpy as np
    from numpy import any as _any

    def prodpow(bases, exponents):
        """
        Examples
        --------
        >>> prodpow([2, 3], np.array([[0, 1], [1, 2]]))
        array([ 3, 18])

        """
        exponents = np.asarray(exponents)
        return np.multiply.reduce(bases**exponents, axis=-1)

except ImportError:  # no NumPy available
    def _any(arg):
        if arg is True:
            return True
        if arg is False:
            return False
        return any(arg)

    def prodpow(bases, exponents):
        """
        Examples
        --------
        >>> prodpow([2, 3], [[0, 1], [1, 2]])
        [3, 18]

        """
        result = []
        for row in exponents:
            res = 1
            for b, e in zip(bases, row):
                res *= b**e
            result.append(res)
        return result


def get_backend(backend):
    if isinstance(backend, str):
        backend = __import__(backend)
    if backend is None:
        try:
            import numpy as backend
        except ImportError:
            import math as backend
    return backend


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


def reducemap(args, reduce_op, map_op):
    return reduce(reduce_op, map(map_op, *args))


def vec_dot_vec(vec1, vec2):
    # return np.dot(vec1, vec2)
    # return np.add.reduce(np.multiply(vec1, vec2))
    return reducemap((vec1, vec2), add, mul)


def mat_dot_vec(iter_mat, iter_vec, iter_term=None):  # pure python (slow)
    if iter_term is None:
        return [vec_dot_vec(row, iter_vec) for row in iter_mat]
    else:
        # daxpy
        return [vec_dot_vec(row, iter_vec) + term for row, term
                in zip(iter_mat, iter_term)]


# def composition_balance(substances, concs, composition_number):
#     if not hasattr(concs, 'ndim') or concs.ndim == 1:
#         res = 0
#     elif concs.ndim == 2:
#         res = np.zeros(concs.shape[0])
#         concs = concs.T
#     else:
#         raise NotImplementedError
#     for s, c in zip(substances, concs):
#         res += s.composition.get(composition_number, 0)*c
#     return res
