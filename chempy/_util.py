# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

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
