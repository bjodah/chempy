# -*- coding: utf-8 -*-

import numpy as np


def allclose(a, b, rtol=1e-8, atol=None):
    d = np.abs(a - b)
    lim = np.abs(a)*rtol
    if atol is not None:
        lim += atol
    return np.all(d < lim)
