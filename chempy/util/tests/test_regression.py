# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from ..regression import least_squares, irls
from ..testing import requires


@requires('numpy')
def test_irls():
    import numpy as np
    x = np.linspace(0, 100)
    y = np.exp(-x/47)
    b, c, info = irls(x, np.log(y))
    assert abs(b[1] + 1/47) < 1e-5
    assert np.all(c < 1e-4)
    assert info['success'] is True
    assert info['niter'] < 3


@requires('numpy')
def test_least_squares():
    import numpy as np
    x, y, w = [0, 1, 2], [0, 1, 2], [1, 1, 1]
    beta, vcv, r2 = least_squares(x, y)
    assert np.allclose(beta, [0, 1])
    assert np.allclose(vcv, 0)
    assert np.allclose(r2, 1)

    wbeta, wvcv, wr2 = least_squares(x, y, w)
    assert np.allclose(wbeta, [0, 1])
    assert np.allclose(wvcv, 0)
    assert np.allclose(wr2, 1)
