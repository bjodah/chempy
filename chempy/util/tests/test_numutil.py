# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from ..testing import requires
from ..numutil import broadcast_stack


@requires('numpy')
def test_broadcast_stack():
    import numpy as np
    a = np.arange(9.)
    b = broadcast_stack(np.arange(5.), np.arange(5., 9.))
    assert a.shape == b.shape and np.all(a == b)

    c = broadcast_stack([0, 1, 2], [3], 4, as_scalars=True)
    d = np.array([
        [0, 3, 4],
        [1, 3, 4],
        [2, 3, 4]
    ])
    assert c.shape == d.shape and np.all(c == d)
