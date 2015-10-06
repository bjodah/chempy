from math import log as ln

from ..debye_huckel import A as A_dh, B as B_dh
from ..debye_huckel import limiting_log_gamma
from .util import allclose

import quantities as pq


def test_A():
    A20 = A_dh(80.1, 293.15, 998.2071)/ln(10)
    assert abs(A20 - 0.50669) < 1e-5

    A20q = A_dh(80.1, 293.15*pq.K, 998.2071*pq.kg/pq.m**3,
                b0=1*pq.mol/pq.kg, constants=pq.constants)
    assert abs(A20q - 0.50669) < 1e-5


def test_B():
    assert abs(1e-10*B_dh(80.1, 293.15, 998.2071, 1) - 0.3282) < 1e-3

    B20q = B_dh(80.1, 293.15*pq.K, 998.2071*pq.kg/pq.m**3,
                b0=pq.mol/pq.kg, constants=pq.constants, units=pq)
    close = allclose(B20q.simplified, 0.3282/pq.angstrom, rtol=1e-3)
    assert close


def test_limiting_log_gamma():
    A20 = A_dh(80.1, 293.15, 998.2071)/ln(10)
    log_gamma = limiting_log_gamma(0.4, -3, A20)
    assert abs(log_gamma + 2.884130) < 1e-4
