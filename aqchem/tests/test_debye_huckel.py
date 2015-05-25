from math import log as ln

from ..debye_huckel import A as A_dh, B as B_dh

import quantities as pq


def test_A():
    A20 = A_dh(80.1, 293.15, 998.2071)/ln(10)
    assert abs(A20 - 0.50669) < 1e-5

    A20q = A_dh(80.1, 293.15*pq.K, 998.2071*pq.kg/pq.m**3,
                b0=1*pq.mol/pq.kg, constants=pq.constants)
    assert abs(A20q - 0.50669) < 1e-5


def test_B():
    assert abs(B_dh(80.1, 293.15) - 0.3282) < 1e-3
