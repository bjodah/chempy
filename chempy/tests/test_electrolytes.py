# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from math import log as ln

from ..electrolytes import A as A_dh, B as B_dh
from ..electrolytes import limiting_log_gamma, _ActivityProductBase, ionic_strength
from ..units import (allclose, units_library,
                     default_constants as consts,
                     default_units as u)
from ..util.testing import requires


def test_ionic_strength():
    assert abs(ionic_strength([0.1, 1.3, 2.1, 0.7],
                              [-1, 2, -3, 4], warn=False) - 17.7) < 1e-14


def test_ActivityProductBase():
    ap = _ActivityProductBase((1, -2, 3), 17, 42)
    assert ap.stoich == (1, -2, 3)
    assert ap.args == (17, 42)


def test_A():
    A20 = A_dh(80.1, 293.15, 998.2071)/ln(10)
    assert abs(A20 - 0.50669) < 1e-5


@requires(units_library)
def test_A__units():
    A20q = A_dh(80.1, 293.15*u.K, 998.2071*u.kg/u.m**3,
                b0=1*u.mol/u.kg, constants=consts)
    assert abs(A20q - 0.50669) < 1e-5


def test_B():
    assert abs(1e-10*B_dh(80.1, 293.15, 998.2071, 1) - 0.3282) < 1e-3


@requires(units_library)
def test_B__units():
    B20q = B_dh(80.1, 293.15*u.K, 998.2071*u.kg/u.m**3,
                b0=u.mol/u.kg, constants=consts, units=u)
    close = allclose(B20q.simplified, 0.3282/u.angstrom, rtol=1e-3)
    assert close


def test_limiting_log_gamma():
    A20 = A_dh(80.1, 293.15, 998.2071)/ln(10)
    log_gamma = limiting_log_gamma(0.4, -3, A20)
    assert abs(log_gamma + 2.884130) < 1e-4
