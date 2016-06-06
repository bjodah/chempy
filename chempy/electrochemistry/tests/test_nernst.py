# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from ..nernst import nernst_potential


def test_nernst_potential():
    # Sodium in cells
    assert abs(1000 * nernst_potential(145, 15, 1, 310) - 60.605) < 1e-4
    # Potassium in cells
    assert abs(1000 * nernst_potential(4, 150, 1, 310) - (-96.8196)) < 1e-4
