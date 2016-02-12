# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


_radii_nm = {  # in nanometer
    'Fe3+': 0.9,
    'SCN-': 0.35,
    'Fe2+': 0.6,
}


def get_radii(key, units=None):
    if units is None:
        return _radii_nm[key]
    else:
        return _radii_nm[key]*units.nm
