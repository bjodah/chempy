# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


_radii_nm = {  # in nanometer
    'Fe/3+': 0.9,
    'SCN/-': 0.35,
    'Fe/2+': 0.6,
}


def get_radii(key, units=None):
    """ Get Debye-Hückel radii for various ions

    For aqueous systems

    Parameters
    ----------
    key: str
        e.g. 'Fe/3+', 'SCN/-'
    units: object (optional)
        :attr:`nm` is accessed.

    Returns
    -------
    radius in nanometers

    """
    if units is None:
        return _radii_nm[key]
    else:
        return _radii_nm[key]*units.nm
