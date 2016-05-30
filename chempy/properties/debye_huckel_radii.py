# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


# From Kielland, Individual Activity Coefficients of Ions in Aqueous Solutions, Vol. 59, 1937, 1675-1678
_radii_nm = {  # in nanometer, from Table I, Page 1676
    'H+': .9,
    'Li+': .6,
    'Na+': .425,
    'K+': .3,
    'Rb+': .25,
    'Cs+': .25,
    'NH4+': .25,
    'Tl+': .25,
    'Ag+': .25,
    'Be+2': .8,
    'Mg+2': .8,
    'Ca+2': .6,
    'Sr+2': .5,
    'Ba+2': .5,
    'Ra+2': .5,
    'Cu+2': .6,
    'Zn+2': .6,
    'Cd+2': .5,
    'Hg+2': .5,
    'Pb+2': .45,
    'Mn+2': .6,
    'Fe+2': .6,
    'Ni+2': .6,
    'Co+2': .6,
    'Sn+2': .6,
    'Al+3': .9,
    'Fe+3': .9,
    'Cr+3': .9,
    'La+3': .9,
    'Ce+3': .9,
    'Pr+3': .9,
    'Nd+3': .9,
    'Sc+3': .9,
    'Sm+3': .9,
    'Y+3': .9,
    'In+3': .9,
    'Sn+4': 1.1,
    'Th+4': 1.1,
    'Zr+4': 1.1,
    'Ce+4': 1.1,
    'F-': .35,
    'Cl-': .3,
    'Br-': .3,
    'I-': .3,
    'ClO3-': .35,
    'ClO4-': .35,
    'NO3-': .3,
    'BrO3-': .35,
    'IO3-': .425,
    'HCO3-': .425,
    'S-2': .5,
    'CO3-2': .45,
    'SO4-2': .4,
    'C2O4-2': .45,
    'SCN-': 0.35,
    'OH-': .35,
}


def get_radii(key, units=None):
    """ Get Debye-Hückel radii for various ions

    For aqueous systems

    Parameters
    ----------
    key: str
        e.g. 'Fe+3', 'SCN-'
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
