# -*- coding: utf-8 -*-


# From Kielland, Individual Activity Coefficients of Ions in Aqueous Solutions, Vol. 59, 1937, 1675-1678
_radii_nm = {  # in nanometer, from Table I, Page 1676
    "H+": 0.9,
    "Li+": 0.6,
    "Na+": 0.425,
    "K+": 0.3,
    "Rb+": 0.25,
    "Cs+": 0.25,
    "NH4+": 0.25,
    "Tl+": 0.25,
    "Ag+": 0.25,
    "Be+2": 0.8,
    "Mg+2": 0.8,
    "Ca+2": 0.6,
    "Sr+2": 0.5,
    "Ba+2": 0.5,
    "Ra+2": 0.5,
    "Cu+2": 0.6,
    "Zn+2": 0.6,
    "Cd+2": 0.5,
    "Hg+2": 0.5,
    "Pb+2": 0.45,
    "Mn+2": 0.6,
    "Fe+2": 0.6,
    "Ni+2": 0.6,
    "Co+2": 0.6,
    "Sn+2": 0.6,
    "Al+3": 0.9,
    "Fe+3": 0.9,
    "Cr+3": 0.9,
    "La+3": 0.9,
    "Ce+3": 0.9,
    "Pr+3": 0.9,
    "Nd+3": 0.9,
    "Sc+3": 0.9,
    "Sm+3": 0.9,
    "Y+3": 0.9,
    "In+3": 0.9,
    "Sn+4": 1.1,
    "Th+4": 1.1,
    "Zr+4": 1.1,
    "Ce+4": 1.1,
    "F-": 0.35,
    "Cl-": 0.3,
    "Br-": 0.3,
    "I-": 0.3,
    "ClO3-": 0.35,
    "ClO4-": 0.35,
    "NO3-": 0.3,
    "BrO3-": 0.35,
    "IO3-": 0.425,
    "HCO3-": 0.425,
    "S-2": 0.5,
    "CO3-2": 0.45,
    "SO4-2": 0.4,
    "C2O4-2": 0.45,
    "SCN-": 0.35,
    "OH-": 0.35,
}


def get_radii(key, units=None):
    """Get Debye-Hückel radii for various ions

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
        return _radii_nm[key] * units.nm
