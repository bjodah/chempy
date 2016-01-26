# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


def electrical_mobility_from_D(D, charge, T, constants=None, units=None):
    """
    Calculates the electrical mobility through Einstein-Smoluchowski relation.

    Parameters
    ----------
    D: float with unit
        Diffusion coefficient
    charge: integer
        charge of the species
    T: float with unit
        Absolute temperature
    constants: object (optional, default: None)
        if None:
            T assumed to be in Kelvin and b0 = 1 mol/kg
        else:
            see source code for what attributes are used.
            Tip: pass quantities.constants
    units: object (optional, default: None)
        attributes accessed: meter, Kelvin and mol

    Returns
    -------
    Electrical mobility

    """

    if constants is None:
        kB = 1.38064852e-23
        e = 1.60217662e-19
        if units is not None:
            kB *= units.joule / units.kelvin / units.mol
            e *= units.coulomb
    else:
        kB = constants.Boltzmann_constant
        e = constants.elementary_charge
    return D*charge*e/(kB*T)
