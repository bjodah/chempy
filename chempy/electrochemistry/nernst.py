# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)
import math


def nernst_potential(ion_conc_out, ion_conc_in, charge, T,
                     constants=None, units=None, backend=math):
    """
    Calculates the Nernst potential using the Nernst equation for a particular
    ion.

    Parameters
    ----------
    ion_conc_out : float with unit
        Extracellular concentration of ion.
    ion_conc_in : float with unit
        Intracellular concentration of ion.
    charge : integer
        Charge of the ion.
    T : float with unit
        Absolute temperature.
    constants : object (optional, default: None)
        Constant attributes accessed:
            F - Faraday constant
            R - Ideal Gas constant
    units : object (optional, default: None)
        Unit attributes: coulomb, joule, kelvin, mol.
    backend : module (optional, default: math)
        Module used to calculate log using `log` method, can be substituted
        with sympy to get symbolic answers.

    Returns
    -------
    Membrane potential.

    """
    if constants is None:
        F = 96485.33289
        R = 8.3144598
        if units is not None:
            F *= units.coulomb / units.mol
            R *= units.joule / units.kelvin / units.mol
    else:
        F = constants.Faraday_constant
        R = constants.molar_gas_constant

    return (R * T) / (charge * F) * backend.log(ion_conc_out / ion_conc_in)
