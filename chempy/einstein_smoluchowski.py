# -*- coding: utf-8 -*-


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
        If None, use numerical SI constants, with units attached when provided.
        Constant attributes accessed:
            Boltzmann_constant
            elementary_charge
        Tip: pass quantities.constants
    units: object (optional, default: None)
        attributes accessed: joule, kelvin and coulomb

    Returns
    -------
    Electrical mobility

    """

    if constants is None:
        kB = 1.38064852e-23
        e = 1.60217662e-19
        if units is not None:
            kB *= units.joule / units.kelvin
            e *= units.coulomb
    else:
        kB = constants.Boltzmann_constant
        e = constants.elementary_charge
    return D * charge * e / (kB * T)
