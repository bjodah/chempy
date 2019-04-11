# -*- coding: utf-8 -*-
"""
This modules collects expressions related to ionic strenght, e.g. the Debye-Hückel expressions.
"""
from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict
import warnings

from ._util import get_backend
from .chemistry import Substance
from .units import allclose


def _get_b0(b0, units=None):
    if units is not None and b0 is 1:
        return b0*units.molal
    else:
        return b0


def ionic_strength(molalities, charges=None, units=None, substances=None,
                   substance_factory=Substance.from_formula, warn=True):
    """ Calculates the ionic strength

    Parameters
    ----------
    molalities: array_like or dict
        Optionally with unit (amount / mass).
        when dict: mapping substance key to molality.
    charges: array_like
        Charge of respective ion, taken for substances when None.
    units: object (optional, default: None)
        Attributes accessed: molal.
    substances: dict, optional
        Mapping of substance keys to Substance instances (used when molalities
        is a dict).
    substance_factory: callback
        Used if `substances` is a string.
    warn: bool
        Issue a warning if molalities violates net charge neutrality.

    Examples
    --------
    >>> ionic_strength([1e-3, 3e-3], [3, -1]) == .5 * (9 + 3) * 1e-3
    True
    >>> ionic_strength({'Mg+2': 6, 'PO4-3': 4})
    30.0

    """
    tot = None
    if charges is None:
        if substances is None:
            substances = ' '.join(molalities.keys())
        if isinstance(substances, str):
            substances = OrderedDict([(k, substance_factory(k)) for k
                                      in substances.split()])
        charges, molalities = zip(*[(substances[k].charge, v) for k, v in molalities.items()])
    if len(molalities) != len(charges):
        raise ValueError("molalities and charges of different lengths")
    for b, z in zip(molalities, charges):
        if tot is None:
            tot = b * z**2
        else:
            tot += b * z**2
    if warn:
        net = None
        for b, z in zip(molalities, charges):
            if net is None:
                net = b * z
            else:
                net += b * z
        if not allclose(net, tot*0, atol=tot*1e-14):
            warnings.warn("Molalities not charge neutral: %s" % str(net))
    return tot/2


class _ActivityProductBase(object):
    """ Baseclass for activity products """

    def __init__(self, stoich, *args):
        self.stoich = stoich
        self.args = args

    def __call__(self, c):
        pass


def A(eps_r, T, rho, b0=1, constants=None, units=None, backend=None):
    """
    Debye Huckel constant A

    Parameters
    ----------
    eps_r: float
        relative permittivity
    T: float with unit
        Temperature (default: assume Kelvin)
    rho: float with unit
        density (default: assume kg/m**3)
    b0: float, optional
        Reference molality, optionally with unit (amount / mass)
        IUPAC defines it as 1 mol/kg. (default: 1).
    units: object (optional, default: None)
        attributes accessed: meter, Kelvin and mol
    constants: object (optional, default: None)
        if None:
            T assumed to be in Kelvin and b0 = 1 mol/kg
        else:
            see source code for what attributes are used.
            Tip: pass quantities.constants

    Notes
    -----
    Remember to divide by ln(10) if you want to use the constant
    with log10 based expression.

    References
    ----------
    Atkins, De Paula, Physical Chemistry, 8th edition

    """
    b0 = _get_b0(b0, units)
    be = get_backend(backend)
    one = be.pi**0
    if constants is None:
        combined = 132871.85866393594
        if units is not None:
            m = units.meter
            K = units.Kelvin
            mol = units.mol
            combined *= (m*K)**(3*one/2) / mol**(one/2)
        return combined*(rho * b0 * T**-3 * eps_r**-3)**0.5
    F = constants.Faraday_constant
    NA = constants.Avogadro_constant
    eps0 = constants.vacuum_permittivity
    kB = constants.Boltzmann_constant
    pi = constants.pi
    A = F**3/(4*pi*NA)*(rho*b0/(2*(eps0*eps_r*kB*NA*T)**3))**(one/2)
    return A


def B(eps_r, T, rho, b0=1, constants=None, units=None, backend=None):
    """
    Extended Debye-Huckel parameter B

    Parameters
    ----------
    eps_r: float
        relative permittivity
    T: float with unit
        temperature
    rho: float with unit
        density
    b0: float with unit
        reference molality
    units: object (optional, default: None)
        attributes accessed: meter, Kelvin and mol
    constants: object (optional, default: None)
        if None:
            T assumed to be in Kelvin, rho in kg/m**3 and b0 = 1 mol/kg
        else:
            attributes accessed: molar_gas_constant, Faraday_constant
            Tip: pass quantities.constants

    Returns
    -------
    Debye Huckel B constant (default in m**-1)

    """
    b0 = _get_b0(b0, units)
    be = get_backend(backend)
    one = be.pi**0
    if constants is None:
        combined = 15903203868.740343
        if units is not None:
            m = units.meter
            K = units.Kelvin
            mol = units.mol
            combined *= (m*K/mol)**(one/2)
        return combined*(rho*b0/(T*eps_r))**0.5
    F = constants.Faraday_constant
    eps0 = constants.vacuum_permittivity
    R = constants.molar_gas_constant
    B = F*(2*rho*b0/(eps_r*eps0*R*T))**(one/2)
    return B


def limiting_log_gamma(I, z, A, I0=1, backend=None):
    """ Debye-Hyckel limiting formula """
    be = get_backend(backend)
    one = be.pi**0
    return -A*z**2*(I/I0)**(one/2)


def extended_log_gamma(I, z, a, A, B, C=0, I0=1, backend=None):
    """ Debye-Huckel extended formula """
    be = get_backend(backend)
    one = be.pi**0
    I_I0 = I/I0
    sqrt_I_I0 = (I_I0)**(one/2)
    return -A*z**2 * sqrt_I_I0/(1 + B*a*sqrt_I_I0) + C*I_I0


def davies_log_gamma(I, z, A, C=-0.3, I0=1, backend=None):
    """ Davies formula """
    be = get_backend(backend)
    one = be.pi**0
    I_I0 = I/I0
    sqrt_I_I0 = (I_I0)**(one/2)
    return -A * z**2 * (sqrt_I_I0/(1 + sqrt_I_I0) + C*I_I0)


def limiting_activity_product(I, stoich, z, T, eps_r, rho, backend=None):
    """ Product of activity coefficients based on DH limiting law. """
    be = get_backend(backend)
    Aval = A(eps_r, T, rho)
    tot = 0
    for idx, nr in enumerate(stoich):
        tot += nr*limiting_log_gamma(I, z[idx], Aval)
    return be.exp(tot)


def extended_activity_product(I, stoich, z, a, T, eps_r, rho, C=0, backend=None):
    be = get_backend(backend)
    Aval = A(eps_r, T, rho)
    Bval = B(eps_r, T, rho)
    tot = 0
    for idx, nr in enumerate(stoich):
        tot += nr*extended_log_gamma(I, z[idx], a[idx], Aval, Bval, C)
    return be.exp(tot)


def davies_activity_product(I, stoich, z, a, T, eps_r, rho, C=-0.3, backend=None):
    be = get_backend(backend)
    Aval = A(eps_r, T, rho)
    tot = 0
    for idx, nr in enumerate(stoich):
        tot += nr*davies_log_gamma(I, z[idx], Aval, C)
    return be.exp(tot)


class LimitingDebyeHuckelActivityProduct(_ActivityProductBase):

    def __call__(self, c):
        z = self.args[0]
        I = ionic_strength(c, z)
        return limiting_activity_product(I, self.stoich, *self.args)


class ExtendedDebyeHuckelActivityProduct(_ActivityProductBase):

    def __call__(self, c):
        z = self.args[0]
        I = ionic_strength(c, z)
        return extended_activity_product(I, self.stoich, *self.args)
