from __future__ import division  # Py2/3 compatibility

from .core import ActivityProduct, ionic_strength

_radii_nm = {  # in nanometer
    'Fe3+': 0.9,
    'SCN-': 0.35,
    'Fe2+': 0.6,
}


def A(eps_r, T, rho, b0=1, constants=None, units=None, one=1):
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
    b0: float with unit
        reference molality (default: assume mol/kg)
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


def B(eps_r, T, rho, b0=1, constants=None, units=None, one=1):
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


def limiting_log_gamma(I, z, A, I0=1, one=1):
    """ Debye-Hyckel limiting formula """
    # `one` allows passing of e.g. one=sympy.S(1)
    return -A*z**2*(I/I0)**(one/2)


def extended_log_gamma(I, z, a, A, B, C=0, I0=1, one=1):
    """ Debye-Huckel extended formula """
    # `one` allows passing of e.g. one=sympy.S(1)
    I_I0 = I/I0
    sqrt_I_I0 = (I_I0)**(one/2)
    return -A*z**2 * sqrt_I_I0/(1 + B*a*sqrt_I_I0) + C*I_I0


def davies_log_gamma(I, z, A, C=-0.3, I0=1, one=1):
    """ Davies formula """
    I_I0 = I/I0
    sqrt_I_I0 = (I_I0)**(one/2)
    return -A * z**2 * (sqrt_I_I0/(1 + sqrt_I_I0) + C*I_I0)


def limiting_activity_product(I, stoich, z, T, eps_r, rho, exp=None):
    if exp is None:
        try:
            from numpy import exp
        except ImportError:
            from math import exp
    Aval = A(eps_r, T, rho)
    tot = 0
    for idx, nr in enumerate(stoich):
        tot += nr*limiting_log_gamma(I, z[idx], Aval)
    return exp(tot)


def extended_activity_product(I, stoich, z, a, T, eps_r, rho, C=0, exp=None):
    if exp is None:
        try:
            from numpy import exp
        except ImportError:
            from math import exp
    Aval = A(eps_r, T, rho)
    Bval = B(eps_r, T, rho)
    tot = 0
    for idx, nr in enumerate(stoich):
        tot += nr*extended_log_gamma(I, z[idx], a[idx], Aval, Bval, C)
    return exp(tot)


def davies_activity_product(I, stoich, z, a, T, eps_r, rho, C=-0.3, exp=None):
    if exp is None:
        try:
            from numpy import exp
        except ImportError:
            from math import exp
    Aval = A(eps_r, T, rho)
    tot = 0
    for idx, nr in enumerate(stoich):
        tot += nr*davies_log_gamma(I, z[idx], Aval, C)
    return exp(tot)


class LimitingDebyeHuckelActivityProduct(ActivityProduct):

    def __call__(self, c):
        z = self.args[0]
        I = ionic_strength(c, z)
        return limiting_activity_product(I, self.stoich, *self.args)


class ExtendedDebyeHuckelActivityProduct(ActivityProduct):

    def __call__(self, c):
        z = self.args[0]
        I = ionic_strength(c, z)
        return extended_activity_product(I, self.stoich, *self.args)
