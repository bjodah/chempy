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
        temperature
    rho: float with unit
        density
    b0: float with unit
        reference molality
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


def B(eps_r, T, one=1, units=None):
    # TODO: proper doc string., find original reference
    param = 50.3
    if units is not None:
        try:
            nm = units.nanometer
        except AttributeError:
            nm = 1e-3 * units.meter
        param *= units.Kelvin**(one/2) / nm
    return param*(eps_r*T)**(-one/2)


def limiting_log_gamma(I, z, A, I0=1, one=1):
    # `one` allows passing of e.g. one=sympy.S(1)
    return -A*z**2*(I/I0)**(one/2)


def extended_log_gamma(I, z, a, A, B, C=0, I0=1, one=1):
    # `one` allows passing of e.g. one=sympy.S(1)
    return -A*z**2*(I/I0)**(one/2)/(1 + B*a*(I/I0)**(one/2)) + C*I/I0


def limiting_activity_product(I, stoich, z, T, eps_r, exp=None):
    if exp is None:
        try:
            from numpy import exp
        except ImportError:
            from math import exp
    Aval = A(eps_r, T)/ln(10)
    tot = 0
    for idx, nr in enumerate(stoich):
        tot += nr*limiting_log_gamma(I, z[idx], Aval)
    return exp(tot)


def extended_activity_product(I, stoich, z, a, T, eps_r, C=0, exp=None):
    if exp is None:
        try:
            from numpy import exp
        except ImportError:
            from math import exp
    if ln is None:
        from math import log as ln
    Aval = A(eps_r, T)
    Bval = B(eps_r, T)
    tot = 0
    for idx, nr in enumerate(stoich):
        tot += nr*extended_log_gamma(I, z[idx], a[idx], A, B, C)
    return exp(tot)


class LimitingDebyeHuckelActivityProduct(ActivityProduct):

    def __call__(self, c):
        I = ionic_strength(c, self.z)
        return limiting_activity_product(I, self.stoich, *self.args)


class ExtendedDebyeHuckelActivityProduct(ActivityProduct):

    def __call__(self, c):
        z = self.args[0]
        I = ionic_strength(c, z)
        return extended_activity_product(I, self.stoich, *self.args)
