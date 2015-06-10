from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

try:
    from numpy import exp as _exp
except ImportError:
    from math import exp as _exp


class Henry(object):
    """
    Henry's gas constant. Note that the reference temperature
    is set by the attribute :py:attr:`T0` which defaults to
    298.15 (Kelvin).

    Parameters
    ----------
    kH0: float
        Henry's constant [M/atm]
    derivative: float
        dln(kH)/d(1/T) [K]
        Equivalent to \Delta_soln H / R
    ref: object
        Reference for origin of parameters
    units: object (optional)
        object with attributes: kelvin
    """

    def __init__(self, kH0, derivative, T0=None, units=None, ref=None):
        if units is None:
            kelvin = 1
        else:
            kelvin = units.kelvin
        self.kH0 = kH0
        self.derivative = derivative
        if T0 is None:
            self.T0 = 298.15 * kelvin
        else:
            self.T0 = T0
        self.ref = ref

    def get_kH_at_T(self, T, exp=None):
        """ Evaluate kH at temperature T """
        if exp is None:
            exp = _exp
        return self.kH0 * exp(
            self.derivative*(1/T - 1/self.T0))

    def get_c_at_T_and_P(self, T, P):
        """
        Calculate what concentration is needed for achieving a given partial
        pressure at a specified temperature
        """
        return P * self.get_kH_at_T(T)

    def get_P_at_T_and_c(self, T, c):
        """
        Calculate the partial pressure for given temperature and concentration
        """
        return c / self.get_kH_at_T(T)
