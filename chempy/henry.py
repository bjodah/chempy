from __future__ import absolute_import, division, print_function


def Henry_H_at_T(T, H, Tderiv, T0=None, units=None, exp=None):
    """ Evaluate Henry's constant H at temperature T

    Parameters
    ----------
    T: float
        Temperature (with units), assumed to be in Kelvin if ``units == None``
    H: float
        Henry's constant
    Tderiv: float (optional)
        dln(H)/d(1/T), assumed to be in Kelvin if ``units == None``.
        default: 298.15 K
    T0: float
        Reference temperature, assumed to be in Kelvin if ``units == None``
    units: object (optional)
        object with attributes: kelvin (e.g. chempy.units.default_units)
    exp: callback (optional)
        callback for calculating the exponential, default: numpy.exp, math.exp

    """
    if exp is None:
        try:
            from numpy import exp
        except ImportError:
            from math import exp

    if units is None:
        K = 1
    else:
        K = units.Kelvin

    if T0 is None:
        T0 = 298.15*K

    return H * exp(Tderiv*(1/T - 1/T0))


class Henry(defaultnamedtuple('Henry', 'Hcp Tderiv T0 ref', [298.15, None])):
    """ Henry's gas constant

    Note that the reference temperature
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

    Examples
    --------
    >>> from chempy.units import to_unitless, default_units as u
    >>> H_CO = Henry(9.7e-6 * u.mol/u.m3/u.Pa, 1300*u.K, 'sander_2015')
    >>> '%.2g' % to_unitless(H_CO(298.15*u.K), u.molar/u.bar)
    1.2e-3

    """

    def __call__(self, T, units=None, exp=None):
        """ Evaluates Henry's constant for provided temperature """
        return Henry_H_at_T(self.Hcp, self.Tderiv, T, units=units)

    def get_c_at_T_and_P(self, T, P, **kwargs):
        """ Convenience method for calculating concentration

        Calculate what concentration is needed to achieve a given partial
        pressure at a specified temperature

        Parameters
        ----------
        T: float
            Temperature
        P: float
            Pressure
        \*\*kwargs:
            Keyword arguments passed on to :meth:`__call__`

        """
        return P * self(T, **kwargs)

    def get_P_at_T_and_c(self, T, c, **kwargs):
        """ Convenience method for calculating concentration

        Calculate the partial pressure for given temperature and concentration


        Parameters
        ----------
        T: float
            Temperature
        P: float
            Pressure
        \*\*kwargs:
            Keyword arguments passed on to :meth:`__call__`
        """
        return c / self(T, **kwargs)
