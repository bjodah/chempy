# -*- coding: utf-8 -*-
"""
Module for dealing with constants `Henry's law
<https://en.wikipedia.org/wiki/Henry's_law>`_.
"""
from __future__ import absolute_import, division, print_function

from ._util import get_backend
from .util.pyutil import defaultnamedtuple, deprecated
from .units import default_units


def Henry_H_at_T(T, H, Tderiv, T0=None, units=None, backend=None):
    """ Evaluate Henry's constant H at temperature T

    Parameters
    ----------
    T: float
        Temperature (with units), assumed to be in Kelvin if ``units == None``
    H: float
        Henry's constant
    Tderiv: float (optional)
        dln(H)/d(1/T), assumed to be in Kelvin if ``units == None``.
    T0: float
        Reference temperature, assumed to be in Kelvin if ``units == None``
        default: 298.15 K
    units: object (optional)
        object with attributes: kelvin (e.g. chempy.units.default_units)
    backend : module (optional)
        module with "exp", default: numpy, math

    """
    be = get_backend(backend)
    if units is None:
        K = 1
    else:
        K = units.Kelvin
    if T0 is None:
        T0 = 298.15*K
    return H * be.exp(Tderiv*(1/T - 1/T0))


class Henry(defaultnamedtuple('Henry', 'Hcp Tderiv T0 ref', [None, None])):
    """ Henry's gas constant

    Note that the reference temperature
    is set by the attribute :py:attr:`T0` which defaults to
    298.15 (Kelvin).

    Parameters
    ----------
    Hcp: float
        Henry's constant [M/atm]
    Tderiv: float
        dln(kH)/d(1/T) [K]
        Equivalent to $\\Delta_{soln}H / R$
    ref: object
        Reference for origin of parameters
    units: object (optional)
        object with attributes: kelvin

    Examples
    --------
    >>> H_H2 = Henry(1.2e-3, 1800, ref='carpenter_1966')
    >>> '%.2g' % H_H2(298.15)
    '0.0012'

    """

    def __call__(self, T, units=None, backend=None):
        """ Evaluates Henry's constant for provided temperature """
        return Henry_H_at_T(T, self.Hcp, self.Tderiv, self.T0,
                            units=units, backend=backend)

    @deprecated('0.3.1', '0.5.0', __call__)
    def get_kH_at_T(self, *args, **kwargs):
        return self(*args, **kwargs)

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
        \\*\\*kwargs:
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
        \\*\\*kwargs:
            Keyword arguments passed on to :meth:`__call__`
        """
        return c / self(T, **kwargs)


class HenryWithUnits(Henry):
    """ Analogous to :class:`Henry`

    Examples
    --------
    >>> from chempy.units import to_unitless, default_units as u
    >>> H_CO = HenryWithUnits(9.7e-6 * u.mol/u.m**3/u.Pa, 1300*u.K, ref='sander_2015')
    >>> '%.2g' % to_unitless(H_CO(298.15 * u.K), u.molar/u.bar)
    '0.00097'

    """
    def __call__(self, T, units=default_units, backend=None):
        """ Evaluates Henry's constant for provided temperature """
        return super(HenryWithUnits, self).__call__(T, units, backend)
