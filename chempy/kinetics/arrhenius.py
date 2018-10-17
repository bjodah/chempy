# -*- coding: utf-8 -*-
"""
Contains functions for the `Arrhenius equation
<https://en.wikipedia.org/wiki/Arrhenius_equation>`_
(:func:`arrhenius_equation`) and a convenience fitting routine
(:func:`fit_arrhenius_equation`).
"""
from __future__ import (absolute_import, division, print_function)

from .._util import get_backend
from ..util.regression import least_squares
from ..util.pyutil import defaultnamedtuple
from ..units import default_constants, default_units, format_string

try:
    import numpy as np
except ImportError:
    np = None


def _fit_linearized(backtransfm, lin_x, lin_y, lin_yerr):
    if len(lin_x) != len(lin_y):
        raise ValueError("k and T needs to be of equal length.")
    if lin_yerr is not None:
        if len(lin_yerr) != len(lin_y):
            raise ValueError("kerr and T needs to be of equal length.")
    lin_p, lin_vcv, lin_r2 = least_squares(lin_x, lin_y, lin_yerr)
    return [cb(lin_p) for cb in backtransfm]


def _fit(T, k, kerr, func, lin_x, lin_y, backtransfm, linearized=False):
    _lin_y = lin_y(T, k)
    if kerr is None:
        lin_yerr = 1
    else:
        lin_yerr = (abs(lin_y(k - kerr, T) - _lin_y) +
                    abs(lin_y(k + kerr, T) - _lin_y))/2

    lopt = _fit_linearized(backtransfm, lin_x(T, k), _lin_y, lin_yerr)
    if linearized:
        return lopt
    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(func, T, k, lopt, kerr)
    return popt, pcov


def _get_R(constants=None, units=None):
    if constants is None:
        R = 8.314472
        if units is not None:
            J = units.Joule
            K = units.Kelvin
            mol = units.mol
            R *= J/mol/K
    else:
        R = constants.molar_gas_constant.simplified
    return R


def arrhenius_equation(A, Ea, T, constants=None, units=None, backend=None):
    """
    Returns the rate coefficient according to the Arrhenius equation

    Parameters
    ----------
    A: float with unit
        frequency factor
    Ea: float with unit
        activation energy
    T: float with unit
        temperature
    constants: object (optional, default: None)
        if None:
            T assumed to be in Kelvin, Ea in J/(K mol)
        else:
            attributes accessed: molar_gas_constant
            Tip: pass quantities.constants
    units: object (optional, default: None)
        attributes accessed: Joule, Kelvin and mol
    backend: module (optional)
        module with "exp", default: numpy, math

    """
    be = get_backend(backend)
    R = _get_R(constants, units)
    try:
        RT = (R*T).rescale(Ea.dimensionality)
    except AttributeError:
        RT = R*T
    return A*be.exp(-Ea/RT)


def fit_arrhenius_equation(T, k, kerr=None, linearized=False, constants=None, units=None):
    """ Curve fitting of the Arrhenius equation to data points

    Parameters
    ----------
    T : float
    k : array_like
    kerr : array_like (optional)
    linearized : bool

    """
    return _fit(T, k, kerr, arrhenius_equation, lambda T, k: 1/T, lambda T, k: np.log(k),
                [lambda p: np.exp(p[0]), lambda p: -p[1]*_get_R(constants, units)], linearized=linearized)


def _fit_arrhenius_equation(T, k, kerr=None, linearized=False):
    """ Curve fitting of the Arrhenius equation to data points

    Parameters
    ----------
    k : array_like
    T : float
    kerr : array_like (optional)
    linearized : bool

    """
    if len(k) != len(T):
        raise ValueError("k and T needs to be of equal length.")
    from math import exp
    import numpy as np
    p = np.polyfit(1/T, np.log(k), 1)
    R = _get_R(constants=None, units=None)
    Ea = -R*p[0]
    A = exp(p[1])
    if linearized:
        return A, Ea
    from scipy.optimize import curve_fit
    if kerr is None:
        weights = None
    else:
        weights = 1/kerr**2
    popt, pcov = curve_fit(arrhenius_equation, T, k, [A, Ea], weights)
    return popt, pcov


class ArrheniusParam(defaultnamedtuple('ArrheniusParam', 'A Ea ref', [None])):
    """ Kinetic data in the form of an Arrhenius parameterisation

    Parameters
    ----------
    Ea: float
        activation energy
    A: float
        preexponential prefactor (Arrhenius type eq.)
    ref: object (default: None)
        arbitrary reference (e.g. string representing citation key)

    Examples
    --------
    >>> k = ArrheniusParam(1e13, 40e3)
    >>> '%.5g' % k(298.15)
    '9.8245e+05'

    """

    @classmethod
    def from_fit_of_data(cls, T, k, kerr=None, **kwargs):
        args, vcv = fit_arrhenius_equation(T, k, kerr)
        return cls(*args, **kwargs)

    def __call__(self, T, constants=None, units=None, backend=None):
        """ Evaluates the arrhenius equation for a specified state

        Parameters
        ----------
        T: float
        constants: module (optional)
        units: module (optional)
        backend: module (default: math)

        See also
        --------
        See :func:`chempy.arrhenius.arrhenius_equation`.

        """
        return arrhenius_equation(self.A, self.Ea, T, constants=constants,
                                  units=units, backend=backend)

    def Ea_over_R(self, constants, units, backend=None):
        return self.Ea/_get_R(constants, units)

    def as_RateExpr(self, unique_keys=None, constants=None, units=None, backend=None):
        from .rates import Arrhenius, MassAction
        args = [self.A, self.Ea_over_R(constants, units)]
        return MassAction(Arrhenius(args, unique_keys))

    def format(self, precision, tex=False):
        try:
            str_A, str_A_unit = format_string(self.A, precision, tex)
            str_Ea, str_Ea_unit = format_string(self.Ea, precision, tex)
        except Exception:
            str_A, str_A_unit = precision.format(self.A), '-'
            str_Ea, str_Ea_unit = precision.format(self.Ea), '-'
        return (str_A, str_A_unit), (str_Ea, str_Ea_unit)

    def equation_as_string(self, precision, tex=False):
        (str_A, str_A_unit), (str_Ea, str_Ea_unit) = self.format(precision, tex)
        if tex:
            return r"{}\exp \left(-\frac{{{}}}{{RT}} \right)".format(
                str_A, str_Ea + ' ' + str_Ea_unit), str_A_unit
        else:
            return "{}*exp(-{}/(R*T))".format(str_A, str_Ea + ' ' + str_Ea_unit), str_A_unit

    def __str__(self):
        return ' '.join(self.equation_as_string('%.5g'))


class ArrheniusParamWithUnits(ArrheniusParam):
    def __call__(self, state, constants=default_constants, units=default_units, backend=None):
        """ See :func:`chempy.arrhenius.arrhenius_equation`. """
        return super(ArrheniusParamWithUnits, self).__call__(state, constants, units, backend)

    def as_RateExpr(self, unique_keys=None, constants=default_constants, units=default_units):
        return super(ArrheniusParamWithUnits, self).as_RateExpr(unique_keys, constants, units)
