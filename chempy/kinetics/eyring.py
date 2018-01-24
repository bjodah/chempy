# -*- coding: utf-8 -*-
"""
Contains functions for the `Eyring equation
<https://en.wikipedia.org/wiki/Eyring_equation>`_.
"""
from __future__ import (absolute_import, division, print_function)

import math

from .._util import get_backend
from ..util.regression import least_squares
from ..util.pyutil import defaultnamedtuple
from ..units import default_units, Backend, default_constants, format_string
from .arrhenius import _get_R, _fit

try:
    import numpy as np
except ImportError:
    np = None


def _get_kB_over_h(constants=None, units=None):
    if constants is None:
        kB_over_h = 2.083664399411865234375e10
        if units is not None:
            s = units.second
            K = units.kelvin
            kB_over_h /= s*K
    else:
        kB_over_h = (constants.Boltzmann_constant /
                     constants.Planck_constant)
    return kB_over_h


def eyring_equation(dH, dS, T, constants=None, units=None, backend=None):
    """
    Returns the rate coefficient according to the Eyring equation

    Parameters
    ----------
    dH: float with unit
        Enthalpy of activation.
    dS: float with unit
        Entropy of activation.
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
    kB_over_h = _get_kB_over_h(constants, units)

    try:
        RT = (R*T).rescale(dH.dimensionality)
    except AttributeError:
        RT = R*T

    try:
        kB_over_h = kB_over_h.simplified
    except AttributeError:
        pass

    return kB_over_h*T*be.exp(dS/R)*be.exp(-dH/RT)


def fit_eyring_equation(T, k, kerr=None, linearized=False, constants=None, units=None):
    """ Curve fitting of the Eyring equation to data points

    Parameters
    ----------
    T : float
    k : array_like
    kerr : array_like (optional)
    linearized : bool

    """
    R = _get_R(constants, units)
    ln_kb_over_h = math.log(_get_kB_over_h(constants, units))
    return _fit(T, k, kerr, eyring_equation, lambda T, k: 1/T, lambda T, k: np.log(k/T),
                [lambda p: -p[1]*R, lambda p: R*(p[0] - ln_kb_over_h)], linearized=linearized)


class EyringParam(defaultnamedtuple('EyringParam', 'dH dS ref', [None])):
    """ Kinetic data in the form of an Eyring parameterisation

    Parameters
    ----------
    dH : float
        Enthalpy of activation.
    dS : float
        Entropy of activation.
    ref: object (default: None)
        arbitrary reference (e.g. citation key or dict with bibtex entries)

    Examples
    --------
    >>> k = EyringParam(72e3, 61.4)
    >>> '%.5g' % k(298.15)
    '2435.4'

    """

    def __call__(self, T, constants=None, units=None, backend=None):
        """ Evaluates the Eyring equation for a specified state

        Parameters
        ----------
        T : float
        constants : module (optional)
        units : module (optional)
        backend : module (default: math)

        See also
        --------
        See :func:`chempy.eyring.eyring_equation`.

        """
        return eyring_equation(self.dH, self.dS, T, constants=constants,
                               units=units, backend=backend)

    def kB_h_times_exp_dS_R(self, constants=None, units=None, backend=math):
        R = _get_R(constants, units)
        kB_over_h = _get_kB_over_h(constants, units)
        return kB_over_h * backend.exp(self.dS/R)

    def dH_over_R(self, constants=None, units=None, backend=None):
        R = _get_R(constants, units)
        return self.dH/R

    def as_RateExpr(self, unique_keys=None, constants=None, units=None, backend=math):
        from .rates import Eyring, MassAction
        args = [self.kB_h_times_exp_dS_R(constants, units, backend),
                self.dH_over_R(constants, units)]
        return MassAction(Eyring(args, unique_keys))

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
            return (
                r"\frac{{k_B T}}{{h}}\exp \left(\frac{{{}}}{{R}} \right)"
                r" \exp \left(-\frac{{{}}}{{RT}} \right)").format(
                    str_A, str_Ea + ' ' + str_Ea_unit), str_A_unit
        else:
            return "kB*T/h*exp({}/R)*exp(-{}/(R*T))".format(
                str_A, str_Ea + ' ' + str_Ea_unit), str_A_unit

    def __str__(self):
        return self.equation_as_string('%.5g')


class EyringParamWithUnits(EyringParam):
    def __call__(self, state, constants=default_constants, units=default_units,
                 backend=None):
        """ See :func:`chempy.eyring.eyring_equation`. """
        if backend is None:
            backend = Backend()
        return super(EyringParamWithUnits, self).__call__(
            state, constants, units, backend)

    def as_RateExpr(self, unique_keys=None, constants=default_constants,
                    units=default_units, backend=None):
        if backend is None:
            backend = Backend()
        return super(EyringParamWithUnits, self).as_RateExpr(
            unique_keys, constants, units, backend)
