# -*- coding: utf-8 -*-
"""
This module calculates the density of aqueous ethanol by interpolation
from tabulated data.
"""

from __future__ import (absolute_import, division, print_function)

import warnings

import pkg_resources
import numpy as np
from scipy.interpolate import interp1d

from .water_density_tanaka_2001 import water_density

_data = np.genfromtxt(pkg_resources.resource_filename(__name__, 'ethanol-water.csv'), delimiter=',')
_T_C = [None, None, None, 10, 20, 25, 30]  # 4degC water as ref., later cols have other ref.


def aqueous_ethanol_density(w, T=None, T0=None, units=None, warn=True):
    """
    Density of aqueous ethanol (kg/m³) as function of temperature (K)
    and mass fraction ethanol (w).

    Parameters
    ----------
    w: float
        Mass fraction (0.0 <= w <= 1.0).
    T: float
        Temperature (in Kelvin) (283.15 <= T <= 298.15) (default: 298.15)
    T0: float
        Value of T for 0 degree Celsius (default: 273.15)
    units: object (optional)
        object with attributes: kelvin, meter, kilogram
    warn: bool (default: True)
        Emit UserWarning when outside T or w range.

    Returns
    -------
    Density (float of kg/m³ if T is float and units is None)

    Examples
    --------
    >>> print('%.4f' % aqueous_ethanol_density(.4))
    0.9315

    References
    ----------
    Lange's Handbook of Chemistry, 10th ed.

    """
    if units is None:
        K = 1
    else:
        K = units.Kelvin
    if T is None:
        T = 298.15*K
    if T0 is None:
        T0 = 273.15*K
    t = T - T0
    if warn:
        if np.any(t < 10) or np.any(t > 25):
            warnings.warn("Temperature is outside range (10-25 degC)")
        if np.any(w < 0) or np.any(w > 1):
            warnings.warn("Mass fraction is outside range (0.0-1.0)")
    if t < 20:
        cols = (3, 4)
    elif t < 25:
        cols = (4, 5)
    else:
        cols = (5, 6)
    t_lo = _T_C[cols[0]]
    t_hi = _T_C[cols[1]]
    t_span = t_hi - t_lo
    x_hi = (t - t_lo)/t_span
    x_lo = 1 - x_hi
    f_lo = interp1d(_data[:, 0], _data[:, cols[0]])
    f_hi = interp1d(_data[:, 1], _data[:, cols[1]])
    y_lo = f_lo(w)
    y_hi = f_hi(w)
    return (x_lo*y_lo + x_hi*y_hi)*water_density(277.15*K, T0, units=units, warn=warn)
