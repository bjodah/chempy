# -*- coding: utf-8 -*-
"""
This module implements the denisty parameterisation of aqueous sulfuric acid
of Myhre et al. from 1998.
"""
from __future__ import (absolute_import, division, print_function)

import numpy as np
import warnings

from ..util import NoConvergence


_data = np.array([[999.8426, 0.03345402, -0.005691304, 0, 0],
                  [547.2659, -5.300445, 0.01187671, 0.0005990008, 0],
                  [5262.95, 37.20445, 0.1201909, -0.004148594, 1.197973E-5],
                  [-62139.58, -287.767, -0.4064638, 0.01119488, 3.607768E-5],
                  [409029.3, 1270.854, 0.326971, -0.01377435, -2.633585E-5],
                  [-1596989, -3062.836, 0.1366499, 0.006373031, 0],
                  [3857411, 4083.714, -0.1927785, 0, 0],
                  [-5808064, -2844.401, 0, 0, 0],
                  [5301976, 809.1053, 0, 0, 0],
                  [-2682616, 0, 0, 0, 0],
                  [576428.8, 0, 0, 0, 0]])


def sulfuric_acid_density(w, T=None, T0=None, units=None, warn=True):
    """
    Density of sulfuric acid (kg/m³) as function of temperature (K)
    and mass fraction acid (w).

    Parameters
    ----------
    w: float
        Acid mass fraction (0.1 <= w <= 0.9)
    T: float
        Temperature (in Kelvin) (273 <= T <= 323) (default: 298.15)
    T0: float
        Value of T for 0 degree Celsius (default: 273.15)
    units: object (optional)
        object with attributes: kelvin, meter, kilogram
    warn: bool (default: True)
        Emit UserWarning when outside T or w range.

    Returns
    -------
    Density of sulfuric acid (float of kg/m³ if T is float and units is None)

    Examples
    --------
    >>> print('%d' % sulfuric_acid_density(.5, 293))
    1396

    References
    ----------
    Cathrine E. L. Myhre , Claus J. Nielsen ,* and Ole W. Saastad
        "Density and Surface Tension of Aqueous H2SO4 at Low Temperature"
        J. Chem. Eng. Data, 1998, 43 (4), pp 617–622
        http://pubs.acs.org/doi/abs/10.1021/je980013g
        DOI: 10.1021/je980013g
    """
    if units is None:
        K = 1
        m = 1
        kg = 1
    else:
        K = units.Kelvin
        m = units.meter
        kg = units.kilogram
    if T is None:
        T = 298.15*K
    m3 = m**3
    if T0 is None:
        T0 = 273.15*K
    t = T - T0
    if warn:
        if np.any(t < 0*K) or np.any(t > 50*K):
            warnings.warn("Temperature is outside range (0-50 degC)")
        if np.any(w < 0.1) or np.any(w > .9):
            warnings.warn("Mass fraction is outside range (0.1-0.9)")
    t_arr = np.array([float(t/K)**j for j in range(5)]).reshape((1, 5))
    w_arr = np.array([w**i for i in range(11)]).reshape((11, 1))
    return np.sum((t_arr*w_arr)*_data)*kg/m3  # Equation (2) in reference


# bibtex format (generated at doi2bib.org):
reference = {
    'doi': '10.1021/je980013g',
    'url': 'http://dx.doi.org/10.1021/je980013g',
    'year ': 1998,
    'month': 'jul',
    'publisher': 'American Chemical Society ({ACS})',
    'volume': 43,
    'number': 4,
    'pages': (617, 622),
    'author': 'Cathrine E. L. Myhre and Claus J. Nielsen and Ole W. Saastad',
    'title': 'Density and Surface Tension of Aqueous H2{SO}4 at Low Temperature',
    'journal': r'Journal of Chemical {\&} Engineering Data'
}


def density_from_concentration(conc, T=None, molar_mass=None,
                               rho_cb=sulfuric_acid_density,
                               units=None, atol=None, maxiter=10, warn=False, **kwargs):
    """ Calculates the density of a solution from its concentration

    Given a function which calculates the density of a solution from the mass
    fraction of the solute, this function calculates (iteratively) the density
    of said solution for a given concentration.

    Parameters
    ----------
    conc : float (optionally with units)
        Concentration (mol / m³).
    T : float (optionally with units)
        Passed to ``rho_cb``.
    molar_mass : float (optionally with units)
        Molar mass of solute.
    rho_cb : callback
        Callback with signature f(w, T, units=None) -> rho
        (default: :func:`sulfuric_acid_density`).
    units : object (optional)
        Object with attributes: meter, kilogram, mol.
    atol : float (optionally with units)
        Convergence criterion for fixed-point iteration
        (default: 1e-3 kg/m³).
    maxiter : int
        Maximum number of iterations (when exceeded a NoConvergence excpetion
        is raised).
    \\*\\*kwargs:
        Keyword arguments passed onto ``rho_cb``.

    Returns
    -------
    Density of sulfuric acid (float of kg/m³ if T is float and units is None)

    Examples
    --------
    >>> print('%d' % density_from_concentration(400, 293))
    1021

    Raises
    ------
    chempy.util.NoConvergence:
        When maxiter is exceeded

    """
    if units is None:
        m = 1
        kg = 1
        mol = 1
    else:
        m = units.meter
        kg = units.kilogram
        mol = units.mol
    kg_per_m3 = kg * m**-3
    if atol is None:
        atol = 1e-3 * kg_per_m3
    if molar_mass is None:
        molar_mass = (1.00794*2 + 32.066 + 4*15.9994)*1e-3 * kg / mol

    if units is not None:
        conc = conc.rescale(mol/m**3)
        molar_mass = molar_mass.rescale(kg/mol)

    rho = 1100 * kg_per_m3
    delta_rho = float('inf') * kg_per_m3

    iter_idx = 0
    while atol < abs(delta_rho):
        # fixed point iteration
        new_rho = rho_cb(conc*molar_mass/rho, T, units=units, warn=warn, **kwargs)
        delta_rho = new_rho - rho
        rho = new_rho
        iter_idx += 1
        if iter_idx > maxiter:
            raise NoConvergence("maxiter exceeded")
    return rho
