# -*- coding: utf-8 -*-

import warnings

try:
    from numpy import any as _any
except ImportError:
    _any = any


def water_self_diffusion_coefficient(T=298.15, units=None, warn=True):
    """
    Temperature-dependent self-diffusion coefficient of water.

    Parameters
    ----------
    T: float
        Temperature (default: in Kelvin)
    units: object (optional)
        object with attributes: Kelvin, meter, kilogram
    warn: bool (default: True)
        Emit UserWarning when outside temperature range.

    References
    ----------
    Temperature-dependent self-diffusion coefficients of water and six selected
        molecular liquids for calibration in accurate 1H NMR PFG measurements
        Manfred Holz, Stefan R. Heila, Antonio Saccob;
        Phys. Chem. Chem. Phys., 2000,2, 4740-4742
        DOI: 10.1039/B005319H
    """
    if units is None:
        K = 1
        m = 1
        s = 1
    else:
        K = units.Kelvin
        m = units.meter
        s = units.second
    D0 = 1.635e-8 * m**2 * s**-1
    dD0 = 2.242e-11 * m**2 * s**-1
    TS = 215.05 * K
    dTS = 1.2 * K
    gamma = 2.063
    dgamma = 0.051
    if warn and (_any(T < 273.15*K) or _any(T > 373.15*K)):
        warnings.warn("Temperature is outside range (0-100 degC)")
    return D0*((T/TS) - 1)**gamma

# generated at doi2bib.org:
bibtex = """
@article{Holz2000,
  doi = {10.1039/b005319h},
  url = {http://dx.doi.org/10.1039/B005319H},
  year  = {2000},
  publisher = {Royal Society of Chemistry ({RSC})},
  volume = {2},
  number = {20},
  pages = {4740--4742},
  author = {Manfred Holz and Stefan R. Heil and Antonio Sacco},
  title = {Temperature-dependent self-diffusion coefficients of water and six
           selected molecular liquids for calibration in accurate 1H {NMR}
           {PFG} measurements},
  journal = {Phys. Chem. Chem. Phys.}
}
"""
