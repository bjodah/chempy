# -*- coding: utf-8 -*-

import warnings

try:
    from numpy import any as _any
except ImportError:
    def _any(arg):
        if arg is True:
            return True
        if arg is False:
            return False
        return any(arg)

# Parameters from paper (SI-units):

gamma = 2.063
dgamma = 0.051
D0 = 1.635e-8
TS = 215.05
dD0 = 2.242e-11
dTS = 1.2
low_t_bound = 273.15  # 0 deg C, (m.p. at ambient pressure)
high_t_bound = 373.15  # 100 deg C, (b.p. at ambient pressure)


def water_self_diffusion_coefficient(T=None, units=None, warn=True,
                                     err_mult=None):
    """
    Temperature-dependent self-diffusion coefficient of water.

    Parameters
    ----------
    T : float
        Temperature (default: in Kelvin)
    units : object (optional)
        object with attributes: Kelvin, meter, kilogram
    warn : bool (default: True)
        Emit UserWarning when outside temperature range.
    err_mult : length 2 array_like (default: None)
        Perturb paramaters D0 and TS with err_mult[0]*dD0 and
        err_mult[1]*dTS respectively, where dD0 and dTS are the
        reported uncertainties in the fitted paramters. Useful
        for estimating error in diffusion coefficient.

    References
    ----------
    Temperature-dependent self-diffusion coefficients of water and six selected
        molecular liquids for calibration in accurate 1H NMR PFG measurements
        Manfred Holz, Stefan R. Heila, Antonio Saccob;
        Phys. Chem. Chem. Phys., 2000,2, 4740-4742
        http://pubs.rsc.org/en/Content/ArticleLanding/2000/CP/b005319h
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
    if T is None:
        T = 298.15*K
    _D0 = D0 * m**2 * s**-1
    _TS = TS * K
    if err_mult is not None:
        _dD0 = dD0 * m**2 * s**-1
        _dTS = dTS * K
        _D0 += err_mult[0]*_dD0
        _TS += err_mult[1]*_dTS
    if warn and (_any(T < low_t_bound*K) or _any(T > high_t_bound*K)):
        warnings.warn("Temperature is outside range (0-100 degC)")
    return _D0*((T/_TS) - 1)**gamma

# bibtex format (generated at doi2bib.org):
reference = {
    'doi': '10.1039/b005319h',
    'url': 'http://dx.doi.org/10.1039/B005319H',
    'year ': 2000,
    'publisher': 'Royal Society of Chemistry ({RSC})',
    'volume': 2,
    'number': 20,
    'pages': (4740, 4742),
    'author': 'Manfred Holz and Stefan R. Heil and Antonio Sacco',
    'title': ('Temperature-dependent self-diffusion coefficients of water and'
              ' six selected molecular liquids for calibration in accurate 1H'
              ' {NMR} {PFG} measurements'),
    'journal': 'Phys. Chem. Chem. Phys.'
}
