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


def water_density(T=None, T0=None, units=None, a=None,
                  just_return_a=False, warn=True):
    """
    Density of water (kg/m3) as function of temperature (K)
    according to VSMOW model between 0 and 40 degree Celsius.
    Fitted using Thiesen's equation.

    Parameters
    ----------
    T : float
        Temperature (in Kelvin) (default: 298.15).
    T0 : float
        Value of T for 0 degree Celsius (default: 273.15).
    units : object (optional)
        Object with attributes: Kelvin, meter, kilogram.
    a : array_like (optional)
        5 parameters to the equation.
    just_return_a : bool (optional, default: False)
        Do not compute rho, just return the parameters ``a``.
    warn : bool (default: True)
        Emit UserWarning when outside temperature range.

    Returns
    -------
    Density of water (float of kg/m3 if T is float and units is None)

    Examples
    --------
    >>> print('%.2f' % water_density(277.13))
    999.97

    References
    ----------
    TANAKA M., GIRARD G., DAVIS R., PEUTO A. and BIGNELL N.,
        "Recommanded table for the density of water between 0 °C and 40 °C
        based on recent experimental reports",
        Metrologia, 2001, 38, 301-309.
        http://iopscience.iop.org/article/10.1088/0026-1394/38/4/3
        doi:10.1088/0026-1394/38/4/3
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
    if a is None:
        a = (-3.983035*K,  # C
             301.797*K,  # C
             522528.9*K*K,  # C**2
             69.34881*K,  # C
             999.974950*kg/m3)
    if just_return_a:
        return a
    if T0 is None:
        T0 = 273.15*K
    t = T - T0
    if warn and (_any(t < 0*K) or _any(t > 40*K)):
        warnings.warn("Temperature is outside range (0-40 degC)")
    return a[4]*(1-((t + a[0])**2*(t + a[1]))/(a[2]*(t + a[3])))

# bibtex format (generated at doi2bib.org):
reference = {
    'doi': '10.1088/0026-1394/38/4/3',
    'url': 'http://dx.doi.org/10.1088/0026-1394/38/4/3',
    'year ': 2001,
    'month': 'aug',
    'publisher': '{IOP} Publishing',
    'volume': 38,
    'number': 4,
    'pages': (301, 309),
    'author': 'M Tanaka and G Girard and R Davis and A Peuto and N Bignell',
    'title': 'Recommended table for the density of water between 0 ~C and 40 ~C based on recent experimental reports',
    'journal': 'Metrologia'
}
