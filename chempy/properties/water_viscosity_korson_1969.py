import warnings
from .._util import _any


A, B, C = 1.1709, 0.001827, 89.93
eta20_cP = 1.0020


def water_viscosity(T=None, eta20=None, units=None, warn=True):
    """ Viscosity of water (cP) as function of temperature (K)

    Parameters
    ----------
    T : float
        Temperature (in Kelvin) (default: 298.15 K)
    eta20 : float
        Viscosity of water at 20 degree Celsius.
    units : object (optional)
        object with attributes: kelvin & centipoise
    warn : bool
        Emit UserWarning when outside temperature range.

    Returns
    -------
    Water viscosity at temperature ``T``.

    """
    if units is None:
        cP = 1
        K = 1
    else:
        cP = units.centipoise
        K = units.kelvin
    if T is None:
        T = 298.15*K
    if eta20 is None:
        eta20 = eta20_cP*cP
    t = T - 273.15*K
    if warn and (_any(t < 0*K) or _any(t > 100*K)):
        warnings.warn("Temperature is outside range (0-100 degC)")
    # equation (5) in the paper says "log" but they seem to mean "log10"
    # when comparing with Table II.
    return eta20 * 10**((A*(20 - t) - B*(t - 20)**2)/(t + C))


reference = dict(
    doi='10.1021/j100721a006',
    url='https://doi.org/10.1021/j100721a006',
    year=1969,
    month='jan',
    publisher='American Chemical Society ({ACS})',
    volume=73,
    number=1,
    pages=(34, 39),
    author='Lawrence Korson and Walter Drost-Hansen and Frank J. Millero',
    title='Viscosity of water at various temperatures',
    journal='The Journal of Physical Chemistry'
)
