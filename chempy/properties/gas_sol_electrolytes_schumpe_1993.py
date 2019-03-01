# -*- coding: utf-8 -*-
"""
This module implements an expression for estimating how gas solubility in water
is affected by dissolved non-reactive salts, the expression and parameters are from Schumpe (1993).
"""

import warnings

Tref_K = 298.2  # Kelvin

p_ion_rM = {  # "rM" === "per molar"
    'H+': 0.0,
    'Li+': 0.0691,
    'Na+': 0.1171,
    'K+': 0.0959,
    'Rb+': 0.0845,
    'Cs+': 0.0660,
    'NH4+': 0.0539,
    'Mg+2': 0.1765,
    'Ca+2': 0.1771,
    'Ba+2': 0.2021,
    'Fe+2': 0.1712,
    'Co+2': 0.1983,
    'Ni+2': 0.2039,
    'Cu+2': 0.1810,
    'Mn+2': 0.1620,
    'Zn+2': 0.1712,
    'Cd+2': 0.2201,
    'Al+3': 0.2253,
    'Fe+3': 0.0996,
    'Cr+3': 0.0595,
    'OH-': 0.0756,
    'F-': 0.1016,
    'Cl-': 0.0334,
    'Br-': 0.0137,
    'I-': 0.0020,
    'NO3-': 0.0050,
    'ClO4-': 0.0502,
    'IO4-': 0.1514,
    'HCO3-': 0.1372,
    'HSO3-': 0.0543,
    'H2PO4-': 0.1025,
    # '0-O-': 0.0765,
    # '0-OCH2COO-': 0.0119,
    'S2O3-2': 0.1109,
    'HPO4-2': 0.1789,
    'CO3-2': 0.1666,
    'SO3-2': 0.1537,
    'SO4-2': 0.1185,
    'PO4-3': 0.2117
}

p_gas_rM = {
    'O2': 0.0,
    'CO2': -0.0183,
    'N2O': -0.0110,
    'C2H2': -0.0174,
    'C2H4': 0.0014,
    'He': -0.036,
    'Ne': -0.020,
    'Ar': -0.009,
    'Kr': 0.003,
    'Xe': 0.005,
    'Rn': 0.015,
    'H2': -0.024,
    'N2': -0.008,
    'NO': 0.004,
    'C2H6': 0.011,
}


def lg_solubility_ratio(electrolytes, gas, units=None, warn=True):
    """ Returns the log10 value of the solubilty ratio

    Implements equation 16, p 156. from Schumpe (1993)

    Parameters
    ----------
    electrolytes : dict
        Mapping substance key (one in ``p_ion_rM``) to concentration.
    gas : str
        Substance key for the gas (one in ``p_gas_rM``).
    units : object (optional)
        object with attribute: molar
    warn : bool (default: True)
        Emit UserWarning when 'F-' among electrolytes.

    """
    if units is None:
        M = 1
    else:
        M = units.molar
    if warn and 'F-' in electrolytes:
        warnings.warn("In Schumpe 1993: data for fluoride uncertain.")
    return sum([(p_gas_rM[gas]/M+p_ion_rM[k]/M)*v for k, v in electrolytes.items()])


reference = dict(
    doi='10.1016/0009-2509(93)80291-W',
    url='http://www.sciencedirect.com/science/article/pii/000925099380291W',
    year=1993,
    publisher='Pergamon',
    author='Schumpe, Adrian',
    title='The estimation of gas solubilities in salt solutions',
    volume=48,
    issn='0009-2509',
    abstract=(
        'The effects of dissolved salts on the solubilities of gases were analysed based on '
        'a comprehensive set of literature data for the temperature of 298.2 K. A new '
        'empirical model was suggested which, at no increase in the number of adjustable '
        'parameters, described the data with a lower standard deviation than previously '
        'suggested models. The parameter values evaluated for the new model allow to '
        'estimate the effects of 20 cations and 19 anions on the solubilities of 15 gases.'),
    pages=(153, 158),
    number=1,
    journaltitle='Chemical Engineering Science',
    shortjournal='Chemical Engineering Science',
)
