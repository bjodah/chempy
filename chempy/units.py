# -*- coding: utf-8 -*-
""" The units module provides the following attributes:

- ``chempy.units.default_units``
- ``chempy.units.default_constants``
- ``chempy.units.SI_base_registry``

Currently `quantities <https://pypi.python.org/pypi/quantities>`_ is used as
the underlying package to handle units. If it is possible you should try to
only use the `chempy.units` module (in case ``ChemPy`` changes this backend).

"""
from __future__ import (absolute_import, division, print_function)

# Currently we use quantities for units. This may change, therefore use this
# file for all units. A requirement is first-class numpy support.

import numpy as np
try:
    import quantities as pq
except ImportError:
    UncertainQuantity = None
    default_constants = None
    default_units = None
    SI_base_registry = None
else:
    UncertainQuantity = pq.UncertainQuantity
    # Let us extend the underlying pq namespace with some common units in
    # chemistry
    default_constants = pq.constants

    class _NameSpace:
        def __init__(self, default):
            self._NameSpace_default = default
            self._NameSpace_attr_store = {}

        def __getattr__(self, attr):
            if attr.startswith('_NameSpace_'):
                return self.__dict__[attr]
            else:
                try:
                    return self._NameSpace_attr_store[attr]
                except KeyError:
                    return getattr(self._NameSpace_default, attr)

        def __setattr__(self, attr, val):
            if attr.startswith('_NameSpace_'):
                self.__dict__[attr] = val
            else:
                self._NameSpace_attr_store[attr] = val

        def as_dict(self):
            result = {k: v for k, v in self._NameSpace_default.__dict__.items()
                      if not k.startswith('_')}
            result.update(self._NameSpace_attr_store)
            return result

    default_units = _NameSpace(pq)
    default_units.decimetre = pq.UnitQuantity(
        'decimetre',  default_units.m / 10.0, u_symbol='dm')
    if not hasattr(default_units, 'molar'):
        default_units.molar = pq.UnitQuantity(
            'M',  default_units.mole / default_units.decimetre ** 3,
            u_symbol='M')
    default_units.per100eV = pq.UnitQuantity(
        'per_100_eV',
        1/(100*default_units.eV*default_constants.Avogadro_constant),
        u_symbol='(100eV)**-1')
    default_units.micromole = pq.UnitQuantity(
        'micromole',  pq.mole/1e6,  u_symbol=u'μmol')
    default_units.kilojoule = pq.UnitQuantity(
        'kilojoule',  1e3*pq.joule,  u_symbol='kJ')
    default_units.perMolar_perSecond = 1/default_units.molar/pq.s
    default_units.per100eV = pq.UnitQuantity(
        'per_100_eV', 1/(100*pq.eV*pq.constants.Avogadro_constant),
        u_symbol='(100eV)**-1')
    default_units.umol = pq.UnitQuantity('micromole',  pq.mole/1e6,
                                         u_symbol=u'μmol')
    default_units.umol_per_J = default_units.umol / pq.joule

    # unit registry data and logic:

    SI_base_registry = {
        'length': default_units.metre,
        'mass': default_units.kilogram,
        'time': default_units.second,
        'current': default_units.ampere,
        'temperature': default_units.kelvin,
        'luminous_intensity': default_units.candela,
        'amount': default_units.mole
    }


def get_derived_unit(registry, key):
    """ Get the unit of a physcial quantity in a provided unit system.

    Parameters
    ----------
    registry: dict (str: unit)
        mapping 'length', 'mass', 'time', 'current', 'temperature',
        'luminous_intensity', 'amount'. If registry is ``None`` the
        function returns 1.0 unconditionally.
    key: str
        one of the registry keys or one of: 'diffusion', 'electrical_mobility',
        'permittivity', 'charge', 'energy', 'concentration', 'density',
        'radiolytic_yield'

    Examples
    --------
    >>> m, s = default_units.meter, default_units.second
    >>> get_derived_unit(SI_base_registry, 'diffusion') == m**2/s
    True

    """
    if registry is None:
        return 1.0
    derived = {
        'diffusion': registry['length']**2/registry['time'],
        'electrical_mobility': (registry['current']*registry['time']**2 /
                                registry['mass']),
        'permittivity': (registry['current']**2*registry['time']**4 /
                         (registry['length']**3*registry['mass'])),
        'charge': registry['current']*registry['time'],
        'energy': registry['mass']*registry['length']**2/registry['time']**2,
        'concentration': registry['amount']/registry['length']**3,
        'density': registry['mass']/registry['length']**3,
    }
    derived['radiolytic_yield'] = registry['amount']/derived['energy']
    try:
        return derived[key]
    except KeyError:
        return registry[key]


def unit_registry_to_human_readable(unit_registry):
    """ Serialization of a unit registry. """
    if unit_registry is None:
        return None
    new_registry = {}
    for k in SI_base_registry:
        if unit_registry[k] is 1:
            new_registry[k] = 1, 1
        else:
            dim_list = list(unit_registry[k].dimensionality)
            if len(dim_list) != 1:
                raise TypeError("Compound units not allowed: {}".format(
                    dim_list))
            u_symbol = dim_list[0].u_symbol
            # u_symbol = unit_registry[k].u_symbol
            new_registry[k] = float(unit_registry[k]), u_symbol
    return new_registry


def unit_registry_from_human_readable(unit_registry):
    """ Deserialization of unit_registry. """
    if unit_registry is None:
        return None
    new_registry = {}
    for k in SI_base_registry:
        factor, u_symbol = unit_registry[k]
        if u_symbol == 1:
            unit_quants = [1]
        else:
            unit_quants = list(pq.Quantity(0, u_symbol).dimensionality.keys())

        if len(unit_quants) != 1:
            raise TypeError("Unkown UnitQuantity: {}".format(unit_registry[k]))
        else:
            new_registry[k] = factor*unit_quants[0]
    return new_registry


# Abstraction of underlying package providing units and dimensional analysis:

def is_unitless(expr):
    """ Returns ``True`` if ``expr`` is unitless, otherwise ``False``

    Examples
    --------
    >>> is_unitless(42)
    True
    >>> is_unitless(42*default_units.kilogram)
    False

    """
    if hasattr(expr, 'dimensionality'):
        return expr.dimensionality == pq.dimensionless.dimensionality
    return True


def unit_of(expr):
    """ Returns the unit of a quantity

    Examples
    --------
    >>> unit_of(42*pq.second) == unit_of(12*pq.second)
    True
    >>> unit_of(42)
    1

    """
    try:
        return expr.units
    except AttributeError:
        return 1


def to_unitless(value, new_unit=None):
    """ Nondimensionalization of a quantity.

    Parameters
    ----------
    value: quantity
    new_unit: unit

    Examples
    --------
    >>> '%.1g' % to_unitless(1*default_units.metre, default_units.nm)
    '1e+09'

    """
    if new_unit is None:
        new_unit = pq.dimensionless
    if isinstance(value, (list, tuple)):
        return np.array([to_unitless(elem, new_unit) for elem in value])
    if isinstance(value, dict):
        new_value = value.copy()
        for k in new_value:
            new_value[k] = to_unitless(new_value[k], new_unit)
        return new_value
    try:
        result = (value*pq.dimensionless/new_unit).rescale(pq.dimensionless)
        if result.ndim == 0:
            return float(result)
        else:
            return np.asarray(result)
    except TypeError:
        return np.array([to_unitless(elem, new_unit) for elem in value])


# NumPy like functions for compatibility:

def allclose(a, b, rtol=1e-8, atol=None):
    """ Analogous to ``numpy.allclose``. """
    try:
        d = abs(a - b)
    except TypeError:
        if len(a) == len(b):
            return all(allclose(_a, _b, rtol, atol) for _a, _b in zip(a, b))
        raise
    lim = abs(a)*rtol
    if atol is not None:
        lim += atol
    try:
        n = len(d)
    except TypeError:
        n = 1

    if n == 1:
        return d < lim
    else:
        return np.all(_d < _lim for _d, _lim in zip(d, lim))


def linspace(start, stop, num=50):
    """ Analogous to ``numpy.linspace``. """
    # work around for quantities v0.10.1 and NumPy
    unit = unit_of(start)
    start_ = to_unitless(start, unit)
    stop_ = to_unitless(stop, unit)
    return np.linspace(start_, stop_, num)*unit
