# -*- coding: utf-8 -*-
"""
"Sandbox" module for exploring API useful for digital labbooks.

Examples
--------
>>> from chempy.units import to_unitless, default_units as u
>>> s1 = Solution(0.1*u.dm3, {'CH3OH': 0.1 * u.molar})
>>> s2 = Solution(0.3*u.dm3, {'CH3OH': 0.4 * u.molar, 'Na+': 2e-3*u.molar, 'Cl-': 2e-3*u.molar})
>>> s3 = s1 + s2
>>> abs(to_unitless(s3.volume - 4e-4 * u.m**3, u.dm3)) < 1e-15
True
>>> s3.concentrations.isclose({'CH3OH': 0.325*u.molar, 'Na+': 1.5e-3*u.molar, 'Cl-': 1.5e-3*u.molar})
True
>>> s4 = s3.dissolve({'CH3OH': 1*u.gram})
>>> abs(s4.concentrations['CH3OH'] - (0.325 + 1/(12.011 + 4*1.008 + 15.999)/.4)*u.molar) < 1e-4
True

"""
from __future__ import (absolute_import, division, print_function)

import copy

from .chemistry import Substance
from .units import (
    get_derived_unit, html_of_unit, is_unitless, SI_base_registry,
    to_unitless, rescale, default_units as u
)
from .util.arithmeticdict import ArithmeticDict, _imul, _itruediv
from .printing import as_per_substance_html_table


class QuantityDict(ArithmeticDict):
    def __init__(self, units, *args, **kwargs):
        self.units = units
        super(QuantityDict, self).__init__(lambda: 0*self.units, *args, **kwargs)
        self._check()

    @classmethod
    def of_quantity(cls, quantity_name, *args, **kwargs):
        instance = cls(get_derived_unit(SI_base_registry, quantity_name), *args, **kwargs)
        instance.quantity_name = quantity_name
        return instance

    def rescale(self, new_units):
        return self.__class__(new_units, {k: rescale(v, new_units) for k, v in self.items()})

    def _repr_html_(self):
        if hasattr(self, 'quantity_name'):
            header = self.quantity_name.capitalize() + ' / '
        else:
            header = ''
        header += html_of_unit(self.units)
        tab = as_per_substance_html_table(to_unitless(self, self.units), header=header)
        return tab._repr_html_()

    def _check(self):
        for k, v in self.items():
            if not is_unitless(v/self.units):
                raise ValueError("entry for %s (%s) is not compatible with %s" % (k, v, self.units))

    def __setitem__(self, key, value):
        if not is_unitless(value/self.units):
            raise ValueError("entry for %s (%s) is not compatible with %s" % (key, value, self.units))
        super(QuantityDict, self).__setitem__(key, value)

    def copy(self):
        return self.__class__(self.units, copy.deepcopy(list(self.items())))

    def __repr__(self):
        return "{}({}, {})".format(self.__class__.__name__,
                                   repr(self.units),
                                   dict(self))

    def __mul__(self, other):
        d = dict(copy.deepcopy(list(self.items())))
        _imul(d, other)
        return self.__class__(self.units * getattr(other, 'units', 1), d)

    def __truediv__(self, other):
        d = dict(copy.deepcopy(list(self.items())))
        _itruediv(d, other)
        return self.__class__(self.units / getattr(other, 'units', 1), d)

    def __floordiv__(self, other):
        a = self.copy()
        if getattr(other, 'units', 1) != 1:
            raise ValueError("Floor division with quantities not defined")
        a //= other
        return a

    def __rtruediv__(self, other):
        """ other / self """
        return self.__class__(getattr(other, 'units', 1)/self.units,
                              {k: other/v for k, v in self.items()})

    def __rfloordiv__(self, other):
        """ other // self """
        return self.__class__(getattr(other, 'units', 1)/self.units,
                              {k: other//v for k, v in self.items()})


class AutoRegisteringSubstanceDict(object):

    def __init__(self, factory=Substance.from_formula):
        self.factory = factory
        self._store = {}

    def __getitem__(self, key):
        if key not in self._store:
            self._store[key] = self.factory(key)
        return self._store[key]


class Solution(object):

    def __init__(self, volume, concentrations, substances=None, solvent=None):
        if not is_unitless(volume/u.dm3):
            raise ValueError("volume need to have a unit (e.g. dm3)")
        self.volume = volume
        self.concentrations = QuantityDict(u.molar, concentrations)
        if substances is None:
            substances = AutoRegisteringSubstanceDict()
        self.substances = substances
        self.solvent = solvent

    def __eq__(self, other):
        if not isinstance(other, Solution):
            return NotImplemented
        return all([getattr(self, k) == getattr(other, k) for k in 'volume concentrations substances solvent'.split()])

    def __add__(self, other):
        if self.solvent != other.solvent:
            raise NotImplementedError("Mixed solvent should be represented as concentrations")
        tot_amount = self.concentrations*self.volume + other.concentrations*other.volume
        tot_vol = self.volume + other.volume
        return Solution(tot_vol, tot_amount / tot_vol, self.substances, self.solvent)

    def dissolve(self, masses):
        contrib = QuantityDict(u.molar, {k: v/self.substances[k].molar_mass()/self.volume for k, v in masses.items()})
        return Solution(self.volume, self.concentrations + contrib, self.substances, self.solvent)

    def withdraw(self, volume):
        if volume > self.volume:
            raise ValueError("Cannot withdraw a volume greater than the solution volume")
        if volume < volume*0:
            raise ValueError("Cannot withdraw a negative volume")
        self.volume -= volume
        return Solution(volume, self.concentrations, self.substances, self.solvent)
