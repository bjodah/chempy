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

from .chemistry import Substance
from .units import default_units as u
from .util.arithmeticdict import ArithmeticDict


class QuantityDict(ArithmeticDict):
    unit = None

    def __init__(self, *args, **kwargs):
        super(QuantityDict, self).__init__(lambda: 0*self.unit, *args, **kwargs)

    def copy(self):
        return self.__class__(self.items())


class ConcDict(QuantityDict):
    unit = u.molar


class AmountDict(QuantityDict):
    unit = u.mol


class MassDict(QuantityDict):
    unit = u.gram


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
        self.volume = volume
        self.concentrations = ConcDict(concentrations)
        if substances is None:
            substances = AutoRegisteringSubstanceDict()
        self.substances = substances
        self.solvent = solvent

    def __add__(self, other):
        if self.solvent != other.solvent:
            raise NotImplementedError("Mixed solvent should be represented as concentrations")
        tot_amount = AmountDict(self.concentrations*self.volume) + AmountDict(other.concentrations*other.volume)
        tot_vol = self.volume + other.volume
        return Solution(tot_vol, ConcDict(tot_amount / tot_vol), self.substances, self.solvent)

    def dissolve(self, masses):
        contrib = ConcDict({k: v/self.substances[k].molar_mass()/self.volume for k, v in masses.items()})
        return Solution(self.volume, self.concentrations + contrib, self.substances, self.solvent)
