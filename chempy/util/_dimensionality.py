# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict
from .pyutil import defaultnamedtuple

dimension_codes = OrderedDict(zip(
    'length mass time current temperature amount'.split(),  # not considering luminous_intensity
    'L M T I Θ N'.split()
))


class DimensionalitySI(defaultnamedtuple('DimensionalitySIBase', dimension_codes.keys(), (0,)*len(dimension_codes))):

    def __mul__(self, other):
        return self.__class__(*(x+y for x, y in zip(self, other)))

    def __truediv__(self, other):
        return self.__class__(*(x-y for x, y in zip(self, other)))

    def __pow__(self, exp):
        return self.__class__(*(x*exp for x in self))

base_registry = {name: DimensionalitySI(**{name: 1}) for name in dimension_codes}
