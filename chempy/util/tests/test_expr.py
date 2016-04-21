# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math
from chempy import Substance
from chempy.util.testing import requires

from ..expr import Expr

def _get_cv():
    class HeatCapacity(Expr):
        parameter_keys = ('temperature',)
        kw = ('substance',)

    class EinsteinSolid(HeatCapacity):
        """ arguments: einstein temperature """
        def __call__(self, variables, args=None, backend=math):
            molar_mass = self.substance.mass
            eps = self.arg(variables, args, 0)  # einstein_temperature
            R = 8.3145
            T = variables['temperature']
            # Canoncial ensemble:
            molar_c_v = 3*R*(eps/(2*R*T))**2 * backend.sinh(eps/(2*R*T))**-2
            return molar_c_v/molar_mass

    Al = Substance.from_formula('Al', other_properties={'DebyeT': 428})
    Be = Substance.from_formula('Be', other_properties={'DebyeT': 1440})
    einT = lambda s: 0.806*s.other_properties['DebyeT']
    return {s.name: EinsteinSolid([einT(s)], substance=s) for s in (Al, Be)}

def test_Expr():
    cv = _get_cv()
    assert abs(cv['Al']({'temperature': 273.15}) - 0.9226900642408176) < 1e-14

@requires('sympy')
def test_Expr_symbolic():
    import sympy
    cv = _get_cv()
    sexpr = cv['Be']({'temperature': sympy.Symbol('T')}, backend=sympy)
    assert sexpr.free_symbols == set([sympy.Symbol('T')])
