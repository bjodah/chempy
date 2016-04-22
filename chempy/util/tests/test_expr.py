# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math

import pytest

from chempy import Substance
from chempy.units import (
    units_library, default_constants, Backend, to_unitless,
    SI_base_registry, default_units as u
)
from chempy.util.testing import requires

from ..expr import Expr, mk_Poly
from ..parsing import parsing_library


def _get_cv(kelvin=1, gram=1, mol=1):
    class HeatCapacity(Expr):
        parameter_keys = ('temperature',)
        kw = {'substance': None}

    class EinsteinSolid(HeatCapacity):
        """ arguments: einstein temperature """
        nargs = 1

        def __call__(self, variables, args=None, backend=math):
            molar_mass = self.substance.mass
            TE = self.arg(variables, args, 0)  # einstein_temperature
            R = variables['R']
            T = variables['temperature']
            # Canoncial ensemble:
            molar_c_v = 3*R*(TE/(2*T))**2 * backend.sinh(TE/(2*T))**-2
            return molar_c_v/molar_mass

    Al = Substance.from_formula('Al', other_properties={'DebyeT': 428*kelvin})
    Be = Substance.from_formula('Be', other_properties={'DebyeT': 1440*kelvin})
    Al.mass *= gram/mol
    Be.mass *= gram/mol

    def einT(s):
        return 0.806*s.other_properties['DebyeT']
    return {s.name: EinsteinSolid([einT(s)], substance=s) for s in (Al, Be)}


@requires(parsing_library)
def test_Expr():
    cv = _get_cv()
    _ref = 0.8108020083055849
    assert abs(cv['Al']({'temperature': 273.15, 'R': 8.3145}) - _ref) < 1e-14


def test_nargs():
    class A(Expr):
        nargs = 1

    with pytest.raises(ValueError):
        A([1, 2])


@requires('sympy')
def test_Expr_symbolic():
    import sympy
    cv = _get_cv()
    R, T = sympy.symbols('R T')
    sexpr = cv['Be']({'temperature': T, 'R': R}, backend=sympy)
    assert sexpr.free_symbols == set([T, R])


@requires(units_library)
def test_Expr_units():
    cv = _get_cv(u.kelvin, u.gram, u.mol)
    R = default_constants.molar_gas_constant.rescale(u.joule/u.mol/u.kelvin)

    def _check(T=273.15*u.kelvin):
        result = cv['Be']({'temperature': T, 'R': R}, backend=Backend())
        ref = 0.7342617587256584*u.joule/u.gram/u.kelvin
        assert abs(to_unitless((result - ref)/ref)) < 1e-10
    _check()
    _check(491.67*u.rankine)


@requires(units_library)
def test_Expr__dedimensionalisation():
    cv = _get_cv(u.kelvin, u.gram, u.mol)
    units, expr = cv['Be']._dedimensionalisation(SI_base_registry)
    assert units == [u.kelvin]
    assert expr.args == [0.806*1440]


def test_mk_Poly():
    Poly = mk_Poly('T', reciprocal=True)
    p = Poly([3, 2, 5, 7, 8, 2, 9])
    assert p.eval_poly({'T': 13}) == 2.57829
