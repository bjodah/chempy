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

from ..expr import Expr, mk_Poly, PiecewisePoly
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
    assert p.parameter_keys == ('T',)


def test_Expr__nargs():

    class Linear(Expr):
        """ Arguments: p0, p1 """
        nargs = 2
        parameter_keys = ('x',)

        def __call__(self, variables, args=None, backend=None):
            p0, p1 = self.all_args(variables, args)
            return p0 + p1*variables['x']

    l1 = Linear([3, 2])
    assert l1(dict(x=5)) == 13
    with pytest.raises(ValueError):
        Linear([3])
    with pytest.raises(ValueError):
        Linear([3, 2, 1])

    l2 = Linear([3, 2], ['a', 'b'])
    with pytest.raises(ValueError):
        Linear([3, 2], ['a'])
    with pytest.raises(ValueError):
        Linear([3, 2], ['a', 'b', 'c'])

    assert l2(dict(x=5)) == 13
    assert l2(dict(x=5), [11, 13]) == 11 + 13*5


def test_PiecewisePoly():
    Poly = mk_Poly('temperature')

    p1 = Poly([0, 1, 0.1])
    assert p1.eval_poly({'temperature': 10}) == 2

    p2 = Poly([0, 3, -.1])
    assert p2.eval_poly({'temperature': 10}) == 2

    pw = PiecewisePoly([[(0, 10), (10, 20)], [p1, p2]])
    assert pw.eval_poly({'temperature': 5}) == 1.5
    assert pw.eval_poly({'temperature': 15}) == 1.5
    assert pw.parameter_keys == ('temperature',)

    with pytest.raises(ValueError):
        pw.eval_poly({'temperature': 21})


@requires('sympy')
def test_PiecewisePoly__sympy():
    import sympy as sp
    Poly = mk_Poly('T')
    p1 = Poly([0, 1, 0.1])
    p2 = Poly([0, 3, -.1])
    pwp = PiecewisePoly([[(0, 10), (10, 20)], [p1, p2]])
    x = sp.Symbol('x')
    res = pwp.eval_poly({'T': x}, backend=sp)
    assert isinstance(res, sp.Piecewise)
    assert res.args[0][0] == 1+0.1*x
    assert res.args[0][1] == sp.And(0 <= x, x <= 10)
    assert res.args[1][0] == 3-0.1*x
    assert res.args[1][1] == sp.And(10 <= x, x <= 20)
