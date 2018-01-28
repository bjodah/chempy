# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math
from operator import add
from functools import reduce

import pytest

from chempy import Substance
from chempy.units import (
    allclose, units_library, default_constants, Backend, to_unitless,
    SI_base_registry, default_units as u
)
from ..testing import requires
from ..pyutil import defaultkeydict
from .._expr import Expr, mk_Poly, mk_PiecewisePoly, create_Piecewise, create_Poly, Log10
from ..parsing import parsing_library


class HeatCapacity(Expr):
    parameter_keys = ('temperature',)


class EinsteinSolid(HeatCapacity):
    parameter_keys = HeatCapacity.parameter_keys + ('molar_gas_constant',)
    argument_names = ('einstein_temperature', 'molar_mass')

    def __call__(self, variables, backend=math):
        TE, molar_mass = self.all_args(variables, backend=backend)  # einstein_temperature
        T, R = self.all_params(variables, backend=backend)
        # Canonical ensemble:
        molar_c_v = 3*R*(TE/(2*T))**2 * backend.sinh(TE/(2*T))**-2
        return molar_c_v/molar_mass


def _get_cv(kelvin=1, gram=1, mol=1):

    Al = Substance.from_formula('Al', data={'DebyeT': 428*kelvin})
    Be = Substance.from_formula('Be', data={'DebyeT': 1440*kelvin})

    def einT(s):
        return 0.806*s.data['DebyeT']
    return {s.name: EinsteinSolid([einT(s), s.mass * gram/mol]) for s in (Al, Be)}


@requires(parsing_library)
def test_Expr():
    cv = _get_cv()
    _ref = 0.8108020083055849
    assert abs(cv['Al']({'temperature': 273.15, 'molar_gas_constant': 8.3145}) - _ref) < 1e-14


def _poly(args, x, backend=None):
    x0, coeffs = args[0], args[1:]
    return reduce(add, [c*(x-x0)**i for i, c in enumerate(coeffs)])


@requires(parsing_library)
def test_Expr__nested_Expr():
    Poly = Expr.from_callback(_poly, parameter_keys=('x',), argument_names=('x0', Ellipsis))
    T = Poly([3, 7, 5])

    cv = _get_cv()
    _ref = 0.8108020083055849
    args = {'temperature': T, 'x': (273.15-7)/5 + 3, 'molar_gas_constant': 8.3145}
    assert abs(cv['Al'](args) - _ref) < 1e-14
    Al2 = cv['Al']/2
    assert abs(Al2(args) - _ref/2) < 1e-14


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
    sexpr = cv['Be']({'temperature': T, 'molar_gas_constant': R}, backend=sympy)
    assert sexpr.free_symbols == set([T, R])


@requires(units_library)
def test_Expr_units():
    cv = _get_cv(u.kelvin, u.gram, u.mol)
    R = default_constants.molar_gas_constant.rescale(u.joule/u.mol/u.kelvin)

    def _check(T=273.15*u.kelvin):
        result = cv['Be']({'temperature': T, 'molar_gas_constant': R}, backend=Backend())
        ref = 0.7342617587256584*u.joule/u.gram/u.kelvin
        assert abs(to_unitless((result - ref)/ref)) < 1e-10
    _check()
    _check(491.67*u.rankine)


@requires(units_library)
def test_Expr_dedimensionalisation__1():
    cv = _get_cv(u.kelvin, u.gram, u.mol)
    units, expr = cv['Be'].dedimensionalisation(SI_base_registry)
    assert units == [u.kelvin, u.kg/u.mol]
    assert abs(expr.args[0] - 0.806*1440) < 1e-14
    assert abs(expr.args[1] - 9.01218e-3) < 1e-7


@requires(units_library)
def test_Expr_dedimensionalisation__2():
    Poly = Expr.from_callback(_poly, parameter_keys=('E',), argument_names=('x0', Ellipsis))
    T = Poly([3*u.J, 7*u.K, 5*u.K/u.J])
    T = Poly([0.7170172084130019*u.cal, 12.6*u.Rankine, 5*u.K/u.J])

    _ref = 0.8108020083055849  # Al at 273.15 K with R=8.3145
    cv_Al = _get_cv(u.kelvin, u.gram, u.mol)['Al']

    assert isinstance(cv_Al, EinsteinSolid)
    assert cv_Al.args[0] == 0.806*428*u.kelvin
    assert abs(cv_Al({'temperature': 273.15*u.K, 'molar_gas_constant': 8.3145*u.J/u.K/u.mol}) -
               _ref*u.J/u.gram/u.kelvin) < 1e-14

    cv_Al_units, Al_dedim = cv_Al.dedimensionalisation(SI_base_registry)
    assert allclose(cv_Al_units, [u.K, u.kg/u.mol])
    assert isinstance(Al_dedim, EinsteinSolid)
    T_units, dT = T.dedimensionalisation(SI_base_registry)
    assert allclose(T_units, [u.J, u.K, u.K/u.J])
    assert allclose(dT.args, [3, 7, 5])
    assert abs(Al_dedim({'temperature': 273.15, 'molar_gas_constant': 8.3145}) - _ref*1000) < 1e-14
    assert abs(Al_dedim({'temperature': dT, 'E': (273.15-7)/5 + 3, 'molar_gas_constant': 8.3145}) - _ref*1000) < 1e-14


@requires(units_library)
def test_Expr_dedimensionalisation__nested():
    Poly = Expr.from_callback(_poly, parameter_keys=('E',), argument_names=('x0', Ellipsis))
    TE = Poly([3*u.J, 7*u.K, 5*u.K/u.J])
    TE = Poly([0.7170172084130019*u.cal, 12.6*u.Rankine, 5*u.K/u.J]) * 0.806*428/273.15

    _ref = 0.8108020083055849  # Al at 273.15 K with R=8.3145
    cv_Al = _get_cv(u.kelvin, u.gram, u.mol)['Al']
    cv_Al_units, Al_dedim = cv_Al.dedimensionalisation(SI_base_registry, {'einstein_temperature': TE})
    assert abs(Al_dedim({'temperature': 273.15, 'E': (273.15-7)/5 + 3, 'molar_gas_constant': 8.3145}) -
               _ref*1000) < 1e-14


def test_Expr__from_callback():
    def two_dim_gauss(args, x, y, backend=None):
        A, x0, y0, sx, sy = args
        xp, yp = x-x0, y-y0
        vx, vy = 2*sx**2, 2*sy**2
        return A*backend.exp(-(xp**2/vx + yp**2/vy))

    TwoDimGauss = Expr.from_callback(two_dim_gauss, parameter_keys=('x', 'y'), nargs=5)
    with pytest.raises(ValueError):
        TwoDimGauss([1, 2])
    args = [3, 2, 1, 4, 5]
    g1 = TwoDimGauss(args)
    ref = two_dim_gauss(args, 6, 7, math)
    assert abs(g1({'x': 6, 'y': 7}) - ref) < 1e-14


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

        def __call__(self, variables, backend=None):
            p0, p1 = self.all_args(variables)
            return p0 + p1*variables['x']

    l1 = Linear([3, 2])
    assert l1(dict(x=5)) == 13
    with pytest.raises(ValueError):
        Linear([3])
    with pytest.raises(ValueError):
        Linear([3, 2, 1])

    l2 = Linear([3, 2], ['a', 'b'])
    # with pytest.raises(ValueError):
    #     Linear([3, 2], ['a'])
    with pytest.raises(ValueError):
        Linear([3, 2], ['a', 'b', 'c'])

    assert l2(dict(x=5)) == 13
    assert l2(dict(x=5, a=11, b=13)) == 11 + 13*5


def test_PiecewisePoly():
    Poly = mk_Poly('temperature')

    p1 = Poly([0, 1, 0.1])
    assert p1.eval_poly({'temperature': 10}) == 2

    p2 = Poly([0, 3, -.1])
    assert p2.eval_poly({'temperature': 10}) == 2

    TPiecewisePoly = mk_PiecewisePoly('temperature')
    tpwp = TPiecewisePoly.from_polynomials([(0, 10), (10, 20)], [p1, p2])
    assert tpwp.eval_poly({'temperature': 5}) == 1.5
    assert tpwp.eval_poly({'temperature': 15}) == 1.5
    assert tpwp.parameter_keys == ('temperature',)

    with pytest.raises(ValueError):
        tpwp.eval_poly({'temperature': 21})


def test_create_Piecewise_Poly():
    PolyT = create_Poly('Tmpr')

    p1 = PolyT([1, 0.1])
    assert p1({'Tmpr': 10}) == 2

    p2 = PolyT([3, -.1])
    assert p2({'Tmpr': 10}) == 2

    PiecewiseT = create_Piecewise('Tmpr')
    pw = PiecewiseT([0, p1, 10, p2, 20])
    assert pw({'Tmpr': 5}) == 1.5
    assert pw({'Tmpr': 15}) == 1.5
    assert pw.parameter_keys == ('Tmpr',)

    with pytest.raises(ValueError):
        pw({'Tmpr': 21})


@requires('sympy')
def test_PiecewisePoly__sympy():
    import sympy as sp
    Poly = mk_Poly('T')
    p1 = Poly([0, 1, 0.1])
    p2 = Poly([0, 3, -.1])

    TPiecewisePoly = mk_PiecewisePoly('temperature')
    tpwp = TPiecewisePoly([2, 2, 0, 10, 2, 10, 20, 0, 1, 0.1, 0, 3, -.1])
    x = sp.Symbol('x')
    res = tpwp.eval_poly({'temperature': x}, backend=sp)
    assert isinstance(res, sp.Piecewise)
    assert res.args[0][0] == 1 + 0.1*x
    assert res.args[0][1] == sp.And(0 <= x, x <= 10)
    assert res.args[1][0] == 3 - 0.1*x
    assert res.args[1][1] == sp.And(10 <= x, x <= 20)

    with pytest.raises(ValueError):
        tpwp.from_polynomials([(0, 10), (10, 20)], [p1, p2])


@requires('sympy')
def test_create_Piecewise_Poly__sympy():
    import sympy as sp
    Poly = create_Poly('Tmpr')
    p1 = Poly([1, 0.1])
    p2 = Poly([3, -.1])

    TPw = create_Piecewise('Tmpr')
    pw = TPw([0, p1, 10, p2, 20])
    x = sp.Symbol('x')
    res = pw({'Tmpr': x}, backend=sp)
    assert isinstance(res, sp.Piecewise)
    assert res.args[0][0] == 1 + 0.1*x
    assert res.args[0][1] == sp.And(0 <= x, x <= 10)
    assert res.args[1][0] == 3 - 0.1*x
    assert res.args[1][1] == sp.And(10 <= x, x <= 20)


@requires('sympy')
def test_create_Piecewise__nan_fallback__sympy():
    import sympy as sp

    TPw = create_Piecewise('Tmpr', nan_fallback=True)
    pw = TPw([0, 42, 10, 43, 20])
    x = sp.Symbol('x')
    res = pw({'Tmpr': x}, backend=sp)
    assert isinstance(res, sp.Piecewise)
    assert res.args[0][0] == 42
    assert res.args[0][1] == sp.And(0 <= x, x <= 10)
    assert res.args[1][0] == 43
    assert res.args[1][1] == sp.And(10 <= x, x <= 20)
    assert res.args[2][0].name.lower() == 'nan'
    assert res.args[2][1] == True  # noqa


def test_BinaryExpr():
    Poly = Expr.from_callback(_poly, parameter_keys=('x',), argument_names=('x0', Ellipsis))
    p1 = Poly([1, 2, 3])
    p2 = Poly([2, 3, 4])
    assert p1({'x': 5}) == 14
    assert p2({'x': 5}) == 15
    assert (p1+p2)({'x': 5}) == 14+15
    assert (p1-p2)({'x': 5}) == 14-15
    assert (p1*p2)({'x': 5}) == 14*15
    assert (p1/p2)({'x': 5}) == 14/15

    assert (p1+2)({'x': 5}) == 14+2
    assert (p1-2)({'x': 5}) == 14-2
    assert (p1*2)({'x': 5}) == 14*2
    assert (p1/2)({'x': 5}) == 14/2

    assert (2+p1)({'x': 5}) == 2+14
    assert (2-p1)({'x': 5}) == 2-14
    assert (2*p1)({'x': 5}) == 2*14
    assert (2/p1)({'x': 5}) == 2/14

    assert p1 + 0 == p1
    assert p1 * 1 == p1
    assert p1 + p2 == p1 + p2
    assert p1 + p2*1 == p1 + p2 + 0

    assert -(-p1) == p1


@pytest.mark.slow
@requires('sympy')
def test_Expr__latex():
    Poly = Expr.from_callback(_poly, parameter_keys=('x',), argument_names=('x0', Ellipsis))
    p = Poly([1, 2, 3, 4])
    import sympy
    t = sympy.Symbol('t')
    ref = sympy.latex((2 + 3*(t-1) + 4*(t-1)**2).simplify())
    assert p.latex({'x': 't'}) == ref

    TE = Poly([3, 7, 5])
    cv_Al = _get_cv()['Al']
    T, E, R, m = sympy.symbols('T E R m')
    _TE = 7 + 5*(E-3)
    ref = sympy.latex(((3*R*(_TE/(2*T))**2 * sympy.sinh(_TE/(2*T))**-2)/m).simplify())
    cv_Al.unique_keys = ('TE_Al', 'm_Al')
    res = cv_Al.latex({'TE_Al': TE, 'temperature': 'T', 'x': 'E', 'molar_gas_constant': 'R', 'm_Al': 'm'})
    assert res == ref

    X = sympy.symbols('X')
    _TE2 = 7 + 5*(X-3)
    ref2 = sympy.latex(((3*R*(_TE2/(2*T))**2 * sympy.sinh(_TE2/(2*T))**-2)/m).simplify())
    res2 = cv_Al.latex({'TE_Al': TE, 'temperature': 'T', 'molar_gas_constant': 'R', 'm_Al': 'm'},
                       default=lambda k: k.upper())
    assert res2 == ref2


class Pressure1(Expr):
    argument_names = ('n',)
    parameter_keys = ('temperature', 'volume', 'R')

    def __call__(self, variables, backend=None):
        n, = self.all_args(variables, backend=backend)
        T, V, R = self.all_params(variables, backend=backend)
        return n*R*T/V


def test_Expr__single_arg():
    p = Pressure1(3)
    assert abs(p({'temperature': 273.15, 'volume': 0.17, 'R': 8.314}) - 3*8.314*273.15/0.17) < 1e-15


@requires(units_library)
def test_Expr__single_arg__units():
    p = Pressure1(3*u.mol)
    variables = {'temperature': 273.15*u.kelvin, 'volume': 170*u.dm3, 'R': 8.314*u.J/u.K/u.mol}
    assert allclose(p(variables), 3*8.314*273.15/0.17*u.Pa)


class Pressure2(Pressure1):
    def args_dimensionality(self):
        return ({'amount': 1},)


@requires(units_library)
def test_Expr__single_arg__units__dimensionality():
    p = Pressure2(unique_keys=('n1',))
    variables = {'temperature': 273.15*u.kelvin, 'volume': 170*u.dm3, 'R': 8.314*u.J/u.K/u.mol}
    assert allclose(p(dict(n1=3*u.mol, **variables)), 3*8.314*273.15/0.17*u.Pa)


def test_Expr__argument_defaults():

    class MyExpr(Expr):
        argument_names = ('a', 'b', 'c')
        argument_defaults = (17, 23)

        def __call__(self, variables={}, backend=math):
            a, b, c = self.all_args(variables, backend=backend)
            return a*b*c

    assert MyExpr([15])() == 15*17*23
    assert MyExpr([15, 17])() == 15*17*23
    assert MyExpr([15, 19])() == 15*19*23
    assert MyExpr([15, 19, 29])() == 15*19*29
    assert MyExpr(dict(zip('abc', [15, 19, 29])))() == 15*19*29


class MyK(Expr):
    argument_names = ('H', 'S')
    parameter_keys = ('T',)
    R = 8.3145

    def __call__(self, variables, backend=math):
        H, S = self.all_args(variables, backend=backend)
        T, = self.all_params(variables, backend=backend)
        return backend.exp(-(H - T*S)/(self.R*T))


def test_Expr__no_args():
    K1 = MyK(unique_keys=('H1', 'S1'))
    K2 = MyK(unique_keys=('H2', 'S2'))
    add = K1 + K2
    T = 298.15
    res = add({'H1': 2, 'H2': 3, 'S1': 5, 'S2': 7, 'T': T})
    RT = 8.3145 * 298.15
    ref = math.exp(-(2 - T*5)/RT) + math.exp(-(3 - T*7)/RT)
    assert abs(res - ref) < 1e-14


@requires(units_library)
def test_Expr__equality():
    K1 = MyK(unique_keys=('H1', 'S1'))
    K2 = MyK(unique_keys=('H2', 'S2'))
    assert K1 != K2
    assert K1 == K1
    K3 = MyK([23e3*u.J/u.mol, 42*u.J/u.mol/u.K])
    K4 = MyK([23e3, 42])
    K5 = MyK([24e3*u.J/u.mol, 42*u.J/u.mol/u.K])
    K6 = MyK([23e3, 43])
    assert K3 == K3
    assert K4 == K4
    assert K5 != K3
    assert K4 != K6
    assert K3 != K4


@requires('sympy')
def test_Expr__no_args__symbolic():
    K1 = MyK(unique_keys=('H1', 'S1'))
    K2 = MyK(unique_keys=('H2', 'S2'))
    add = K1 + K2
    import sympy
    v = defaultkeydict(sympy.Symbol)
    res = add(v, backend=sympy)
    R = 8.3145
    expr1 = sympy.exp(-(v['H1'] - v['T']*v['S1'])/R/v['T'])
    expr2 = sympy.exp(-(v['H2'] - v['T']*v['S2'])/R/v['T'])
    ref = expr1 + expr2
    assert (res - ref).simplify() == 0


class MyK2(Expr):
    argument_names = ('H', 'S', 'Cp', 'Tref')
    argument_defaults = (0, 298.15)
    parameter_keys = ('T')
    R = 8.3145

    def __call__(self, variables, backend=math):
        H, S, Cp, Tref = self.all_args(variables, backend=backend)
        T, = self.all_params(variables, backend=backend)
        _H = H + Cp*(T-Tref)
        _S = S + Cp*backend.log(T/Tref)
        return backend.exp(-(_H - T*_S)/(self.R*T))


def test_Expr__no_args__arg_defaults():
    K1 = MyK2(unique_keys=('H1', 'S1', 'Cp1'))
    K2 = MyK2(unique_keys=('H2', 'S2'))
    add = K1 + K2
    assert add.all_unique_keys() == set(['H1', 'H2', 'S1', 'S2', 'Cp1'])

    T = 293.15
    res = add({'H1': 2, 'H2': 3, 'S1': 5, 'S2': 7, 'T': T, 'Cp1': 13})
    RT = 8.3145 * T
    H1p = 2 + 13*(T - 298.15)
    S1p = 5 + 13*math.log(T/298.15)
    ref = math.exp(-(H1p - T*S1p)/RT) + math.exp(-(3 - T*7)/RT)
    assert abs(res - ref) < 1e-14


def test_create_Piecewise():
    PW = create_Piecewise('T')
    Ha, Sa, Hb, Sb, Ta, Tb = 40e3, -60, 37e3, -42, 293.15, 303.15
    a = MyK([Ha, Sa])
    b = MyK([Hb, Sb])
    pw = PW([273.15, a, 298.15, b, 323.15])  # 0, 25, 50 *C
    res_a = pw({'T': Ta})  # 20 *C
    res_b = pw({'T': Tb})  # 30 *C
    ref_a = math.exp(-(Ha - Ta*Sa)/(MyK.R*Ta))
    ref_b = math.exp(-(Hb - Tb*Sb)/(MyK.R*Tb))
    assert abs(res_a - ref_a) < 1e-14
    assert abs(res_b - ref_b) < 1e-14


def test_create_Poly():
    PolyT = create_Poly('T')
    p = PolyT([1, 2, 3, 4, 5])
    assert p({'T': 11}) == 1 + 2*11 + 3*11**2 + 4*11**3 + 5*11**4


def test_pow():
    PolyT = create_Poly('T')
    p = PolyT([1, 2, 3])
    p_to_two = p**2
    res = p_to_two({'T': 3})
    ref = (1 + 2*3 + 3*9)**2
    assert abs(res - ref) < 1e-12
    two_to_p = 2**p
    res = two_to_p({'T': 3})
    ref = 2**(1 + 2*3 + 3*9)
    assert abs(res - ref) < 1e-12


def test_functions():
    PolyT = create_Poly('T')
    p = PolyT([1, 2, 3])
    lgp = Log10(p)
    res = lgp({'T': 3})
    ref = math.log10(1 + 2*3 + 3*9)
    assert abs(res - ref) < 1e-12


def test_str_arg():
    PolyT = create_Poly('x')
    p = PolyT([1, 2, 3])
    x = Log10('T')
    res = p({'x': x, 'T': 1000})
    ref = 1 + 2*3 + 3*9
    assert abs(res - ref) < 1e-12


def test_implicit_str():
    PolyT = create_Poly('x')
    p = PolyT([1, 2, 3])
    expr1 = p / 'u'
    assert abs(expr1({'x': 3, 'u': 5}) - (1 + 2*3 + 3*9)/5) < 1e-12
    expr2 = 'u' / p
    assert abs(expr2({'x': 3, 'u': 5}) - 5/(1 + 2*3 + 3*9)) < 1e-12
    assert expr1.all_parameter_keys() == set(['x'])
    assert expr2.all_parameter_keys() == set(['x'])
