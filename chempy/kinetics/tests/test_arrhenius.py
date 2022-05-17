# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function)

import math

from chempy.chemistry import Reaction
from chempy.util.testing import requires
from chempy.units import (
    units_library, allclose, SymPyDeDim, to_unitless,
    default_units as u, default_constants as dc,
    UncertainQuantity as UQ

)
from ..arrhenius import (
    arrhenius_equation, ArrheniusParam, ArrheniusParamWithUnits
)


def test_arrhenius_equation():
    assert abs(arrhenius_equation(3, 831.4472, 100) - 3/2.7182818) < 1e-7

_A1, _Ea1, _T1, _k1 = 1e10, 42e3, 273.15, 1e10 * math.exp(-42e3/(8.3145*273.15))


def test_ArrheniusParam():
    k = ArrheniusParam(_A1, _Ea1)(_T1)
    assert abs((k - _k1)/_k1) < 1e-4


@requires('numpy')
def test_ArrheniusParam__from_rateconst_at_T():
    ap = ArrheniusParam.from_rateconst_at_T(_Ea1, (_T1, _k1))
    assert abs((ap.A - _A1)/_A1) < 1e-4


def _get_ref2_units():
    A__s = 1e10
    act_J__mol = 42e3
    freezing_K = 273.15

    class ValueHolder:
        A = A__s / u.s
        Ea = act_J__mol * u.J/u.mol
        T = freezing_K * u.K
        k = A__s/u.s * math.exp(-act_J__mol/(8.3145*freezing_K))

    return ValueHolder()


@requires(units_library)
def test_ArrheniusParamWithUnits():
    _2 = _get_ref2_units()
    ap = ArrheniusParamWithUnits(_2.A, _2.Ea)

    k = ap(_2.T)
    assert abs((k - _2.k)/_2.k) < 1e-4

    r = Reaction({'H2O2': 1}, {'OH': 2}, ap)
    ratc = r.rate_expr().rate_coeff({'temperature': _2.T})
    assert allclose(ratc, _2.k, rtol=1e-4)


@requires(units_library)
def test_ArrheniusParamWithUnits__from_rateconst_at_T():
    _2 = _get_ref2_units()
    apu = ArrheniusParamWithUnits.from_rateconst_at_T(_2.Ea, (_2.T, _2.k))
    assert allclose(apu(_2.T), _2.k)


def test_ArrheniusParamWithUnits__UncertainQuantity():
    ap = ArrheniusParamWithUnits(UQ(1e10, 1/u.molar/u.s, 0.2e10), UQ(32e3, u.joule/u.mole, 2e3))
    sym_consts = namedtuple("SymConstants", ["molar_gas_constant"])(
        molar_gas_constant=np.array(sympy.Symbol('R_J_Kmol'), dtype='object')*u.J/u.K/u.mol
    )
    sym_vars = {
        'temperature': np.array(sympy.Symbol('T_K'), dtype='object')*u.K
        'A': np.array(1.0*sympy.Symbol(f'[{sk}]'), dtype=object)*u.molar,
    }
    const_vars = [
        (dc, {'temperature': 298*u.K, 'A': 2*u.molar}),
        (sym_consts, sym_vars)
    ]
    for consts, variables in const_vars:
        re = ap.as_RateExpr(constants=consts, units=u)
        rx = Reaction({'A': 2}, {'B': 1})
        sympy_dedim = SymPyDeDim()
        rate = re(variables, backend=sympy_dedim, reaction=rx)
        ref = 2*2*1e10*math.exp(-32e3/298/8.314511)
        assert allclose(to_unitless(rate, u.molar/u.s), ref, rtol=1e-3)
