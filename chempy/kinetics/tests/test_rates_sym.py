# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math

import pytest

from chempy.symbolic import (
    evaluate, concentrations as c_reg,
    parameters as p_reg
)

from chempy import Reaction, ReactionSystem
from chempy.units import units_library, default_units as u
from chempy.util.testing import requires
from ..rates_sym import RateExpr

class SpecialFraction(RateExpr):
    """
    Example from Atkins, De Paula, Physical Chemistry
    H2 + Br2 -> 2HBr
    ½ dHBr/dt = k[H2][Br2]**(3/2) / ([Br2] + k'[HBr])
    """
    def __call__(self, variables, args=None, backend=math):
        args = args or self.args
        two = 2 * backend.pi**0
        k = variables.get(self.arg_keys[0], args[0])
        kp = variables.get(self.arg_keys[0], args[1])
        H2, Br2, HBr = variables['H2'], variables['Br2'], variables['HBr']
        return k*H2*Br2**(3/two) / (Br2 + kp*HBr)

def _get_SpecialFraction_rateexpr():
    ratex = SpecialFraction([11, 13], rxn)
    rxn = Reaction({'H2': 1, 'Br2': 1}, {'HBr': 2}, ratex)


def _get_h2_br2_rsys(k, kprime):
    # rateexpr = _get_h2_br2_rateexpr(c_reg, k, kprime)
    # rxn = Reaction({'H2': 1, 'Br2': 1}, {'HBr': 2}, rateexpr)
    def _rateexpr(conc, state, rxn_params, backend=math):
        k, kp = rxn_params.get('k', 1.13), rxn_params.get("k'", 2.14)
        return k*c['H2']*c['Br2']**(3/be.S(2)) / (c['Br2'] + kp*c['HBr'])

    rxn = Reaction({'H2': 1, 'Br2': 1}, {'HBr': 2}, _rateexpr)
    rsys = ReactionSystem([rxn], 'H2 Br2 HBr')
    return rsys


@requires('sympy')
def test_symbolic_param():
    k = p_reg("k", positive=True)
    kprime = p_reg("k'", positive=True)
    r = _get_h2_br2_rateexpr(c_reg, k, kprime)


@requires('numpy', 'sympy')
def test_symbolic_expr():
    rsys = _get_h2_br2_rsys(11, 13)
    p = rsys.rxns[0].param

    conc = {'H2': 2, 'Br2': 3, 'HBr': 5}

    def _check(_k, kprime, c, params=None):
        ref = k*c['H2']*c['Br2']**1.5/(c['Br2'] + kprime*c['HBr'])
        val = evaluate(q, c, [c_reg])
        assert abs(val - ref) < 1e-15

    _check(11, 13, conc)


# @requires(units_library)
# def test_Quotient__units():
#     _k, _kprime = 3.5 * u.s**-1 * u.molar**-0.5, 9.2
#     rsys = _get_h2_br2_rsys(_k, _kprime)
#     q = rsys.rxns[0].param
#     conc = {'H2': 2*u.molar, 'Br2': 3*u.molar, 'HBr': 5*u.molar}

#     def _check(k, kprime, c, params=None):
#         ref = k*c['H2']*c['Br2']**1.5/(c['Br2'] + kprime*c['HBr'])
#         val = q.eval(rsys, 0, rsys.as_per_substance_array(c, unit=u.molar),
#                      params=params)
#         assert abs(val - ref) < 1e-15
#     _check(_k, _kprime, conc)
#     with pytest.raises(ValueError):
#         _check(_k, _kprime*u.second, conc)
#     _check(_k, _kprime, conc, [_k, 1, _kprime])
#     alt_k, alt_kprime = 2*u.s**-1*u.molar**-0.5, 5
#     _check(alt_k, alt_kprime, conc, [alt_k, 1, alt_kprime])

#     q2 = q.rebuild()
#     assert q2.get_params() == [_k, 1, _kprime]

#     q3 = q.rebuild([1, 2, 3])
#     assert q3.get_params() == [1, 2, 3]


def _get_ExpReciprocalT_rsys_1(A=1e10, Ea_over_R=40e3/8.3145):
    args = (A, -Ea_over_R)
    ert = ExpReciprocalT(args)
    assert tuple(ert.get_params()) == args

    rxn = Reaction({'A': 1}, {'B': 1}, ert)
    rsys = ReactionSystem([rxn], 'A B')
    return rsys


def test_ExpReciprocalT():
    rsys = _get_ExpReciprocalT_rsys_1()
    ert = rsys.rxns[0].param
    ref = 1e10 * math.exp(-40e3/(8.3145*298.15))
    res = ert.eval(rsys, 0, [3.0], global_p={'T': 298.15})
    assert abs((res - ref)/ref) < 1e-14
