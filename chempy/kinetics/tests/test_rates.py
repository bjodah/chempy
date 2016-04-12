# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import pytest

from chempy import Reaction, ReactionSystem
from chempy.units import default_units as u
from ..rates import Quotient, Sum, GeneralPow


def _get_h2_br2_param(k, kprime):
    # Example from Atkins, De Paula, Physical Chemistry
    # H2 + Br2 -> 2HBr
    # ½ dHBr/dt = k[H2][Br2]**(3/2) / ([Br2] + k'[HBr])
    return Quotient([
        GeneralPow([k], {'H2': 1, 'Br2': 3/2}),
        Sum([GeneralPow([1], {'Br2': 1}),
             GeneralPow([kprime], {'HBr': 1})])
    ])


def _get_h2_br2_rsys(k, kprime):
    param = _get_h2_br2_param(k, kprime)
    rxn = Reaction({'H2': 1, 'Br2': 1}, {'HBr': 2}, param)
    rsys = ReactionSystem([rxn], 'H2 Br2 HBr')
    return rsys


def test_Quotient():
    rsys = _get_h2_br2_rsys(11, 13)
    q = rsys.rxns[0].param

    assert q.get_params() == [11, 1, 13]

    conc = {'H2': 2, 'Br2': 3, 'HBr': 5}

    def _check(k, kprime, c, params=None):
        ref = k*c['H2']*c['Br2']**1.5/(c['Br2'] + kprime*c['HBr'])
        val = q.eval(rsys, 0, rsys.as_per_substance_array(c), params=params)
        assert abs(val - ref) < 1e-15

    _check(11, 13, conc)
    _check(11, 13, conc, [11, 1, 13])
    _check(2, 7, conc, [2, 1, 7])

    q2 = q.rebuild()
    assert q2.get_params() == [11, 1, 13]


def test_Quotient__units():
    _k, _kprime = 3.5 * u.s**-1 * u.molar**-0.5, 9.2
    rsys = _get_h2_br2_rsys(_k, _kprime)
    q = rsys.rxns[0].param
    conc = {'H2': 2*u.molar, 'Br2': 3*u.molar, 'HBr': 5*u.molar}

    def _check(k, kprime, c, params=None):
        ref = k*c['H2']*c['Br2']**1.5/(c['Br2'] + kprime*c['HBr'])
        val = q.eval(rsys, 0, rsys.as_per_substance_array(c, unit=u.molar), params=params)
        assert abs(val - ref) < 1e-15
    _check(_k, _kprime, conc)
    with pytest.raises(ValueError):
        _check(_k, _kprime*u.second, conc)
    _check(_k, _kprime, conc, [_k, 1, _kprime])
    alt_k, alt_kprime = 2*u.s**-1*u.molar**-0.5, 5
    _check(alt_k, alt_kprime, conc, [alt_k, 1, alt_kprime])

    q2 = q.rebuild()
    assert q2.get_params() == [_k, 1, _kprime]

    q3 = q.rebuild([1, 2, 3])
    assert q3.get_params() == [1, 2, 3]
