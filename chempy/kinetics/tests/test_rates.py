# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math

from chempy import Reaction, ReactionSystem
from chempy.units import to_unitless, units_library, default_units as u
from chempy.util.testing import requires
from ..rates import RateExpr, MassAction, ArrheniusMassAction, Radiolytic


class SpecialFraction(RateExpr):
    """
    Example from Atkins, De Paula, Physical Chemistry
    H2 + Br2 -> 2HBr
    ½ dHBr/dt = k[H2][Br2]**(3/2) / ([Br2] + k'[HBr])
    """
    def __call__(self, variables, args=None, backend=math):
        args = args or self.args
        two = 2 * backend.pi**0
        k = self.arg(variables, args, 0)
        kp = self.arg(variables, args, 1)
        H2, Br2, HBr = variables['H2'], variables['Br2'], variables['HBr']
        return k*H2*Br2**(3/two) / (Br2 + kp*HBr)


def _get_SpecialFraction_rsys(k, kprime):
    ratex = SpecialFraction([k, kprime])
    rxn = Reaction({'H2': 1, 'Br2': 1}, {'HBr': 2}, ratex)
    return ReactionSystem([rxn], 'H2 Br2 HBr')


@requires('numpy')
def test_specialfraction_rateexpr():
    rsys = _get_SpecialFraction_rsys(11, 13)
    p = rsys.rxns[0].param

    conc = {'H2': 2, 'Br2': 3, 'HBr': 5}

    def _check(k, kprime, c):
        ref = k*c['H2']*c['Br2']**1.5/(c['Br2'] + kprime*c['HBr'])
        assert abs(p(c) - ref) < 1e-15

    _check(11, 13, conc)


def test_MassAction():
    ma = MassAction([3.14])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, ma, {'B': 1})
    assert abs(ma({'A': 11, 'B': 13, 'C': 17}) - 3.14*13*11**2) < 1e-14
    assert abs(ma({'A': 11, 'B': 13, 'C': 17}, [2.72]) - 2.72*13*11**2) < 1e-12


def test_ArrheniusMassAction():
    A, Ea_over_R = 1.2e11, 40e3/8.3145
    ama = ArrheniusMassAction([A, Ea_over_R])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, ama, {'B': 1})
    T_ = 'temperature'

    def ref(v):
        return 1.2e11*math.exp(-Ea_over_R/v[T_])*v['B']*v['A']**2

    for params in [(11., 13., 17., 311.2),
                   (12, 8, 5, 270)]:
        var = dict(zip(['A', 'B', 'C', T_], params))
        r = ref(var)
        assert abs((ama(var) - r)/r) < 1e-14


@requires(units_library)
def test_specialfraction_rateexpr__units():
    _k, _kprime = 3.5 * u.s**-1 * u.molar**-0.5, 9.2
    rsys = _get_SpecialFraction_rsys(_k, _kprime)
    p = rsys.rxns[0].param

    conc = {'H2': 2*u.molar, 'Br2': 3000*u.mol/u.m**3, 'HBr': 5*u.molar}

    def _check(k, kprime, c, args=None):
        ref = k*c['H2']*c['Br2']**1.5/(c['Br2'] + kprime*c['HBr'])
        assert abs(p(c, args) - ref) < 1e-15*u.molar/u.second

    _check(_k, _kprime, conc)
    alt_k, alt_kprime = 2*u.s**-1*u.molar**-0.5, 5
    _check(alt_k, alt_kprime, conc, [alt_k, alt_kprime])


@requires(units_library)
def test_ArrheniusMassAction__units():
    A, Ea_over_R = 1.2e11/u.molar**2/u.second, 40e3/8.3145*u.kelvin
    ama = ArrheniusMassAction([A, Ea_over_R])
    Reaction({'A': 2, 'B': 1}, {'C': 1}, ama, {'B': 1})
    T_ = 'temperature'

    def ref(v):
        return 1.2e11/u.molar**2/u.second*math.exp(
            -Ea_over_R/v[T_])*v['B']*v['A']**2

    for params in [(11.*u.molar, 13.*u.molar, 17.*u.molar, 311.2*u.kelvin),
                   (12*u.molar, 8*u.molar, 5*u.molar, 270*u.kelvin)]:
        var = dict(zip(['A', 'B', 'C', T_], params))
        r = ref(var)
        assert abs((ama(var) - r)/r) < 1e-14


def test_Radiolytic():
    r = Radiolytic([2.1e-7])
    res = r({'doserate': 0.15, 'density': 0.998})
    assert abs(res - 0.15*0.998*2.1e-7) < 1e-15


@requires(units_library)
def test_Radiolytic__units():

    def _check(r):
        res = r({'doserate': 0.15*u.gray/u.second, 'density': 0.998*u.kg/u.decimetre**3})
        ref = 0.15*0.998*2.1e-7*u.molar/u.second
        assert abs(to_unitless((res - ref)/ref)) < 1e-15

    _check(Radiolytic([2.1e-7*u.mol/u.joule]))
    _check(Radiolytic([2.0261921896167396*u.per100eV]))
