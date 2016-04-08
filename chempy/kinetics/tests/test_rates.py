# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import numpy as np

from chempy import Reaction, ReactionSystem
from ..rates import Quotient, Sum, GeneralPow
from ..ode import get_odesys


def test_MassAction():
    pass


def test_h2_br2():
    # Example from Atkins, De Paula, Physical Chemistry
    # H2 + Br2 -> 2HBr
    # ½ dHBr/dt = k[H2][Br2]**(3/2) / ([Br2] + k'[HBr])
    k = 3.142
    kprime = 2.718
    param = Quotient([
        GeneralPow([k], {'H2': 1, 'Br2': 3/2}),
        Sum([GeneralPow([1], {'Br2': 1}),
             GeneralPow([kprime], {'HBr': 1})])
    ])
    substances = 'H2 Br2 HBr'
    rxn = Reaction({'H2': 1, 'Br2': 1}, {'HBr': 2}, param)
    rsys = ReactionSystem([rxn], substances)

    odesys = get_odesys(rsys, include_params=True)
    c0 = {'H2': 13, 'Br2': 17, 'HBr': 19}
    r = k*c0['H2']*c0['Br2']**(3/2)/(c0['Br2'] + kprime*c0['HBr'])
    ref = rsys.as_per_substance_array({'H2': -r, 'Br2': -r, 'HBr': 2*r})
    res = odesys.f_cb(0, rsys.as_per_substance_array(c0))
    assert np.allclose(res, ref)
