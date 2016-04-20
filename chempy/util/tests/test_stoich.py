# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

try:
    import numpy as np
except ImportError:
    np = None
else:
    from ..stoich import (
        get_coeff_mtx, decompose_yields,
    )
from chempy.units import units_library, allclose, _sum
from ..testing import requires


@requires('numpy')
def test_get_coeff_mtx():
    r = [
        ({'A': 1},       {'B': 1}),
        ({'A': 1, 'B': 1}, {'C': 2})
    ]
    A = get_coeff_mtx('ABC', r)
    Aref = np.array([
        [-1, -1],
        [1, -1],
        [0,  2]
    ])
    assert np.allclose(A, Aref)


@requires('numpy')
def test_decompose_yields_1():
    from chempy import Reaction

    gamma_yields = {
        'OH-': 0.5,
        'H2O2': 0.7,
        'OH': 2.7,
        'H2': 0.45,
        'H': 0.66,
        'H+': 3.1,
        'HO2': 0.02,
        'e-(aq)': 2.6,
    }

    rxns = [
        Reaction({'H2O': 1}, {'H+': 1, 'OH-': 1}),
        Reaction({'H2O': 1}, {'H+': 1, 'e-(aq)': 1, 'OH': 1}),
        Reaction({'H2O': 1}, {'H': 2, 'H2O2': 1}, inact_reac={'H2O': 1}),
        Reaction({'H2O': 1}, {'H2': 1, 'H2O2': 1}, inact_reac={'H2O': 1}),
        Reaction({'H2O': 1}, {'H2': 1, 'OH': 2}, inact_reac={'H2O': 1}),
        Reaction({'H2O': 1}, {'H2': 3, 'HO2': 2}, inact_reac={'H2O': 3}),
    ]

    k = decompose_yields(gamma_yields, rxns)
    k_ref = [0.5, 2.6, 0.33, 0.37, 0.05, 0.01]

    assert np.allclose(k, k_ref)

    G_H2O = sum(rxn.net_stoich(['H2O'])[0]*k[i] for
                i, rxn in enumerate(rxns))

    assert abs(G_H2O+4.64) < 1e-3


@requires(units_library)
def test_decompose_yields__units_1():
    from chempy import Reaction
    from chempy.units import default_units as u
    gamma_yields = {
        'OH-': 0.5*u.per100eV,
        'H2O2': 0.7*u.per100eV,
        'OH': 2.7*u.per100eV,
        'H2': 0.45*u.per100eV,
        'H': 0.66*u.per100eV,
        'H+': 3.1*u.per100eV,
        'HO2': 0.02*u.per100eV,
        'e-(aq)': 2.6*u.per100eV,
    }

    rxns = [
        Reaction({'H2O': 1}, {'H+': 1, 'OH-': 1}),
        Reaction({'H2O': 1}, {'H+': 1, 'e-(aq)': 1, 'OH': 1}),
        Reaction({'H2O': 1}, {'H': 2, 'H2O2': 1}, inact_reac={'H2O': 1}),
        Reaction({'H2O': 1}, {'H2': 1, 'H2O2': 1}, inact_reac={'H2O': 1}),
        Reaction({'H2O': 1}, {'H2': 1, 'OH': 2}, inact_reac={'H2O': 1}),
        Reaction({'H2O': 1}, {'H2': 3, 'HO2': 2}, inact_reac={'H2O': 3}),
    ]

    k = decompose_yields(gamma_yields, rxns)
    k_ref = [0.5, 2.6, 0.33, 0.37, 0.05, 0.01]*u.per100eV

    assert allclose(k, k_ref)

    G_H2O = [rxn.net_stoich(['H2O'])[0]*k[i] for i, rxn in enumerate(rxns)]
    ref = 4.64*u.per100eV
    assert abs((_sum(G_H2O)+ref)/ref) < 1e-3


@requires('numpy')
def test_decompose_yields_2():
    from chempy import Reaction
    yields = {'B': 3.0, 'C': 24.0}
    rxns = [
        Reaction({'A': 1}, {'B': 1, 'C': 1}, inact_reac={'A': 1}),
        Reaction({'A': 1}, {'C': 3})
    ]
    k = decompose_yields(yields, rxns)
    k_ref = [3, 7]

    rtol = 1e-12
    for a, b in zip(k, k_ref):
        assert abs(a-b) < abs(a*rtol)
