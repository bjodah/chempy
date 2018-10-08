# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from itertools import chain

import pytest

from ..util.testing import requires
from ..util.parsing import parsing_library
from ..units import default_units, units_library
from ..chemistry import Substance, Reaction
from ..reactionsystem import ReactionSystem


@requires(parsing_library, 'numpy')
def test_ReactionSystem():
    import numpy as np
    kw = dict(substance_factory=Substance.from_formula)
    r1 = Reaction.from_string('H2O -> H+ + OH-', 'H2O H+ OH-', name='r1')
    rs = ReactionSystem([r1], 'H2O H+ OH-', **kw)
    r2 = Reaction.from_string('H2O -> 2 H+ + OH-', 'H2O H+ OH-', name='r2')
    with pytest.raises(ValueError):
        ReactionSystem([r2], 'H2O H+ OH-', **kw)
    with pytest.raises(ValueError):
        ReactionSystem([r1, r1], 'H2O H+ OH-', **kw)
    assert rs.as_substance_index('H2O') == 0
    assert rs.as_substance_index(0) == 0
    varied, varied_keys = rs.per_substance_varied({'H2O': 55.4, 'H+': 1e-7, 'OH-': 1e-7},
                                                  {'H+': [1e-8, 1e-9, 1e-10, 1e-11], 'OH-': [1e-3, 1e-2]})
    assert varied_keys == ('H+', 'OH-')
    assert len(varied.shape) == 3
    assert varied.shape[:-1] == (4, 2)
    assert varied.shape[-1] == 3
    assert np.all(varied[..., 0] == 55.4)
    assert np.all(varied[:, 1, 2] == 1e-2)

    assert rs['r1'] is r1
    rs.rxns.append(r2)
    assert rs['r2'] is r2
    with pytest.raises(KeyError):
        rs['r3']
    rs.rxns.append(Reaction({}, {}, 0, name='r2', checks=()))
    with pytest.raises(ValueError):
        rs['r2']

    empty_rs = ReactionSystem([])
    rs2 = empty_rs + rs
    assert rs2 == rs
    rs3 = rs + empty_rs
    assert rs3 == rs


@requires(parsing_library)
def test_ReactionSystem__missing_substances_from_keys():
    r1 = Reaction({'H2O'}, {'H+', 'OH-'})
    with pytest.raises(ValueError):
        ReactionSystem([r1], substances={'H2O': Substance.from_formula('H2O')})
    kw = dict(missing_substances_from_keys=True, substance_factory=Substance.from_formula)
    rs = ReactionSystem([r1], substances={'H2O': Substance.from_formula('H2O')}, **kw)
    assert rs.substances['OH-'].composition == {0: -1, 1: 1, 8: 1}


@requires(parsing_library)
def test_ReactionSystem__check_balance():
    rs1 = ReactionSystem.from_string('\n'.join(['2 NH3 -> N2 + 3 H2', 'N2H4 -> N2 + 2 H2']))
    assert rs1.check_balance(strict=True)
    rs2 = ReactionSystem.from_string('\n'.join(['2 A -> B', 'B -> 2A']),
                                     substance_factory=Substance)
    assert not rs2.check_balance(strict=True)
    assert rs2.composition_balance_vectors() == ([], [])


def test_ReactionSystem__per_reaction_effect_on_substance():
    rs = ReactionSystem([Reaction({'H2': 2, 'O2': 1}, {'H2O': 2})])
    assert rs.per_reaction_effect_on_substance('H2') == {0: -2}
    assert rs.per_reaction_effect_on_substance('O2') == {0: -1}
    assert rs.per_reaction_effect_on_substance('H2O') == {0: 2}


def test_ReactionSystem__rates():
    rs = ReactionSystem([Reaction({'H2O'}, {'H+', 'OH-'}, 11)])
    assert rs.rates({'H2O': 3, 'H+': 5, 'OH-': 7}) == {'H2O': -11*3, 'H+': 11*3, 'OH-': 11*3}


def test_ReactionSystem__rates__cstr():
    k = 11
    rs = ReactionSystem([Reaction({'H2O2': 2}, {'O2': 1, 'H2O': 2}, k)])
    c0 = {'H2O2': 3, 'O2': 5, 'H2O': 53}
    fr = 7
    fc = {'H2O2': 13, 'O2': 17, 'H2O': 23}
    r = k*c0['H2O2']**2
    ref = {
        'H2O2': -2*r + fr*fc['H2O2'] - fr*c0['H2O2'],
        'O2': r + fr*fc['O2'] - fr*c0['O2'],
        'H2O': 2*r + fr*fc['H2O'] - fr*c0['H2O']
    }
    variables = dict(chain(c0.items(), [('fc_'+key, val) for key, val in fc.items()], [('fr', fr)]))
    for fck in (['fc_'+key for key in rs.substances], 'fc_'):
        assert rs.rates(variables, cstr_fr_fc=('fr', fck)) == ref


@requires('numpy')
def test_ReactionSystem__html_tables():
    r1 = Reaction({'A': 2}, {'A'}, name='R1')
    r2 = Reaction({'A'}, {'A': 2}, name='R2')
    rs = ReactionSystem([r1, r2])
    ut, unc = rs.unimolecular_html_table()
    assert unc == {0}
    from chempy.printing import html
    assert html(ut) == u'<table><tr><td>A</td><td ><a title="1: A → 2 A">R2</a></td></tr></table>'

    bt, bnc = rs.bimolecular_html_table()
    assert bnc == {1}
    assert html(bt) == u'<table><th></th><th>A</th>\n<tr><td>A</td><td ><a title="0: 2 A → A">R1</a></td></tr></table>'


@requires(parsing_library, 'numpy')
def test_ReactionSystem__substance_factory():
    r1 = Reaction.from_string('H2O -> H+ + OH-', 'H2O H+ OH-')
    rs = ReactionSystem([r1], 'H2O H+ OH-',
                        substance_factory=Substance.from_formula)
    assert rs.net_stoichs(['H2O']) == [-1]
    assert rs.net_stoichs(['H+']) == [1]
    assert rs.net_stoichs(['OH-']) == [1]
    assert rs.substances['H2O'].composition[8] == 1
    assert rs.substances['OH-'].composition[0] == -1
    assert rs.substances['H+'].charge == 1


@requires(units_library)
def test_ReactionSystem__as_per_substance_array_dict():
    mol = default_units.mol
    m = default_units.metre
    M = default_units.molar
    rs = ReactionSystem([], [Substance('H2O')])
    c = rs.as_per_substance_array({'H2O': 1*M}, unit=M)
    assert c.dimensionality == M.dimensionality
    assert abs(c[0]/(1000*mol/m**3) - 1) < 1e-16

    c = rs.as_per_substance_array({'H2O': 1})
    with pytest.raises(KeyError):
        c = rs.as_per_substance_array({'H': 1})

    assert rs.as_per_substance_dict([42]) == {'H2O': 42}


@requires(parsing_library)
def test_ReactionSystem__add():
    rs1 = ReactionSystem.from_string('\n'.join(['2 H2O2 -> O2 + 2 H2O', 'H2 + O2 -> H2O2']))
    rs2 = ReactionSystem.from_string('\n'.join(['2 NH3 -> N2 + 3 H2']))
    rs3 = rs1 + rs2
    assert rs1 == rs1
    assert rs1 != rs2
    assert rs3 != rs1
    assert len(rs1.rxns) == 2 and len(rs2.rxns) == 1 and len(rs3.rxns) == 3
    for k in 'H2O2 O2 H2O H2 NH3 N2'.split():
        assert k in rs3.substances
    rs1 += rs2
    assert len(rs1.rxns) == 3 and len(rs2.rxns) == 1
    assert rs1 == rs3

    rs4 = ReactionSystem.from_string("H2O -> H+ + OH-; 1e-4")
    rs4 += [Reaction({'H+', 'OH-'}, {'H2O'}, 1e10)]
    res = rs4.rates({'H2O': 1, 'H+': 1e-7, 'OH-': 1e-7})
    for k in 'H2O H+ OH-'.split():
        assert abs(res[k]) < 1e-16

    rs5 = ReactionSystem.from_string("H3O+ -> H+ + H2O")
    rs6 = rs4 + rs5
    rs7 = rs6 + (Reaction.from_string("H+ + H2O -> H3O+"),)
    assert len(rs7.rxns) == 4


@requires(parsing_library)
def test_ReactionSystem__from_string():
    rs = ReactionSystem.from_string('-> H + OH; Radiolytic(2.1e-7)', checks=())
    assert rs.rxns[0].reac == {}
    assert rs.rxns[0].prod == {'H': 1, 'OH': 1}
    assert rs.rxns[0].param.args == [2.1e-7]
    ref = 2.1e-7 * 0.15 * 998
    assert rs.rates({'doserate': .15, 'density': 998}) == {'H': ref, 'OH': ref}

    r2, = ReactionSystem.from_string("H2O + H2O + H+ -> H3O+ + H2O").rxns
    assert r2.reac == {'H2O': 2, 'H+': 1}
    assert r2.prod == {'H2O': 1, 'H3O+': 1}

    r3, = ReactionSystem.from_string('(H2O) -> e-(aq) + H+ + OH; Radiolytic(2.1e-7*mol/J)').rxns
    assert len(r3.reac) == 0 and r3.inact_reac == {'H2O': 1}
    assert r3.prod == {'e-(aq)': 1, 'H+': 1, 'OH': 1}
    from chempy.kinetics.rates import Radiolytic
    mol, J = default_units.mol, default_units.J
    assert r3.param == Radiolytic(2.1e-7*mol/J)
    assert r3.param != Radiolytic(2.0e-7*mol/J)
    assert r3.param != Radiolytic(2.1e-7)
    assert r3.order() == 0


@requires(parsing_library)
def test_ReactionSystem__from_string__string_rate_const():
    rsys = ReactionSystem.from_string("H+ + OH- -> H2O; 'kf'")
    r2, = rsys.rxns
    assert r2.reac == {'OH-': 1, 'H+': 1}
    assert r2.prod == {'H2O': 1}
    r2str = r2.string(rsys.substances, with_param=True)
    assert r2str.endswith('; kf')


@requires('numpy')
def test_ReactionSystem__upper_conc_bounds():
    rs = ReactionSystem.from_string('\n'.join(['2 NH3 -> N2 + 3 H2', 'N2H4 -> N2 +   2  H2']))
    c0 = {'NH3': 5, 'N2': 7, 'H2': 11, 'N2H4': 2}
    _N = 5 + 14 + 4
    _H = 15 + 22 + 8
    ref = {
        'NH3': min(_N, _H/3),
        'N2': _N/2,
        'H2': _H/2,
        'N2H4': min(_N/2, _H/4),
    }
    res = rs.as_per_substance_dict(rs.upper_conc_bounds(c0))
    assert res == ref


@requires('numpy')
def test_ReactionSystem__upper_conc_bounds__a_substance_no_composition():
    rs = ReactionSystem.from_string("""
    H2O -> e-(aq) + H2O+
    H2O+ + e-(aq) -> H2O
    """)
    c0 = {'H2O': 55.0, 'e-(aq)': 2e-3, 'H2O+': 3e-3}
    _O = 55 + 3e-3
    _H = 2*55 + 2*3e-3
    ref = {
        'H2O': min(_O, _H/2),
        'e-(aq)': float('inf'),
        'H2O+': min(_O, _H/2),
    }
    res = rs.as_per_substance_dict(rs.upper_conc_bounds(c0))
    assert res == ref


@requires(parsing_library)
def test_ReactionSystem__identify_equilibria():
    rsys = ReactionSystem.from_string("""
    2 H2 +  O2 -> 2 H2O     ; 1e-3
           H2O -> H+ + OH-  ; 1e-4/55.35
      H+ + OH- -> H2O       ; 1e10
         2 H2O -> 2 H2 + O2
    """)
    assert rsys.identify_equilibria() == [(0, 3), (1, 2)]


@requires(parsing_library)
def test_ReactionSystem__categorize_substances():
    rsys1 = ReactionSystem.from_string("""
    2 H2 +  O2 -> 2 H2O     ; 1e-3
           H2O -> H+ + OH-  ; 1e-4/55.35
      H+ + OH- -> H2O       ; 1e10
         2 H2O -> 2 H2 + O2
    """)
    assert all(not s for s in rsys1.categorize_substances().values())

    rsys2 = ReactionSystem.from_string('\n'.join(['2 NH3 -> N2 + 3 H2', 'N2H4 -> N2 +   2  H2']))
    assert rsys2.categorize_substances() == dict(accumulated={'N2', 'H2'}, depleted={'NH3', 'N2H4'},
                                                 unaffected=set(), nonparticipating=set())

    rsys3 = ReactionSystem.from_string("H+ + OH- -> H2O; 'kf'")
    assert rsys3.categorize_substances() == dict(accumulated={'H2O'}, depleted={'H+', 'OH-'},
                                                 unaffected=set(), nonparticipating=set())

    rsys4 = ReactionSystem([Reaction({'H2': 2, 'O2': 1}, {'H2O': 2})], 'H2 O2 H2O N2 Ar')
    assert rsys4.categorize_substances() == dict(accumulated={'H2O'}, depleted={'H2', 'O2'},
                                                 unaffected=set(), nonparticipating={'N2', 'Ar'})

    rsys5 = ReactionSystem.from_string("""
    A -> B; MassAction(unique_keys=('k1',))
    B + C -> A + C; MassAction(unique_keys=('k2',))
    2 B -> B + C; MassAction(unique_keys=('k3',))
    """, substance_factory=lambda formula: Substance(formula))
    assert rsys5.categorize_substances() == dict(accumulated={'C'}, depleted=set(),
                                                 unaffected=set(), nonparticipating=set())

    rsys6 = ReactionSystem.from_string("""H2O2 + Fe+3 + (H2O2) -> 2 H2O + O2 + Fe+3""")
    assert rsys6.rxns[0].order() == 2  # the additional H2O2 within parenthesis
    assert rsys6.categorize_substances() == dict(accumulated={'H2O', 'O2'}, depleted={'H2O2'},
                                                 unaffected={'Fe+3'}, nonparticipating=set())


@requires(parsing_library)
def test_ReactionSystem__split():
    a = """
    2 H2 +  O2 -> 2 H2O     ; 1e-3
           H2O -> H+ + OH-  ; 1e-4/55.35
      H+ + OH- -> H2O       ; 1e10
        2 H2O  -> 2 H2 + O2"""
    b = """
        2 N    -> N2"""
    c = """
        2 ClBr -> Cl2 + Br2
    """
    rsys1 = ReactionSystem.from_string(a+b+c)
    res = rsys1.split()
    ref = list(map(ReactionSystem.from_string, [a, b, c]))
    for rs in chain(res, ref):
        rs.sort_substances_inplace()
    res1a, res1b, res1c = res
    ref1a, ref1b, ref1c = ref
    assert res1a == ref1a
    assert res1b == ref1b
    assert res1c == ref1c
    assert res1c != ref1a
    assert rsys1.categorize_substances() == dict(
        accumulated={'N2', 'Cl2', 'Br2'}, depleted={'N', 'ClBr'},
        unaffected=set(), nonparticipating=set())


def test_ReactionSystem__subset():
    r1 = Reaction({'NH3': 2}, {'N2': 1, 'H2': 3})
    r2 = Reaction({'N2H4': 1}, {'N2': 1, 'H2': 2})
    rs1 = ReactionSystem([r1, r2])
    rs2 = rs1.subset(lambda r: 'N2H4' in r.keys())
    assert len(rs1.rxns) == 2 and len(rs2.rxns) == 1
    assert rs2 == ReactionSystem([r2])
