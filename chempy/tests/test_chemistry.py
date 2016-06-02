# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from operator import attrgetter

import pytest

from ..util.testing import requires
from ..util.parsing import parsing_library
from ..units import default_units, units_library, to_unitless
from ..chemistry import (
    equilibrium_quotient, Substance, Species, Reaction, ReactionSystem,
    Equilibrium, balance_stoichiometry
)


@requires('numpy')
def test_equilibrium_quotient():
    assert abs(equilibrium_quotient([2.3, 3.7, 5.1], (-1, -1, 1)) -
               5.1/2.3/3.7) < 1e-14


@requires(parsing_library)
def test_Substance():
    s = Substance.from_formula('H+')
    assert s.composition == {0: 1, 1: 1}
    assert s.charge == 1
    assert abs(s.mass - 1.008) < 1e-3


def test_Substance__2():
    H2O = Substance(name='H2O',  charge=0, latex_name=r'$\mathrm{H_{2}O}$',
                    other_properties={'pKa': 14})  # will_be_missing_in='0.5.0', use data=...
    OH_m = Substance(name='OH-',  charge=-1, latex_name=r'$\mathrm{OH^{-}}$')
    assert sorted([OH_m, H2O], key=attrgetter('name')) == [H2O, OH_m]


@requires(parsing_library)
def test_Substance__from_formula():
    H2O = Substance.from_formula('H2O')
    assert H2O.composition == {1: 2, 8: 1}
    assert H2O.latex_name == 'H_{2}O'
    assert H2O.unicode_name == u'H₂O'
    assert H2O.html_name == u'H<sub>2</sub>O'


@requires(parsing_library)
def test_Species():
    s = Species.from_formula('H2O')
    assert s.phase_idx == 0
    mapping = {'(aq)': 0, '(s)': 1, '(g)': 2}
    assert Species.from_formula('CO2(g)').phase_idx == 3
    assert Species.from_formula('CO2(g)', mapping).phase_idx == 2
    assert Species.from_formula('CO2(aq)', mapping).phase_idx == 0
    assert Species.from_formula('NaCl(s)').phase_idx == 1
    assert Species.from_formula('NaCl(s)', phase_idx=7).phase_idx == 7
    assert Species.from_formula('CO2(aq)', mapping, phase_idx=7).phase_idx == 7


def test_Solute():
    from ..chemistry import Solute
    from ..util.pyutil import ChemPyDeprecationWarning
    with pytest.warns(ChemPyDeprecationWarning):
        w = Solute('H2O')
    assert w.name == 'H2O'


def test_Reaction():
    substances = s_Hp, s_OHm, s_H2O = (
        Substance('H+', composition={0: 1, 1: 1}),
        Substance('OH-', composition={0: -1, 1: 1, 8: 1}),
        Substance('H2O', composition={0: 0, 1: 2, 8: 1}),
    )
    substance_names = Hp, OHm, H2O = [s.name for s in substances]
    substance_dict = {n: s for n, s in zip(substance_names, substances)}
    r1 = Reaction({Hp: 1, OHm: 1}, {H2O: 1})
    assert sum(r1.composition_violation(substance_dict)) == 0
    assert r1.charge_neutrality_violation(substance_dict) == 0

    r2 = Reaction({Hp: 1, OHm: 1}, {H2O: 2})
    assert sum(r2.composition_violation(substance_dict)) != 0
    assert r2.charge_neutrality_violation(substance_dict) == 0

    r3 = Reaction({Hp: 2, OHm: 1}, {H2O: 2})
    assert sum(r3.composition_violation(substance_dict)) != 0
    assert r3.charge_neutrality_violation(substance_dict) != 0

    assert r3.keys() == {Hp, OHm, H2O}

    with pytest.raises(ValueError):
        Reaction({Hp: -1, OHm: -1}, {H2O: -1})

    assert r1 == Reaction({'H+', 'OH-'}, {'H2O'})

    r4 = Reaction({Hp, OHm}, {H2O}, 7)
    ref = {Hp: -3*5*7, OHm: -3*5*7, H2O: 3*5*7}
    r4.rate({Hp: 5, OHm: 3}) == ref
    r5 = r4.copy()
    assert r5 == r4
    assert r5 != r1


@requires(parsing_library)
def test_Reaction_parsing():
    r4 = Reaction({'H+': 2, 'OH-': 1}, {'H2O': 2}, 42.0)
    assert Reaction.from_string(str(r4), 'H+ OH- H2O') == r4
    assert Reaction.from_string(str(r4), None) == r4
    r5 = Reaction.from_string('2 H2O2 -> O2 + 2 H2O; 1e-7/molar/second', 'H2O O2 H2O2')
    assert to_unitless(r5.param, 1/default_units.molar/default_units.second) == 1e-7
    r6 = Reaction.from_string('->', checks=())
    assert r6.reac == {} and r6.prod == {}


@requires(parsing_library, units_library)
def test_Substance__molar_mass():
    mw_water = Substance.from_formula('H2O').molar_mass(default_units)
    q = mw_water / ((15.9994 + 2*1.008)*default_units.gram/default_units.mol)
    assert abs(q - 1) < 1e-3


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


def test_ReactionSystem__rates():
    rs = ReactionSystem([Reaction({'H2O'}, {'H+', 'OH-'}, 11)])
    assert rs.rates({'H2O': 3, 'H+': 5, 'OH-': 7}) == {'H2O': -11*3, 'H+': 11*3, 'OH-': 11*3}


def test_ReactionSystem__html_tables():
    r1 = Reaction({'A': 2}, {'A'}, name='R1')
    r2 = Reaction({'A'}, {'A': 2}, name='R2')
    rs = ReactionSystem([r1, r2])
    ut, unc = rs.unimolecular_html_table()
    assert unc == [r1]
    assert ut == u'<table><tr><td>A</td><td ><a title="A → 2 A">R2</a></td></tr></table>'

    bt, bnc = rs.bimolecular_html_table()
    assert bnc == [r2]
    assert bt == u'<table><th></th><th>A</th>\n<tr><td>A</td><td ><a title="2 A → A">R1</a></td></tr></table>'


@requires(parsing_library)
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


@requires(units_library)
def test_Equilibrium__as_reactions():
    s = default_units.second
    M = default_units.molar
    H2O, Hp, OHm = map(Substance, 'H2O H+ OH-'.split())
    eq = Equilibrium({'H2O': 1}, {'H+': 1, 'OH-': 1}, 1e-14)
    rate = 1.31e11/M/s
    fw, bw = eq.as_reactions(kb=rate, units=default_units)
    assert abs((bw.param - rate)/rate) < 1e-15
    assert abs((fw.param / M)/bw.param - 1e-14)/1e-14 < 1e-15


@requires(parsing_library)
def test_Reaction__from_string():
    r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", 'H2O H+ OH-'.split())
    assert r.reac == {'H2O': 1} and r.prod == {'H+': 1, 'OH-': 1}

    with pytest.raises(ValueError):
        Reaction.from_string("H2O -> H+ + OH-; 1e-4", 'H2O H OH-'.split())

    r2 = Reaction.from_string("H2O -> H+ + OH-; 1e-4; ref='important_paper'")
    assert r2.ref == 'important_paper'

    with pytest.raises(ValueError):
        Reaction.from_string("H2O -> H2O")
    Reaction.from_string("H2O -> H2O; None; checks=()")


@requires(parsing_library)
def test_ReactioN__latex():
    keys = 'H2O H2 O2'.split()
    subst = {k: Substance.from_formula(k) for k in keys}
    r2 = Reaction.from_string("2 H2O -> 2 H2 + O2", subst)
    assert r2.latex(subst) == r'2 H_{2}O \rightarrow 2 H_{2} + O_{2}'


@requires(parsing_library)
def test_ReactioN__unicode():
    keys = u'H2O H2 O2'.split()
    subst = {k: Substance.from_formula(k) for k in keys}
    r2 = Reaction.from_string("2 H2O -> 2 H2 + O2", subst)
    assert r2.unicode(subst) == u'2 H₂O → 2 H₂ + O₂'
    assert r2.html(subst) == \
        u'2 H<sub>2</sub>O &rarr; 2 H<sub>2</sub> + O<sub>2</sub>'


def test_Reaction__idempotency():
    with pytest.raises(ValueError):
        Reaction({'A': 1}, {'A': 1})
    with pytest.raises(ValueError):
        Reaction({}, {})
    with pytest.raises(ValueError):
        Reaction({'A': 1}, {'B': 1}, inact_reac={'B': 1}, inact_prod={'A': 1})


@requires('sympy')
def test_Equilibrium__eliminate():
    e1 = Equilibrium({'A': 1, 'B': 2}, {'C': 3})
    e2 = Equilibrium({'D': 5, 'B': 7}, {'E': 11})
    coeff = Equilibrium.eliminate([e1, e2], 'B')
    assert coeff == [7, -2]

    e3 = coeff[0]*e1 + coeff[1]*e2
    assert e3.net_stoich('B') == (0,)

    e4 = e1*coeff[0] + coeff[1]*e2
    assert e4.net_stoich('B') == (0,)

    assert (-e1).reac == {'C': 3}
    assert (e2*-3).reac == {'E': 33}


def test_Equilibrium__cancel():
    # 2B + C -> E
    e1 = Equilibrium({'A': 26, 'B': 20, 'C': 7}, {'D': 4, 'E': 7})
    e2 = Equilibrium({'A': 13, 'B': 3}, {'D': 2})
    coeff = e1.cancel(e2)
    assert coeff == -2


@requires('sympy')
def test_balance_stoichiometry():
    # 4 NH4ClO4 -> 2 N2 + 4 HCl + 6H2O + 5O2
    # 4 Al + 3O2 -> 2Al2O3
    # ---------------------------------------
    # 6 NH4ClO4 + 10 Al + -> 3 N2 + 6 HCl + 9 H2O + 5 Al2O3
    reac, prod = balance_stoichiometry({'NH4ClO4', 'Al'},
                                       {'Al2O3', 'HCl', 'H2O', 'N2'})
    assert reac == {'NH4ClO4': 6, 'Al': 10}
    assert prod == {'Al2O3': 5, 'HCl': 6, 'H2O': 9, 'N2': 3}

    with pytest.raises(ValueError):
        reac, prod = balance_stoichiometry({'C7H5(NO2)3', 'Al', 'NH4NO3'}, {'CO', 'H2O', 'N2'})

    r2, p2 = balance_stoichiometry({'Na2CO3'}, {'Na2O', 'CO2'})
    assert r2 == {'Na2CO3': 1}
    assert p2 == {'Na2O': 1, 'CO2': 1}

    r3, p3 = balance_stoichiometry({'C2H6', 'O2'}, {'H2O', 'CO2'})
    assert r3 == {'C2H6': 2, 'O2': 7}
    assert p3 == {'CO2': 4, 'H2O': 6}
    with pytest.raises(ValueError):
        reac, prod = balance_stoichiometry({'C2H6', 'O2'}, {'H2O', 'CO2', 'CO'})


@requires(parsing_library)
def test_ReactionSystem__from_string():
    rs = ReactionSystem.from_string('-> H + OH; Radiolytic(2.1e-7)', checks=())
    assert rs.rxns[0].reac == {}
    assert rs.rxns[0].prod == {'H': 1, 'OH': 1}
    assert rs.rxns[0].param.args == [2.1e-7]
    ref = 2.1e-7 * 0.15 * 998
    assert rs.rates({'doserate': .15, 'density': 998}) == {'H': ref, 'OH': ref}
