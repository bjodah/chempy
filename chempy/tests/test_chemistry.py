# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from functools import reduce
from operator import attrgetter, add

import pytest

from ..util.arithmeticdict import ArithmeticDict
from ..util.testing import requires
from ..util.parsing import parsing_library
from ..units import default_units, units_library, to_unitless
from ..chemistry import (
    equilibrium_quotient, Substance, Species, Reaction,
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
    H2O = Substance(name='H2O',  charge=0, latex_name=r'\mathrm{H_{2}O}',
                    data={'pKa': 14})  # will_be_missing_in='0.8.0', use data=...
    OH_m = Substance(name='OH-',  charge=-1, latex_name=r'\mathrm{OH^{-}}')
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

    uranyl_ads = Species.from_formula('UO2+2(ads)', phases={'(aq)': 0, '(ads)': 1})
    assert uranyl_ads.composition == {0: 2, 92: 1, 8: 2}
    assert uranyl_ads.phase_idx == 1


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
    assert r1.composition_violation(substance_dict, ['H+']) == [0]
    viol, cks = r1.composition_violation(substance_dict, True)
    assert viol == [0]*3 and sorted(cks) == [0, 1, 8]
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

    lhs5, rhs5 = {'H+': 1, 'OH-': 1}, {'H2O': 1}
    r5 = Reaction(lhs5, rhs5)
    assert r5.reac == lhs5 and r5.prod == rhs5


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

    with pytest.raises(ValueError):
        Reaction({'H2O': 2}, {'H2O2': 2, 'O2': -2})

    r4 = Reaction({'H+': 2, 'OH-': 1}, {'H2O': 2}, 42.0)
    assert Reaction.from_string(str(r4), 'H+ OH- H2O') == r4
    assert Reaction.from_string(str(r4), None) == r4


@requires(parsing_library, units_library)
def test_Reaction_from_string__units():
    r5 = Reaction.from_string('2 H2O2 -> O2 + 2 H2O; 1e-7/molar/second', 'H2O O2 H2O2')
    assert to_unitless(r5.param, 1/default_units.molar/default_units.second) == 1e-7
    r6 = Reaction.from_string('->', checks=())
    assert r6.reac == {} and r6.prod == {}

    r7 = Reaction.from_string('2 A -> B; 2e-3*metre**3/mol/hour', None)
    assert r7.reac == {'A': 2} and r7.prod == {'B': 1}
    assert r7.param == 2e-3*default_units.metre**3/default_units.mol/default_units.hour

    with pytest.raises(ValueError):
        Reaction.from_string('2 A -> B; 2e-3/hour', None)

    r8 = Reaction.from_string('A -> B; "k"')
    assert r8.rate_expr().args is None
    assert r8.rate_expr().unique_keys == ('k',)
    r9 = Reaction.from_string('A -> B; 42.0')
    assert r9.rate_expr().args == [42.0]
    assert r9.rate_expr().unique_keys is None

    Reaction.from_string("H+ + OH- -> H2O; 1e10/M/s", 'H2O H+ OH-'.split())
    with pytest.raises(ValueError):
        Reaction.from_string("H2O -> H+ + OH-; 1e-4/M/s", 'H2O H+ OH-'.split())


@requires(parsing_library, units_library)
def test_Substance__molar_mass():
    mw_water = Substance.from_formula('H2O').molar_mass(default_units)
    q = mw_water / ((15.9994 + 2*1.008)*default_units.gram/default_units.mol)
    assert abs(q - 1) < 1e-3


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
def test_ReactioN__latex():
    keys = 'H2O H2 O2'.split()
    subst = {k: Substance.from_formula(k) for k in keys}
    r2 = Reaction.from_string("2 H2O -> 2 H2 + O2", subst)
    assert r2.latex(subst) == r'2 H_{2}O \rightarrow 2 H_{2} + O_{2}'


@requires(parsing_library)
def test_Reaction__unicode():
    keys = u'H2O H2 O2'.split()
    subst = {k: Substance.from_formula(k) for k in keys}
    r2 = Reaction.from_string("2 H2O -> 2 H2 + O2", subst)
    assert r2.unicode(subst) == u'2 H₂O → 2 H₂ + O₂'


@requires(parsing_library)
def test_Reaction__html():
    keys = 'H2O H2 O2'.split()
    subst = {k: Substance.from_formula(k) for k in keys}
    r2 = Reaction.from_string("2 H2O -> 2 H2 + O2", subst)
    assert r2.html(subst) == \
        '2 H<sub>2</sub>O &rarr; 2 H<sub>2</sub> + O<sub>2</sub>'
    assert r2.html(subst, Reaction_coeff_fmt=lambda s: '<b>{0}</b>'.format(s)) == \
        '<b>2</b> H<sub>2</sub>O &rarr; <b>2</b> H<sub>2</sub> + O<sub>2</sub>'
    assert r2.html(subst, Reaction_formula_fmt=lambda s: '<b>{0}</b>'.format(s)) == \
        '2 <b>H<sub>2</sub>O</b> &rarr; 2 <b>H<sub>2</sub></b> + <b>O<sub>2</sub></b>'


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


@requires(parsing_library, units_library)
def test_Equilibrium__from_string():
    Equilibrium.from_string('H2O = H+ + OH-')
    Equilibrium.from_string('H2O = H+ + OH-; 1e-14')
    Equilibrium.from_string('H2O = H+ + OH-; 1e-14*molar')
    with pytest.raises(ValueError):
        Equilibrium.from_string('H+ + OH- = H2O; 1e-14*molar')


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

    r3, p3 = balance_stoichiometry({'C2H6', 'O2'}, {'H2O', 'CO2'})
    assert r3 == {'C2H6': 2, 'O2': 7}
    assert p3 == {'CO2': 4, 'H2O': 6}

    r4, p4 = balance_stoichiometry({'C7H5(NO2)3', 'NH4NO3'}, {'CO', 'H2O', 'N2'})
    assert r4 == {'C7H5(NO2)3': 2, 'NH4NO3': 7}
    assert p4 == {'CO': 14, 'H2O': 19, 'N2': 10}

    a5, b5 = {"C3H5NO", "CH4", "NH3", "H2O"}, {"C2H6", "CH4O", "CH5N", "CH3N"}
    formulas = list(set.union(a5, b5))
    substances = dict(zip(formulas, map(Substance.from_formula, formulas)))
    compositions = {k: ArithmeticDict(int, substances[k].composition) for k in formulas}
    r5, p5 = balance_stoichiometry(a5, b5)
    compo_reac = dict(reduce(add, [compositions[k]*v for k, v in r5.items()]))
    compo_prod = dict(reduce(add, [compositions[k]*v for k, v in p5.items()]))
    assert compo_reac == compo_prod

    a6, b6 = map(lambda x: set(x.split()), 'CuSCN KIO3 HCl;CuSO4 KCl HCN ICl H2O'.split(';'))
    r6, p6 = balance_stoichiometry(a6, b6)
    assert r6 == dict(CuSCN=4, KIO3=7, HCl=14)
    assert p6 == dict(CuSO4=4, KCl=7, HCN=4, ICl=7, H2O=5)


@requires('sympy')
def test_balance_stoichiometry__ordering():
    reac, prod = 'CuSCN KIO3 HCl'.split(), 'CuSO4 KCl HCN ICl H2O'.split()
    rxn = Reaction(*balance_stoichiometry(reac, prod))
    res = rxn.string()
    ref = '4 CuSCN + 7 KIO3 + 14 HCl -> 4 CuSO4 + 7 KCl + 4 HCN + 7 ICl + 5 H2O'
    assert res == ref


@requires('sympy')
def test_balance_stoichiometry__simple():
    r2, p2 = balance_stoichiometry({'Na2CO3'}, {'Na2O', 'CO2'})
    assert r2 == {'Na2CO3': 1}
    assert p2 == {'Na2O': 1, 'CO2': 1}


@requires('sympy')
@pytest.mark.parametrize('underdet', [False, None, True])
def test_balance_stoichiometry__impossible(underdet):
    from pulp.solvers import PulpSolverError
    with pytest.raises((ValueError, PulpSolverError)):
        r1, p1 = balance_stoichiometry({'CO'}, {'CO2'}, underdetermined=underdet)


@requires('sympy', 'pulp')
def test_balance_stoichiometry__underdetermined():
    from pulp.solvers import PulpSolverError

    with pytest.raises(ValueError):
        balance_stoichiometry({'C2H6', 'O2'}, {'H2O', 'CO2', 'CO'}, underdetermined=False)
    reac, prod = balance_stoichiometry({'C2H6', 'O2'}, {'H2O', 'CO2', 'CO'})

    r1 = {'C7H5O3-', 'O2', 'C21H27N7O14P2-2', 'H+'}
    p1 = {'C7H5O4-', 'C21H26N7O14P2-', 'H2O'}  # see https://github.com/bjodah/chempy/issues/67
    bal1 = balance_stoichiometry(r1, p1, underdetermined=None)
    assert bal1 == ({'C21H27N7O14P2-2': 1, 'H+': 1, 'C7H5O3-': 1, 'O2': 1},
                    {'C21H26N7O14P2-': 1, 'H2O': 1, 'C7H5O4-': 1})

    with pytest.raises(ValueError):
        balance_stoichiometry({'C3H4O3', 'H3PO4'}, {'C3H6O3'}, underdetermined=None)

    for underdet in [False, True, None]:
        with pytest.raises((ValueError, PulpSolverError)):
            balance_stoichiometry({'C3H6O3'}, {'C3H4O3'}, underdetermined=underdet)

    with pytest.raises(ValueError):  # https://github.com/bjodah/chempy/pull/86#issuecomment-375421609
        balance_stoichiometry({'C21H36N7O16P3S', 'C3H4O3'}, {'H2O', 'C5H8O3', 'C24H38N7O18P3S'})


@requires('sympy', 'pulp')
def test_balance_stoichiometry__very_underdetermined():
    r3 = set('O2 Fe Al Cr'.split())
    p3 = set('FeO Fe2O3 Fe3O4 Al2O3 Cr2O3 CrO3'.split())
    bal3 = balance_stoichiometry(r3, p3, underdetermined=None)
    ref3 = {'Fe': 7, 'Al': 2, 'Cr': 3, 'O2': 9}, {k: 2 if k == 'FeO' else 1 for k in p3}
    substances = {k: Substance.from_formula(k) for k in r3 | p3}
    assert all(viol == 0 for viol in Reaction(*ref3).composition_violation(substances))
    assert sum(bal3[0].values()) + sum(bal3[1].values()) <= sum(ref3[0].values()) + sum(ref3[1].values())
    assert bal3 == ref3


@requires('sympy', 'pulp')
def test_balance_stoichiometry__underdetermined__canoncial():
    # This tests for canoncial representation of the underdetermined system
    # where all coefficients are integer and >= 1. It is however of limited
    # practical use (and hence marked ``xfail``) since underdetermined systems
    # have infinite number of solutions. It should however be possible to rewrite
    # the logic so that such canoncial results are returned from balance_stoichiometry
    r2 = {'O2', 'O3', 'C', 'NO', 'N2O', 'NO2', 'N2O4'}
    p2 = {'CO', 'CO2', 'N2'}
    bal2 = balance_stoichiometry(r2, p2, underdetermined=None)
    ref2 = ({'O2': 1, 'O3': 1, 'C': 7, 'NO': 1, 'N2O': 1, 'NO2': 1, 'N2O4': 1},
            {'CO': 1, 'CO2': 6, 'N2': 3})
    substances = {k: Substance.from_formula(k) for k in r2 | p2}
    assert all(viol == 0 for viol in Reaction(*ref2).composition_violation(substances))
    assert sum(bal2[0].values()) + sum(bal2[1].values()) <= sum(ref2[0].values()) + sum(ref2[1].values())
    assert bal2 == ref2


@requires('sympy', 'pulp')
def test_balance_stoichiometry__substances__underdetermined():
    substances = {s.name: s for s in [
        Substance('eggs_6pack', composition=dict(eggs=6)),
        Substance('milk_carton', composition=dict(cups_of_milk=4)),
        Substance('flour_bag', composition=dict(spoons_of_flour=30)),
        Substance('pancake', composition=dict(eggs=1, cups_of_milk=1, spoons_of_flour=2)),
        Substance('waffle', composition=dict(eggs=2, cups_of_milk=2, spoons_of_flour=3)),
    ]}
    ur1 = {'eggs_6pack', 'milk_carton', 'flour_bag'}
    up1 = {'pancake', 'waffle'}
    br1, bp1 = balance_stoichiometry(ur1, up1, substances=substances, underdetermined=None)
    ref_r1 = {'eggs_6pack': 6, 'flour_bag': 2, 'milk_carton': 9}
    ref_p1 = {'pancake': 12, 'waffle': 12}
    assert all(viol == 0 for viol in Reaction(ref_r1, ref_p1).composition_violation(substances))
    assert all(v > 0 for v in br1.values()) and all(v > 0 for v in bp1.values())
    assert bp1 == ref_p1
    assert br1 == ref_r1


@requires('sympy')
def test_balance_stoichiometry__missing_product_atom():
    with pytest.raises(ValueError):  # No Al on product side
        balance_stoichiometry({'C7H5(NO2)3', 'Al', 'NH4NO3'}, {'CO', 'H2O', 'N2'})
