# -*- coding: utf-8 -*-

import pytest
from pyparsing import ParseException

from ..parsing import (
    formula_to_composition,
    formula_to_html,
    formula_to_latex,
    formula_to_unicode,
    parsing_library,
    to_reaction,
)

from ..testing import requires


@requires(parsing_library)
@pytest.mark.parametrize(
    "species, composition",
    [
        ("H2O*", {1: 2, 8: 1}),
        ("H2O'", {1: 2, 8: 1}),
        ("H2O''", {1: 2, 8: 1}),
        ("H2O'''", {1: 2, 8: 1}),
        ("Na", {11: 1}),
        ("Na*", {11: 1}),
        ("Na'", {11: 1}),
        ("Na''", {11: 1}),
        ("Na'''", {11: 1}),
        ("Na(g)", {11: 1}),
        ("Na*(g)", {11: 1}),
        ("Na'(g)", {11: 1}),
        ("Na''(g)", {11: 1}),
        ("Na'''(g)", {11: 1}),
        ("Na+", {0: 1, 11: 1}),
        ("Na*+", {0: 1, 11: 1}),
        ("Na'+", {0: 1, 11: 1}),
        ("Na''+", {0: 1, 11: 1}),
        ("Na'''+", {0: 1, 11: 1}),
        ("Na+(g)", {0: 1, 11: 1}),
        ("Na*+(g)", {0: 1, 11: 1}),
        ("Na'+(g)", {0: 1, 11: 1}),
        ("Na''+(g)", {0: 1, 11: 1}),
        ("Na'''+(g)", {0: 1, 11: 1}),
        ("Na2CO3*", {11: 2, 6: 1, 8: 3}),
        ("Na2CO3..7H2O*(s)", {11: 2, 6: 1, 8: 10, 1: 14}),
    ],
)
def test_formula_to_composition_primes(species, composition):
    """Should parse special species."""
    assert formula_to_composition(species) == composition


@pytest.mark.parametrize(
    "species, composition",
    [
        ("CO2(g)", {6: 1, 8: 2}),
        ("CO2(l)", {6: 1, 8: 2}),
        ("CO2(s)", {6: 1, 8: 2}),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_state_in_suffixes(species, composition):
    """Should parse species with state in suffixes."""
    assert (
        formula_to_composition(
            species,
            suffixes=("(g)", "(l)", "(s)"),
        )
        == composition
    )


@pytest.mark.parametrize(
    "species, composition",
    [
        ("CO2(aq)", {6: 1, 8: 2}),
        ("H2O(aq)", {1: 2, 8: 1}),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_state_not_in_suffixes(species, composition):
    """Should parse species without state in suffixes."""
    assert (
        formula_to_composition(
            species,
            suffixes=("(g)", "(l)", "(s)"),
        )
        == composition
    )


@pytest.mark.parametrize(
    "species, composition",
    [
        ("Li@C60", {3: 1, 6: 60}),
        ("Li@C60Cl", {3: 1, 6: 60, 17: 1}),
        ("(Li@C60)+", {0: 1, 3: 1, 6: 60}),
        ("Na@C60", {11: 1, 6: 60}),
        ("(Na@C60)+", {0: 1, 11: 1, 6: 60}),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_caged(species, composition):
    """Should parse cage species."""
    assert formula_to_composition(species) == composition


@pytest.mark.parametrize(
    "species",
    [
        ("ch3oh"),
        ("Ch3oh"),
        ("Ch3OH"),
        ("ch3OH(l)"),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_fail(species):
    """Should raise an exception."""
    with pytest.raises(ParseException):
        formula_to_composition(species)


@pytest.mark.parametrize(
    "species, composition",
    [
        ("Cl-", {0: -1, 17: 1}),
        ("Fe(SCN)2+", {0: 1, 6: 2, 7: 2, 16: 2, 26: 1}),
        ("Fe(SCN)2+1", {0: 1, 6: 2, 7: 2, 16: 2, 26: 1}),
        ("Fe+3", {0: 3, 26: 1}),
        ("NH4+", {0: 1, 1: 4, 7: 1}),
        ("Na+", {0: 1, 11: 1}),
        ("Na+1", {0: 1, 11: 1}),
        ("OH-", {0: -1, 1: 1, 8: 1}),
        ("SO4-2(aq)", {0: -2, 8: 4, 16: 1}),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_ions(species, composition):
    assert formula_to_composition(species) == composition


@pytest.mark.parametrize(
    "species",
    [
        # Ions.
        ("Cl/-"),
        ("Cl/-(aq)"),
        ("Fe(SCN)2/+"),
        ("Fe(SCN)2/+(aq)"),
        ("Fe/3+"),
        ("Fe/3+(aq)"),
        ("Na/+"),
        ("Na/+(aq)"),
        # Electrons.
        ("e/-"),
        ("e/-(aq)"),
        # Radicals.
        (".NO3/2-"),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_deprecated_charge(species):
    with pytest.raises(ValueError):
        formula_to_composition(species)


@pytest.mark.parametrize(
    "species",
    [
        ("Na+Cl-"),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_bad_charges(species):
    with pytest.raises(ValueError):
        formula_to_composition(species)


@pytest.mark.parametrize(
    "species, composition",
    [
        # With and without water of hydration.
        ("BaCl2", {17: 2, 56: 1}),
        ("BaCl2(s)", {17: 2, 56: 1}),
        ("BaCl2..2H2O(s)", {1: 4, 8: 2, 17: 2, 56: 1}),
        ("Na2CO3..7H2O(s)", {11: 2, 6: 1, 8: 10, 1: 14}),
        ("NaCl", {11: 1, 17: 1}),
        ("NaCl(s)", {11: 1, 17: 1}),
        ("Ni", {28: 1}),
        ("NI", {7: 1, 53: 1}),
        ("KF", {9: 1, 19: 1}),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_ionic_compounds(species, composition):
    assert formula_to_composition(species) == composition


@pytest.mark.parametrize(
    "species, composition",
    [
        # With and without water of hydration.
        ("Al2(SO4)3", {8: 12, 13: 2, 16: 3}),
        ("Al2(SO4)3(s)", {8: 12, 13: 2, 16: 3}),
        ("Al2(SO4)3(aq)", {8: 12, 13: 2, 16: 3}),
        ("K4[Fe(CN)6]", {6: 6, 7: 6, 19: 4, 26: 1}),
        ("K4[Fe(CN)6](s)", {6: 6, 7: 6, 19: 4, 26: 1}),
        ("K4[Fe(CN)6](aq)", {6: 6, 7: 6, 19: 4, 26: 1}),
        (
            "[Fe(H2O)6][Fe(CN)6]..19H2O",
            {
                1: 50,
                6: 6,
                7: 6,
                8: 25,
                26: 2,
            },
        ),
        (
            "[Fe(H2O)6][Fe(CN)6]..19H2O(s)",
            {
                1: 50,
                6: 6,
                7: 6,
                8: 25,
                26: 2,
            },
        ),
        (
            "[Fe(H2O)6][Fe(CN)6]..19H2O(aq)",
            {
                1: 50,
                6: 6,
                7: 6,
                8: 25,
                26: 2,
            },
        ),
        (
            "[Fe(CN)6]-3",
            {
                0: -3,
                6: 6,
                7: 6,
                26: 1,
            },
        ),
        (
            "[Fe(CN)6]-3(aq)",
            {
                0: -3,
                6: 6,
                7: 6,
                26: 1,
            },
        ),
        (
            "Ag[NH3]+",
            {
                0: 1,
                1: 3,
                7: 1,
                47: 1,
            },
        ),
        (
            "[Ni(NH3)6]+2",
            {
                0: 2,
                1: 18,
                7: 6,
                28: 1,
            },
        ),
        (
            "[PtCl6]-2",
            {
                0: -2,
                17: 6,
                78: 1,
            },
        ),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_complexes(species, composition):
    assert formula_to_composition(species) == composition


@pytest.mark.parametrize(
    "species",
    [
        ("[Fe(CN)6)-3"),
        ("(Fe(CN)6]-3"),
        ("[Fe(CN]6]-3"),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_bad_complexes(species):
    with pytest.raises(ParseException):
        formula_to_composition(species)


@pytest.mark.parametrize(
    "species, composition",
    [
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6",
            {
                6: 6,
                8: 18,
                12: 5.395,
                20: 2.832,
                26: 0.6285,
            },
        ),
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6(s)",
            {
                6: 6,
                8: 18,
                12: 5.395,
                20: 2.832,
                26: 0.6285,
            },
        ),
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6..8H2O(s)",
            {
                1: 16,
                6: 6,
                8: 26,
                12: 5.395,
                20: 2.832,
                26: 0.6285,
            },
        ),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_fractional_subscripts(species, composition):
    assert formula_to_composition(species) == composition


@pytest.mark.parametrize(
    "species, composition",
    [
        ("e-", {0: -1}),
        ("e-1", {0: -1}),
        ("e-(aq)", {0: -1}),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_solvated_electrons(species, composition):
    assert formula_to_composition(species) == composition


@pytest.mark.parametrize(
    "species, composition",
    [
        ("H2O", {1: 2, 8: 1}),
        ("((H2O)2OH)12", {1: 60, 8: 36}),
        ("PCl5", {15: 1, 17: 5}),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_covalent_compounds(species, composition):
    assert formula_to_composition(species) == composition


@pytest.mark.parametrize(
    "species, composition",
    [
        ("CH4(g)", {1: 4, 6: 1}),
        ("CH3CH3(g)", {1: 6, 6: 2}),
        # Many ways to write benzene.
        ("C6H6(l)", {1: 6, 6: 6}),
        ("(CH)6(l)", {1: 6, 6: 6}),
        ("CHCHCHCHCHCH(l)", {1: 6, 6: 6}),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_organic_compounds(species, composition):
    assert formula_to_composition(species) == composition


@pytest.mark.parametrize(
    "species, composition",
    [
        (".NO2(g)", {7: 1, 8: 2}),
        (".NH2", {1: 2, 7: 1}),
        ("ONOOH", {1: 1, 7: 1, 8: 3}),
        (".ONOO", {7: 1, 8: 3}),
        (".NO3-2", {0: -2, 7: 1, 8: 3}),
        # Deprecated charges as of 0.8.0.
        # (".NO3/2-", {0: -2, 7: 1, 8: 3}),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_radicals(species, composition):
    assert formula_to_composition(species) == composition


@requires(parsing_library)
def test_formula_to_composition_structural_formulas():
    with pytest.raises(ValueError):
        formula_to_composition("F-F")


@pytest.mark.parametrize(
    "species, composition",
    [
        ("alpha-FeOOH(s)", {1: 1, 8: 2, 26: 1}),
        ("epsilon-Zn(OH)2(s)", {1: 2, 8: 2, 30: 1}),
    ],
)
@requires(parsing_library)
def test_formula_to_composition_crystal_phases(species, composition):
    assert formula_to_composition(species) == composition


@requires(parsing_library)
def test_to_reaction():
    from chempy.chemistry import Reaction, Equilibrium

    rxn = to_reaction(
        "H+ + OH- -> H2O; 1.4e11; ref={'doi': '10.1039/FT9908601539'}",
        "H+ OH- H2O".split(),
        "->",
        Reaction,
    )
    assert rxn.__class__ == Reaction

    assert rxn.reac["H+"] == 1
    assert rxn.reac["OH-"] == 1
    assert rxn.prod["H2O"] == 1
    assert rxn.param == 1.4e11
    assert rxn.ref["doi"].startswith("10.")

    eq = to_reaction(
        "H+ + OH- = H2O; 1e-14; ref='rt, [H2O] == 1 M'",
        "H+ OH- H2O".split(),
        "=",
        Equilibrium,
    )
    assert eq.__class__ == Equilibrium

    assert eq.reac["H+"] == 1
    assert eq.reac["OH-"] == 1
    assert eq.prod["H2O"] == 1
    assert eq.ref.startswith("rt")

    for s in [
        "2 e-(aq) + (2 H2O) -> H2 + 2 OH- ; 1e6 ; ",
        "2 * e-(aq) + (2 H2O) -> 1 * H2 + 2 * OH- ; 1e6 ; ",
    ]:
        rxn2 = to_reaction(s, "e-(aq) H2 OH- H2O".split(), "->", Reaction)
        assert rxn2.__class__ == Reaction
        assert rxn2.reac["e-(aq)"] == 2
        assert rxn2.inact_reac["H2O"] == 2
        assert rxn2.prod["H2"] == 1
        assert rxn2.prod["OH-"] == 2
        assert rxn2.param == 1e6

    r1 = to_reaction("-> H2O", None, "->", Reaction)
    assert r1.reac == {}
    assert r1.prod == {"H2O": 1}
    assert r1.param is None

    r2 = to_reaction("H2O ->", None, "->", Reaction)
    assert r2.reac == {"H2O": 1}
    assert r2.prod == {}
    assert r2.param is None

    from chempy.kinetics.rates import MassAction

    ma = MassAction([3.14])
    r3 = to_reaction("H+ + OH- -> H2O", None, "->", Reaction, param=ma)
    assert r3.param.args == [3.14]

    rxn3 = to_reaction(
        "H2O + H2O -> H3O+ + OH-", "H3O+ OH- H2O".split(), "->", Reaction
    )
    assert rxn3.reac == {"H2O": 2} and rxn3.prod == {"H3O+": 1, "OH-": 1}

    rxn4 = to_reaction(
        "2 e-(aq) + (2 H2O) + (2 H+) -> H2 + 2 H2O",
        "e-(aq) H2 H2O H+".split(),
        "->",
        Reaction,
    )
    assert (
        rxn4.reac == {"e-(aq)": 2}
        and rxn4.inact_reac == {"H2O": 2, "H+": 2}
        and rxn4.prod == {"H2": 1, "H2O": 2}
    )


@pytest.mark.parametrize(
    "species, latex",
    [
        ("H2O", "H_{2}O"),
        # ("C6H6/+", "C_{6}H_{6}^{+}"),
        ("C6H6+", "C_{6}H_{6}^{+}"),
        # ("C18H38/2+", "C_{18}H_{38}^{2+}"),
        # ("C18H38/+2", "C_{18}H_{38}^{2+}"),
        ("C18H38+2", "C_{18}H_{38}^{2+}"),
        ("NaCl", "NaCl"),
        ("NaCl(s)", "NaCl(s)"),
        ("e-(aq)", "e^{-}(aq)"),
        ("Ca+2(aq)", "Ca^{2+}(aq)"),
        (".NO2(g)", r"^\bullet NO_{2}(g)"),
        (".NH2", r"^\bullet NH_{2}"),
        ("ONOOH", "ONOOH"),
        (".ONOO", r"^\bullet ONOO"),
        # (".NO3/2-", r"^\bullet NO_{3}^{2-}"),
        (".NO3-2", r"^\bullet NO_{3}^{2-}"),
        ("alpha-FeOOH(s)", r"\alpha-FeOOH(s)"),
        ("epsilon-Zn(OH)2(s)", r"\varepsilon-Zn(OH)_{2}(s)"),
        ("Na2CO3..7H2O(s)", r"Na_{2}CO_{3}\cdot 7H_{2}O(s)"),
        ("Na2CO3..1H2O(s)", r"Na_{2}CO_{3}\cdot H_{2}O(s)"),
        ("K4[Fe(CN)6]", r"K_{4}[Fe(CN)_{6}]"),
        ("K4[Fe(CN)6](s)", r"K_{4}[Fe(CN)_{6}](s)"),
        ("K4[Fe(CN)6](aq)", r"K_{4}[Fe(CN)_{6}](aq)"),
        ("[Fe(H2O)6][Fe(CN)6]..19H2O", r"[Fe(H_{2}O)_{6}][Fe(CN)_{6}]\cdot 19H_{2}O"),
        (
            "[Fe(H2O)6][Fe(CN)6]..19H2O(s)",
            r"[Fe(H_{2}O)_{6}][Fe(CN)_{6}]\cdot 19H_{2}O(s)",
        ),
        (
            "[Fe(H2O)6][Fe(CN)6]..19H2O(aq)",
            r"[Fe(H_{2}O)_{6}][Fe(CN)_{6}]\cdot 19H_{2}O(aq)",
        ),
        ("[Fe(CN)6]-3", r"[Fe(CN)_{6}]^{3-}"),
        ("[Fe(CN)6]-3(aq)", r"[Fe(CN)_{6}]^{3-}(aq)"),
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6",
            r"Ca_{2.832}Fe_{0.6285}Mg_{5.395}(CO_{3})_{6}",
        ),
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6(s)",
            r"Ca_{2.832}Fe_{0.6285}Mg_{5.395}(CO_{3})_{6}(s)",
        ),
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6..8H2O(s)",
            r"Ca_{2.832}Fe_{0.6285}Mg_{5.395}(CO_{3})_{6}\cdot 8H_{2}O(s)",
        ),
    ],
)
@requires(parsing_library)
def test_formula_to_latex(species, latex):
    assert formula_to_latex(species) == latex


@pytest.mark.parametrize(
    "species, latex",
    [
        # Parentheses.
        ("Fe(CN)6-3", "Fe(CN)_{6}^{3-}"),
        ("((H2O)2OH)12", "((H_{2}O)_{2}OH)_{12}"),
        # Square brackets.
        ("Fe[CN]6-3", "Fe[CN]_{6}^{3-}"),
        ("[(H2O)2OH]12", "[(H_{2}O)_{2}OH]_{12}"),
        # Curly braces.
        ("Fe{CN}6-3", r"Fe\{CN\}_{6}^{3-}"),
        ("{(H2O)2OH}12", r"\{(H_{2}O)_{2}OH\}_{12}"),
    ],
)
@requires(parsing_library)
def test_formula_to_latex_braces(species, latex):
    assert formula_to_latex(species) == latex


@pytest.mark.parametrize(
    "species, latex",
    [
        ("Li@C60", r"Li@C_{60}"),
        ("(Li@C60)+", r"(Li@C_{60})^{+}"),
        ("Na@C60", r"Na@C_{60}"),
        ("(Na@C60)+", r"(Na@C_{60})^{+}"),
    ],
)
@requires(parsing_library)
def test_formula_to_latex_caged(species, latex):
    """Should produce LaTeX for cage species."""
    assert formula_to_latex(species) == latex


@pytest.mark.parametrize(
    "species, unicode",
    [
        ("NH4+", u"NH₄⁺"),
        ("H2O", u"H₂O"),
        # ("C6H6/+", u"C₆H₆⁺"),
        ("C6H6+", u"C₆H₆⁺"),
        # ("Fe(CN)6/3-", u"Fe(CN)₆³⁻"),
        ("Fe(CN)6-3", u"Fe(CN)₆³⁻"),
        # ("C18H38/2+", u"C₁₈H₃₈²⁺"),
        # ("C18H38/+2", u"C₁₈H₃₈²⁺"),
        ("C18H38+2", u"C₁₈H₃₈²⁺"),
        ("((H2O)2OH)12", u"((H₂O)₂OH)₁₂"),
        ("[(H2O)2OH]12", u"[(H₂O)₂OH]₁₂"),
        ("{(H2O)2OH}12", u"{(H₂O)₂OH}₁₂"),
        ("NaCl", u"NaCl"),
        ("NaCl(s)", u"NaCl(s)"),
        ("e-(aq)", u"e⁻(aq)"),
        ("Ca+2(aq)", u"Ca²⁺(aq)"),
        (".NO2(g)", u"⋅NO₂(g)"),
        (".NH2", u"⋅NH₂"),
        ("ONOOH", u"ONOOH"),
        (".ONOO", u"⋅ONOO"),
        # (".NO3/2-", u"⋅NO₃²⁻"),
        (".NO3-2", u"⋅NO₃²⁻"),
        ("alpha-FeOOH(s)", u"α-FeOOH(s)"),
        ("epsilon-Zn(OH)2(s)", u"ε-Zn(OH)₂(s)"),
        ("Na2CO3..7H2O(s)", u"Na₂CO₃·7H₂O(s)"),
        ("Na2CO3..1H2O(s)", u"Na₂CO₃·H₂O(s)"),
        ("K4[Fe(CN)6]", r"K₄[Fe(CN)₆]"),
        ("K4[Fe(CN)6](s)", r"K₄[Fe(CN)₆](s)"),
        ("K4[Fe(CN)6](aq)", r"K₄[Fe(CN)₆](aq)"),
        ("[Fe(H2O)6][Fe(CN)6]..19H2O", r"[Fe(H₂O)₆][Fe(CN)₆]·19H₂O"),
        ("[Fe(H2O)6][Fe(CN)6]..19H2O(s)", r"[Fe(H₂O)₆][Fe(CN)₆]·19H₂O(s)"),
        ("[Fe(H2O)6][Fe(CN)6]..19H2O(aq)", r"[Fe(H₂O)₆][Fe(CN)₆]·19H₂O(aq)"),
        ("[Fe(CN)6]-3", r"[Fe(CN)₆]³⁻"),
        ("[Fe(CN)6]-3(aq)", r"[Fe(CN)₆]³⁻(aq)"),
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6",
            r"Ca₂.₈₃₂Fe₀.₆₂₈₅Mg₅.₃₉₅(CO₃)₆",
        ),
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6(s)",
            r"Ca₂.₈₃₂Fe₀.₆₂₈₅Mg₅.₃₉₅(CO₃)₆(s)",
        ),
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6..8H2O(s)",
            r"Ca₂.₈₃₂Fe₀.₆₂₈₅Mg₅.₃₉₅(CO₃)₆·8H₂O(s)",
        ),
        (
            "Zn(NO3)2..6H2O",
            r"Zn(NO₃)₂·6H₂O",
        ),
    ],
)
@requires(parsing_library)
def test_formula_to_unicode(species, unicode):
    assert formula_to_unicode(species) == unicode


@pytest.mark.parametrize(
    "species, unicode",
    [
        ("Li@C60", r"Li@C₆₀"),
        ("(Li@C60)+", r"(Li@C₆₀)⁺"),
        ("Na@C60", r"Na@C₆₀"),
        ("(Na@C60)+", r"(Na@C₆₀)⁺"),
    ],
)
@requires(parsing_library)
def test_formula_to_unicode_caged(species, unicode):
    """Should produce LaTeX for cage species."""
    assert formula_to_unicode(species) == unicode


@pytest.mark.parametrize(
    "species, html",
    [
        ("H2O", "H<sub>2</sub>O"),
        # ("C6H6/+", "C<sub>6</sub>H<sub>6</sub><sup>+</sup>"),
        ("C6H6+", "C<sub>6</sub>H<sub>6</sub><sup>+</sup>"),
        # ("Fe(CN)6/3-", "Fe(CN)<sub>6</sub><sup>3-</sup>"),
        ("Fe(CN)6-3", "Fe(CN)<sub>6</sub><sup>3-</sup>"),
        # ("C18H38/2+", "C<sub>18</sub>H<sub>38</sub><sup>2+</sup>"),
        # ("C18H38/+2", "C<sub>18</sub>H<sub>38</sub><sup>2+</sup>"),
        ("C18H38+2", "C<sub>18</sub>H<sub>38</sub><sup>2+</sup>"),
        ("((H2O)2OH)12", "((H<sub>2</sub>O)<sub>2</sub>OH)<sub>12</sub>"),
        ("[(H2O)2OH]12", "[(H<sub>2</sub>O)<sub>2</sub>OH]<sub>12</sub>"),
        ("{(H2O)2OH}12", "{(H<sub>2</sub>O)<sub>2</sub>OH}<sub>12</sub>"),
        ("NaCl", "NaCl"),
        ("NaCl(s)", "NaCl(s)"),
        ("e-(aq)", "e<sup>-</sup>(aq)"),
        ("Ca+2(aq)", "Ca<sup>2+</sup>(aq)"),
        (".NO2(g)", r"&sdot;NO<sub>2</sub>(g)"),
        (".NH2", r"&sdot;NH<sub>2</sub>"),
        ("ONOOH", "ONOOH"),
        (".ONOO", r"&sdot;ONOO"),
        # (".NO3/2-", r"&sdot;NO<sub>3</sub><sup>2-</sup>"),
        (".NO3-2", r"&sdot;NO<sub>3</sub><sup>2-</sup>"),
        ("alpha-FeOOH(s)", r"&alpha;-FeOOH(s)"),
        ("epsilon-Zn(OH)2(s)", (r"&epsilon;-Zn(OH)<sub>2</sub>(s)")),
        ("Na2CO3..7H2O(s)", "Na<sub>2</sub>CO<sub>3</sub>&sdot;7H<sub>2</sub>O(s)"),
        ("Na2CO3..1H2O(s)", "Na<sub>2</sub>CO<sub>3</sub>&sdot;H<sub>2</sub>O(s)"),
        ("K4[Fe(CN)6]", r"K<sub>4</sub>[Fe(CN)<sub>6</sub>]"),
        ("K4[Fe(CN)6](s)", r"K<sub>4</sub>[Fe(CN)<sub>6</sub>](s)"),
        ("K4[Fe(CN)6](aq)", r"K<sub>4</sub>[Fe(CN)<sub>6</sub>](aq)"),
        (
            "[Fe(H2O)6][Fe(CN)6]..19H2O",
            r"[Fe(H<sub>2</sub>O)<sub>6</sub>][Fe(CN)<sub>6</sub>]&sdot;19H<sub>2</sub>O",
        ),
        (
            "[Fe(H2O)6][Fe(CN)6]..19H2O(s)",
            r"[Fe(H<sub>2</sub>O)<sub>6</sub>][Fe(CN)<sub>6</sub>]&sdot;19H<sub>2</sub>O(s)",
        ),
        (
            "[Fe(H2O)6][Fe(CN)6]..19H2O(aq)",
            r"[Fe(H<sub>2</sub>O)<sub>6</sub>][Fe(CN)<sub>6</sub>]&sdot;19H<sub>2</sub>O(aq)",
        ),
        ("[Fe(CN)6]-3", r"[Fe(CN)<sub>6</sub>]<sup>3-</sup>"),
        ("[Fe(CN)6]-3(aq)", r"[Fe(CN)<sub>6</sub>]<sup>3-</sup>(aq)"),
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6",
            r"Ca<sub>2.832</sub>Fe<sub>0.6285</sub>Mg<sub>5.395</sub>(CO<sub>3</sub>)<sub>6</sub>",
        ),
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6(s)",
            r"Ca<sub>2.832</sub>Fe<sub>0.6285</sub>Mg<sub>5.395</sub>(CO<sub>3</sub>)<sub>6</sub>(s)",
        ),
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6..8H2O(s)",
            r"Ca<sub>2.832</sub>Fe<sub>0.6285</sub>Mg<sub>5.395</sub>(CO<sub>3</sub>)<sub>6</sub>&sdot;8H<sub>2</sub>O(s)",
        ),
        (
            "Ca2.832Fe0.6285Mg5.395(CO3)6..8H2O(s)",
            r"Ca<sub>2.832</sub>Fe<sub>0.6285</sub>Mg<sub>5.395</sub>(CO<sub>3</sub>)<sub>6</sub>&sdot;8H<sub>2</sub>O(s)",
        ),
        (
            "Zn(NO3)2..6H2O",
            r"Zn(NO<sub>3</sub>)<sub>2</sub>&sdot;6H<sub>2</sub>O",
        ),
    ],
)
@requires(parsing_library)
def test_formula_to_html(species, html):
    assert formula_to_html(species) == html


@pytest.mark.parametrize(
    "species, html",
    [
        ("Li@C60", r"Li@C<sub>60</sub>"),
        ("(Li@C60)+", r"(Li@C<sub>60</sub>)<sup>+</sup>"),
        ("Na@C60", r"Na@C<sub>60</sub>"),
        ("(Na@C60)+", r"(Na@C<sub>60</sub>)<sup>+</sup>"),
    ],
)
@requires(parsing_library)
def test_formula_to_html_caged(species, html):
    """Should produce HTML for cage species."""
    assert formula_to_html(species) == html


def test_composition_dot_as_crystal_water_chempy08x():
    ref = {30: 1, 7: 2, 8: 12, 1: 12}
    assert formula_to_composition('Zn(NO3)2{}6H2O'.format('\u00B7')) == ref
    assert formula_to_composition('Zn(NO3)2..6H2O') == ref
