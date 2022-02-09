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
from ..parsing import _get_formula_parser
from ..testing import requires


@requires(parsing_library)
def test_formula_to_composition_primes():
    """Should parse special species."""
    assert formula_to_composition("H2O*") == {1: 2, 8: 1}
    assert formula_to_composition("H2O'") == {1: 2, 8: 1}
    assert formula_to_composition("H2O''") == {1: 2, 8: 1}
    assert formula_to_composition("H2O'''") == {1: 2, 8: 1}
    assert formula_to_composition("Na*") == {11: 1}
    assert formula_to_composition("Na'") == {11: 1}
    assert formula_to_composition("Na''") == {11: 1}
    assert formula_to_composition("Na*(g)") == {11: 1}
    assert formula_to_composition("Na*+") == {0: 1, 11: 1}
    assert formula_to_composition("Na2CO3*") == {11: 2, 6: 1, 8: 3}
    assert formula_to_composition("Na2CO3..7H2O*(s)") == {11: 2, 6: 1, 8: 10, 1: 14}


@requires(parsing_library)
def test_formula_to_composition_state_in_suffixes():
    """Should parse species with state in suffixes."""
    assert formula_to_composition(
        "CO2(g)",
        suffixes=("(g)", "(l)", "(s)"),
    ) == {6: 1, 8: 2}
    assert formula_to_composition(
        "CO2(l)",
        suffixes=("(g)", "(l)", "(s)"),
    ) == {6: 1, 8: 2}
    assert formula_to_composition(
        "CO2(s)",
        suffixes=("(g)", "(l)", "(s)"),
    ) == {6: 1, 8: 2}


@requires(parsing_library)
def test_formula_to_composition_state_not_in_suffixes():
    """Should parse species without state in suffixes."""
    assert formula_to_composition(
        "CO2(aq)",
        suffixes=("(g)", "(l)", "(s)"),
    ) == {6: 1, 8: 2}


@requires(parsing_library)
def test_formula_to_composition_caged():
    """Should parse cage species."""
    assert formula_to_composition("Li@C60") == {3: 1, 6: 60}
    assert formula_to_composition("Li@C60Cl") == {3: 1, 6: 60, 17: 1}
    assert formula_to_composition("(Li@C60)+") == {0: 1, 3: 1, 6: 60}
    assert formula_to_composition("Na@C60") == {11: 1, 6: 60}
    assert formula_to_composition("(Na@C60)+") == {0: 1, 11: 1, 6: 60}


@requires(parsing_library)
def test_formula_to_composition_fail():
    """Should raise an exception."""
    with pytest.raises(ParseException):
        formula_to_composition("ch3oh")

    with pytest.raises(ParseException):
        print(_get_formula_parser().parseString("Ch3OH"))
        formula_to_composition("Ch3OH")

    with pytest.raises(ParseException):
        formula_to_composition("ch3oh")

    with pytest.raises(ParseException):
        formula_to_composition("Ch3OH(l)")


@requires(parsing_library)
def test_formula_to_composition_ions():
    assert formula_to_composition("Cl-") == {0: -1, 17: 1}
    assert formula_to_composition("Fe(SCN)2+") == {0: 1, 6: 2, 7: 2, 16: 2, 26: 1}
    assert formula_to_composition("Fe(SCN)2+1") == {0: 1, 6: 2, 7: 2, 16: 2, 26: 1}
    assert formula_to_composition("Fe+3") == {0: 3, 26: 1}
    assert formula_to_composition("NH4+") == {0: 1, 1: 4, 7: 1}
    assert formula_to_composition("Na+") == {0: 1, 11: 1}
    assert formula_to_composition("Na+1") == {0: 1, 11: 1}
    assert formula_to_composition("OH-") == {0: -1, 1: 1, 8: 1}
    assert formula_to_composition("SO4-2(aq)") == {0: -2, 8: 4, 16: 1}
    # Deprecated charges as of 0.8.0.
    # assert formula_to_composition("Fe/3+") == {0: 3, 26: 1}
    # assert formula_to_composition("Na/+") == {0: 1, 11: 1}
    # assert formula_to_composition("Cl/-") == {0: -1, 17: 1}
    # assert formula_to_composition("Fe(SCN)2/+") == {0: 1, 6: 2, 7: 2, 16: 2, 26: 1}


@requires(parsing_library)
def test_formula_to_composition_deprecated_charge():
    with pytest.raises(ValueError):
        # Ions.
        formula_to_composition("Cl/-")
        formula_to_composition("Cl/-(aq)")
        formula_to_composition("Fe(SCN)2/+")
        formula_to_composition("Fe(SCN)2/+(aq)")
        formula_to_composition("Fe/3+")
        formula_to_composition("Fe/3+(aq)")
        formula_to_composition("Na/+")
        formula_to_composition("Na/+(aq)")
        # Electrons.
        formula_to_composition("e/-")
        formula_to_composition("e/-(aq)")
        # Radicals.
        formula_to_composition(".NO3/2-")


@requires(parsing_library)
def test_formula_to_composition_bad_charges():
    with pytest.raises(ValueError):
        formula_to_composition("Na+Cl-")


@requires(parsing_library)
def test_formula_to_composition_ionic_compounds():
    # With and without water of hydration.
    assert formula_to_composition("BaCl2") == {17: 2, 56: 1}
    assert formula_to_composition("BaCl2(s)") == {17: 2, 56: 1}
    assert formula_to_composition("BaCl2..2H2O(s)") == {1: 4, 8: 2, 17: 2, 56: 1}
    assert formula_to_composition("Na2CO3..7H2O(s)") == {11: 2, 6: 1, 8: 10, 1: 14}
    assert formula_to_composition("NaCl") == {11: 1, 17: 1}
    assert formula_to_composition("NaCl(s)") == {11: 1, 17: 1}
    assert formula_to_composition("Ni") == {28: 1}
    assert formula_to_composition("NI") == {7: 1, 53: 1}
    assert formula_to_composition("KF") == {9: 1, 19: 1}


@requires(parsing_library)
def test_formula_to_composition_complexes():
    # With and without water of hydration.
    assert formula_to_composition("Al2(SO4)3") == {8: 12, 13: 2, 16: 3}
    assert formula_to_composition("Al2(SO4)3(s)") == {8: 12, 13: 2, 16: 3}
    assert formula_to_composition("Al2(SO4)3(aq)") == {8: 12, 13: 2, 16: 3}
    assert formula_to_composition("K4[Fe(CN)6]") == {6: 6, 7: 6, 19: 4, 26: 1}
    assert formula_to_composition("K4[Fe(CN)6](s)") == {6: 6, 7: 6, 19: 4, 26: 1}
    assert formula_to_composition("K4[Fe(CN)6](aq)") == {6: 6, 7: 6, 19: 4, 26: 1}
    assert formula_to_composition("[Fe(H2O)6][Fe(CN)6]..19H2O") == {
        1: 50,
        6: 6,
        7: 6,
        8: 25,
        26: 2,
    }
    assert formula_to_composition("[Fe(H2O)6][Fe(CN)6]..19H2O(s)") == {
        1: 50,
        6: 6,
        7: 6,
        8: 25,
        26: 2,
    }
    assert formula_to_composition("[Fe(H2O)6][Fe(CN)6]..19H2O(aq)") == {
        1: 50,
        6: 6,
        7: 6,
        8: 25,
        26: 2,
    }
    assert formula_to_composition("[Fe(CN)6]-3") == {
        0: -3,
        6: 6,
        7: 6,
        26: 1,
    }
    assert formula_to_composition("[Fe(CN)6]-3(aq)") == {
        0: -3,
        6: 6,
        7: 6,
        26: 1,
    }
    assert formula_to_composition("Ag[NH3]+") == {
        0: 1,
        1: 3,
        7: 1,
        47: 1,
    }
    assert formula_to_composition("[Ni(NH3)6]+2") == {
        0: 2,
        1: 18,
        7: 6,
        28: 1,
    }
    assert formula_to_composition("[PtCl6]-2") == {
        0: -2,
        17: 6,
        78: 1,
    }


@requires(parsing_library)
def test_formula_to_composition_bad_complexes():
    with pytest.raises(ParseException):
        formula_to_composition("[Fe(CN)6)-3")

    with pytest.raises(ParseException):
        formula_to_composition("(Fe(CN)6]-3")

    with pytest.raises(ParseException):
        formula_to_composition("[Fe(CN]6]-3")


@requires(parsing_library)
def test_formula_to_composition_fractional_subscripts():
    assert formula_to_composition("Ca2.832Fe0.6285Mg5.395(CO3)6") == {
        6: 6,
        8: 18,
        12: 5.395,
        20: 2.832,
        26: 0.6285,
    }
    assert formula_to_composition("Ca2.832Fe0.6285Mg5.395(CO3)6(s)") == {
        6: 6,
        8: 18,
        12: 5.395,
        20: 2.832,
        26: 0.6285,
    }
    assert formula_to_composition("Ca2.832Fe0.6285Mg5.395(CO3)6..8H2O(s)") == {
        1: 16,
        6: 6,
        8: 26,
        12: 5.395,
        20: 2.832,
        26: 0.6285,
    }


@requires(parsing_library)
def test_formula_to_composition_solvated_electrons():
    assert formula_to_composition("e-") == {0: -1}
    assert formula_to_composition("e-1") == {0: -1}
    assert formula_to_composition("e-(aq)") == {0: -1}
    # Deprecated charges as of 0.8.0.
    # assert formula_to_composition("e/-") == {0: -1}
    # assert formula_to_composition("e/-(aq)") == {0: -1}


@requires(parsing_library)
def test_formula_to_composition_covalent_compounds():
    assert formula_to_composition("H2O") == {1: 2, 8: 1}
    assert formula_to_composition("((H2O)2OH)12") == {1: 60, 8: 36}
    assert formula_to_composition("PCl5") == {15: 1, 17: 5}


@requires(parsing_library)
def test_formula_to_composition_organic_compounds():
    assert formula_to_composition("CH4(g)") == {1: 4, 6: 1}
    assert formula_to_composition("CH3CH3(g)") == {1: 6, 6: 2}
    # Many ways to write benzene.
    assert formula_to_composition("C6H6(l)") == {1: 6, 6: 6}
    assert formula_to_composition("(CH)6(l)") == {1: 6, 6: 6}
    assert formula_to_composition("CHCHCHCHCHCH(l)") == {1: 6, 6: 6}


@requires(parsing_library)
def test_formula_to_composition_radicals():
    assert formula_to_composition(".NO2(g)") == {7: 1, 8: 2}
    assert formula_to_composition(".NH2") == {1: 2, 7: 1}
    assert formula_to_composition("ONOOH") == {1: 1, 7: 1, 8: 3}
    assert formula_to_composition(".ONOO") == {7: 1, 8: 3}
    assert formula_to_composition(".NO3-2") == {0: -2, 7: 1, 8: 3}
    # Deprecated charges as of 0.8.0.
    # assert formula_to_composition(".NO3/2-") == {0: -2, 7: 1, 8: 3}


@requires(parsing_library)
def test_formula_to_composition_structural_formulas():
    with pytest.raises(ValueError):
        formula_to_composition("F-F")


@requires(parsing_library)
def test_formula_to_composition_crystal_phases():
    assert formula_to_composition("alpha-FeOOH(s)") == {1: 1, 8: 2, 26: 1}
    assert formula_to_composition("epsilon-Zn(OH)2(s)") == {1: 2, 8: 2, 30: 1}


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


@requires(parsing_library)
def test_formula_to_latex():
    assert formula_to_latex("H2O") == "H_{2}O"
    # assert formula_to_latex("C6H6/+") == "C_{6}H_{6}^{+}"
    assert formula_to_latex("C6H6+") == "C_{6}H_{6}^{+}"
    # assert formula_to_latex("C18H38/2+") == "C_{18}H_{38}^{2+}"
    # assert formula_to_latex("C18H38/+2") == "C_{18}H_{38}^{2+}"
    assert formula_to_latex("C18H38+2") == "C_{18}H_{38}^{2+}"
    assert formula_to_latex("NaCl") == "NaCl"
    assert formula_to_latex("NaCl(s)") == "NaCl(s)"
    assert formula_to_latex("e-(aq)") == "e^{-}(aq)"
    assert formula_to_latex("Ca+2(aq)") == "Ca^{2+}(aq)"
    assert formula_to_latex(".NO2(g)") == r"^\bullet NO_{2}(g)"
    assert formula_to_latex(".NH2") == r"^\bullet NH_{2}"
    assert formula_to_latex("ONOOH") == "ONOOH"
    assert formula_to_latex(".ONOO") == r"^\bullet ONOO"
    # assert formula_to_latex(".NO3/2-") == r"^\bullet NO_{3}^{2-}"
    assert formula_to_latex(".NO3-2") == r"^\bullet NO_{3}^{2-}"
    assert formula_to_latex("alpha-FeOOH(s)") == r"\alpha-FeOOH(s)"
    assert formula_to_latex("epsilon-Zn(OH)2(s)") == (r"\varepsilon-Zn(OH)_{2}(s)")
    assert formula_to_latex("Na2CO3..7H2O(s)") == r"Na_{2}CO_{3}\cdot 7H_{2}O(s)"
    assert formula_to_latex("Na2CO3..1H2O(s)") == r"Na_{2}CO_{3}\cdot H_{2}O(s)"
    assert formula_to_latex("K4[Fe(CN)6]") == r"K_{4}[Fe(CN)_{6}]"
    assert formula_to_latex("K4[Fe(CN)6](s)") == r"K_{4}[Fe(CN)_{6}](s)"
    assert formula_to_latex("K4[Fe(CN)6](aq)") == r"K_{4}[Fe(CN)_{6}](aq)"
    assert (
        formula_to_latex("[Fe(H2O)6][Fe(CN)6]..19H2O")
        == r"[Fe(H_{2}O)_{6}][Fe(CN)_{6}]\cdot 19H_{2}O"
    )
    assert (
        formula_to_latex("[Fe(H2O)6][Fe(CN)6]..19H2O(s)")
        == r"[Fe(H_{2}O)_{6}][Fe(CN)_{6}]\cdot 19H_{2}O(s)"
    )
    assert (
        formula_to_latex("[Fe(H2O)6][Fe(CN)6]..19H2O(aq)")
        == r"[Fe(H_{2}O)_{6}][Fe(CN)_{6}]\cdot 19H_{2}O(aq)"
    )
    assert formula_to_latex("[Fe(CN)6]-3") == r"[Fe(CN)_{6}]^{3-}"
    assert formula_to_latex("[Fe(CN)6]-3(aq)") == r"[Fe(CN)_{6}]^{3-}(aq)"
    assert (
        formula_to_latex("Ca2.832Fe0.6285Mg5.395(CO3)6")
        == r"Ca_{2.832}Fe_{0.6285}Mg_{5.395}(CO_{3})_{6}"
    )
    assert (
        formula_to_latex("Ca2.832Fe0.6285Mg5.395(CO3)6(s)")
        == r"Ca_{2.832}Fe_{0.6285}Mg_{5.395}(CO_{3})_{6}(s)"
    )
    assert (
        formula_to_latex("Ca2.832Fe0.6285Mg5.395(CO3)6..8H2O(s)")
        == r"Ca_{2.832}Fe_{0.6285}Mg_{5.395}(CO_{3})_{6}\cdot 8H_{2}O(s)"
    )


@requires(parsing_library)
def test_formula_to_latex_braces():
    # Parentheses.
    assert formula_to_latex("Fe(CN)6-3") == "Fe(CN)_{6}^{3-}"
    assert formula_to_latex("((H2O)2OH)12") == "((H_{2}O)_{2}OH)_{12}"

    # Square brackets.
    assert formula_to_latex("Fe[CN]6-3") == "Fe[CN]_{6}^{3-}"
    assert formula_to_latex("[(H2O)2OH]12") == "[(H_{2}O)_{2}OH]_{12}"

    # Curly braces.
    assert formula_to_latex("Fe{CN}6-3") == r"Fe\{CN\}_{6}^{3-}"
    assert formula_to_latex("{(H2O)2OH}12") == r"\{(H_{2}O)_{2}OH\}_{12}"


@requires(parsing_library)
def test_formula_to_latex_caged():
    """Should produce LaTeX for cage species."""
    assert formula_to_latex("Li@C60") == r"Li@C_{60}"
    assert formula_to_latex("(Li@C60)+") == r"(Li@C_{60})^{+}"
    assert formula_to_latex("Na@C60") == r"Na@C_{60}"
    assert formula_to_latex("(Na@C60)+") == r"(Na@C_{60})^{+}"


@requires(parsing_library)
def test_formula_to_unicode():
    assert formula_to_unicode("NH4+") == u"NH₄⁺"
    assert formula_to_unicode("H2O") == u"H₂O"
    # assert formula_to_unicode("C6H6/+") == u"C₆H₆⁺"
    assert formula_to_unicode("C6H6+") == u"C₆H₆⁺"
    # assert formula_to_unicode("Fe(CN)6/3-") == u"Fe(CN)₆³⁻"
    assert formula_to_unicode("Fe(CN)6-3") == u"Fe(CN)₆³⁻"
    # assert formula_to_unicode("C18H38/2+") == u"C₁₈H₃₈²⁺"
    # assert formula_to_unicode("C18H38/+2") == u"C₁₈H₃₈²⁺"
    assert formula_to_unicode("C18H38+2") == u"C₁₈H₃₈²⁺"
    assert formula_to_unicode("((H2O)2OH)12") == u"((H₂O)₂OH)₁₂"
    assert formula_to_unicode("[(H2O)2OH]12") == u"[(H₂O)₂OH]₁₂"
    assert formula_to_unicode("{(H2O)2OH}12") == u"{(H₂O)₂OH}₁₂"
    assert formula_to_unicode("NaCl") == u"NaCl"
    assert formula_to_unicode("NaCl(s)") == u"NaCl(s)"
    assert formula_to_unicode("e-(aq)") == u"e⁻(aq)"
    assert formula_to_unicode("Ca+2(aq)") == u"Ca²⁺(aq)"
    assert formula_to_unicode(".NO2(g)") == u"⋅NO₂(g)"
    assert formula_to_unicode(".NH2") == u"⋅NH₂"
    assert formula_to_unicode("ONOOH") == u"ONOOH"
    assert formula_to_unicode(".ONOO") == u"⋅ONOO"
    # assert formula_to_unicode(".NO3/2-") == u"⋅NO₃²⁻"
    assert formula_to_unicode(".NO3-2") == u"⋅NO₃²⁻"
    assert formula_to_unicode("alpha-FeOOH(s)") == u"α-FeOOH(s)"
    assert formula_to_unicode("epsilon-Zn(OH)2(s)") == u"ε-Zn(OH)₂(s)"
    assert formula_to_unicode("Na2CO3..7H2O(s)") == u"Na₂CO₃·7H₂O(s)"
    assert formula_to_unicode("Na2CO3..1H2O(s)") == u"Na₂CO₃·H₂O(s)"
    assert formula_to_unicode("K4[Fe(CN)6]") == r"K₄[Fe(CN)₆]"
    assert formula_to_unicode("K4[Fe(CN)6](s)") == r"K₄[Fe(CN)₆](s)"
    assert formula_to_unicode("K4[Fe(CN)6](aq)") == r"K₄[Fe(CN)₆](aq)"
    assert (
        formula_to_unicode("[Fe(H2O)6][Fe(CN)6]..19H2O") == r"[Fe(H₂O)₆][Fe(CN)₆]·19H₂O"
    )
    assert (
        formula_to_unicode("[Fe(H2O)6][Fe(CN)6]..19H2O(s)")
        == r"[Fe(H₂O)₆][Fe(CN)₆]·19H₂O(s)"
    )
    assert (
        formula_to_unicode("[Fe(H2O)6][Fe(CN)6]..19H2O(aq)")
        == r"[Fe(H₂O)₆][Fe(CN)₆]·19H₂O(aq)"
    )
    assert formula_to_unicode("[Fe(CN)6]-3") == r"[Fe(CN)₆]³⁻"
    assert formula_to_unicode("[Fe(CN)6]-3(aq)") == r"[Fe(CN)₆]³⁻(aq)"


@requires(parsing_library)
def test_formula_to_unicode_caged():
    """Should produce LaTeX for cage species."""
    assert formula_to_unicode("Li@C60") == r"Li@C₆₀"
    assert formula_to_unicode("(Li@C60)+") == r"(Li@C₆₀)⁺"
    assert formula_to_unicode("Na@C60") == r"Na@C₆₀"
    assert formula_to_unicode("(Na@C60)+") == r"(Na@C₆₀)⁺"


@requires(parsing_library)
def test_formula_to_html():
    assert formula_to_html("H2O") == "H<sub>2</sub>O"
    # assert formula_to_html("C6H6/+") == "C<sub>6</sub>H<sub>6</sub><sup>+</sup>"
    assert formula_to_html("C6H6+") == "C<sub>6</sub>H<sub>6</sub><sup>+</sup>"
    # assert formula_to_html("Fe(CN)6/3-") == "Fe(CN)<sub>6</sub><sup>3-</sup>"
    assert formula_to_html("Fe(CN)6-3") == "Fe(CN)<sub>6</sub><sup>3-</sup>"
    # assert formula_to_html("C18H38/2+") == "C<sub>18</sub>H<sub>38</sub><sup>2+</sup>"
    # assert formula_to_html("C18H38/+2") == "C<sub>18</sub>H<sub>38</sub><sup>2+</sup>"
    assert formula_to_html("C18H38+2") == "C<sub>18</sub>H<sub>38</sub><sup>2+</sup>"
    assert (
        formula_to_html("((H2O)2OH)12")
        == "((H<sub>2</sub>O)<sub>2</sub>OH)<sub>12</sub>"
    )
    assert (
        formula_to_html("[(H2O)2OH]12")
        == "[(H<sub>2</sub>O)<sub>2</sub>OH]<sub>12</sub>"
    )
    assert (
        formula_to_html("{(H2O)2OH}12")
        == "{(H<sub>2</sub>O)<sub>2</sub>OH}<sub>12</sub>"
    )
    assert formula_to_html("NaCl") == "NaCl"
    assert formula_to_html("NaCl(s)") == "NaCl(s)"
    assert formula_to_html("e-(aq)") == "e<sup>-</sup>(aq)"
    assert formula_to_html("Ca+2(aq)") == "Ca<sup>2+</sup>(aq)"
    assert formula_to_html(".NO2(g)") == r"&sdot;NO<sub>2</sub>(g)"
    assert formula_to_html(".NH2") == r"&sdot;NH<sub>2</sub>"
    assert formula_to_html("ONOOH") == "ONOOH"
    assert formula_to_html(".ONOO") == r"&sdot;ONOO"
    # assert formula_to_html(".NO3/2-") == r"&sdot;NO<sub>3</sub><sup>2-</sup>"
    assert formula_to_html(".NO3-2") == r"&sdot;NO<sub>3</sub><sup>2-</sup>"
    assert formula_to_html("alpha-FeOOH(s)") == r"&alpha;-FeOOH(s)"
    assert formula_to_html("epsilon-Zn(OH)2(s)") == (r"&epsilon;-Zn(OH)<sub>2</sub>(s)")
    assert (
        formula_to_html("Na2CO3..7H2O(s)")
        == "Na<sub>2</sub>CO<sub>3</sub>&sdot;7H<sub>2</sub>O(s)"
    )
    assert (
        formula_to_html("Na2CO3..1H2O(s)")
        == "Na<sub>2</sub>CO<sub>3</sub>&sdot;H<sub>2</sub>O(s)"
    )
    assert formula_to_html("K4[Fe(CN)6]") == r"K<sub>4</sub>[Fe(CN)<sub>6</sub>]"
    assert formula_to_html("K4[Fe(CN)6](s)") == r"K<sub>4</sub>[Fe(CN)<sub>6</sub>](s)"
    assert (
        formula_to_html("K4[Fe(CN)6](aq)") == r"K<sub>4</sub>[Fe(CN)<sub>6</sub>](aq)"
    )
    assert (
        formula_to_html("[Fe(H2O)6][Fe(CN)6]..19H2O")
        == r"[Fe(H<sub>2</sub>O)<sub>6</sub>][Fe(CN)<sub>6</sub>]&sdot;19H<sub>2</sub>O"
    )
    assert (
        formula_to_html("[Fe(H2O)6][Fe(CN)6]..19H2O(s)")
        == r"[Fe(H<sub>2</sub>O)<sub>6</sub>][Fe(CN)<sub>6</sub>]&sdot;19H<sub>2</sub>O(s)"
    )
    assert (
        formula_to_html("[Fe(H2O)6][Fe(CN)6]..19H2O(aq)")
        == r"[Fe(H<sub>2</sub>O)<sub>6</sub>][Fe(CN)<sub>6</sub>]&sdot;19H<sub>2</sub>O(aq)"
    )
    assert formula_to_html("[Fe(CN)6]-3") == r"[Fe(CN)<sub>6</sub>]<sup>3-</sup>"
    assert (
        formula_to_html("[Fe(CN)6]-3(aq)") == r"[Fe(CN)<sub>6</sub>]<sup>3-</sup>(aq)"
    )


@requires(parsing_library)
def test_formula_to_html_caged():
    """Should produce HTML for cage species."""
    assert formula_to_html("Li@C60") == r"Li@C<sub>60</sub>"
    assert formula_to_html("(Li@C60)+") == r"(Li@C<sub>60</sub>)<sup>+</sup>"
    assert formula_to_html("Na@C60") == r"Na@C<sub>60</sub>"
    assert formula_to_html("(Na@C60)+") == r"(Na@C<sub>60</sub>)<sup>+</sup>"
