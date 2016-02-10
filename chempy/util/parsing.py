# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function)

from collections import defaultdict
import re
from pyparsing import (
    Forward, Group, OneOrMore, Optional, ParseResults,
    Regex, Suppress, Word, nums
)


def _get_formula_parser():
    """ Create a forward pyparsing parser for chemical formulae

    BNF for simple chemical formula (no nesting)

        integer :: '0'..'9'+
        element :: 'A'..'Z' 'a'..'z'*
        term :: element [integer]
        formula :: term+


    BNF for nested chemical formula

        integer :: '0'..'9'+
        element :: 'A'..'Z' 'a'..'z'*
        term :: (element | '(' formula ')') [integer]
        formula :: term+

    Notes
    -----
    The code in this function is from an answer on StackOverflow:
        http://stackoverflow.com/a/18555142/790973
        written by:
            Paul McGuire, http://stackoverflow.com/users/165216/paul-mcguire
        in answer to the question formulated by:
            Thales MG, http://stackoverflow.com/users/2708711/thales-mg
        the code is licensed under 'CC-WIKI'.
        (see: http://blog.stackoverflow.com/2009/06/attribution-required/)

    """

    LPAR, RPAR = map(Suppress, "()")
    integer = Word(nums)

    # add parse action to convert integers to ints, to support doing addition
    # and multiplication at parse time
    integer.setParseAction(lambda t: int(t[0]))

    # element = Word(alphas.upper(), alphas.lower())
    # or if you want to be more specific, use this Regex
    element = Regex(
        r"A[cglmrstu]|B[aehikr]?|C[adeflmorsu]?|D[bsy]|E[rsu]|F[emr]?|"
        "G[ade]|H[efgos]?|I[nr]?|Kr?|L[airu]|M[dgnot]|N[abdeiop]?|"
        "Os?|P[abdmortu]?|R[abefghnu]|S[bcegimnr]?|T[abcehilm]|"
        "Uu[bhopqst]|U|V|W|Xe|Yb?|Z[nr]")

    # forward declare 'formula' so it can be used in definition of 'term'
    formula = Forward()

    term = Group((element | Group(LPAR + formula + RPAR)("subgroup")) +
                 Optional(integer, default=1)("mult"))

    # define contents of a formula as one or more terms
    formula << OneOrMore(term)

    # add parse actions for parse-time processing

    # parse action to multiply out subgroups
    def multiplyContents(tokens):
        t = tokens[0]
        # if these tokens contain a subgroup, then use multiplier to
        # extend counts of all elements in the subgroup
        if t.subgroup:
            mult = t.mult
            for term in t.subgroup:
                term[1] *= mult
            return t.subgroup
    term.setParseAction(multiplyContents)

    # add parse action to sum up multiple references to the same element
    def sumByElement(tokens):
        elementsList = [t[0] for t in tokens]

        # construct set to see if there are duplicates
        duplicates = len(elementsList) > len(set(elementsList))

        # if there are duplicate element names, sum up by element and
        # return a new nested ParseResults
        if duplicates:
            ctr = defaultdict(int)
            for t in tokens:
                ctr[t[0]] += t[1]
            return ParseResults([ParseResults([k, v]) for k, v in ctr.items()])
    formula.setParseAction(sumByElement)

    return formula

_parser = _get_formula_parser()

symbols = (
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
    'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
    'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',
    'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W',
    'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
    'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
    'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo'
)

names = (
    'Hydrogen', 'Helium', 'Lithium', 'Beryllium',
    'Boron', 'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium',
    'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur',
    'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium', 'Titanium',
    'Vanadium', 'Chromium', 'Manganese', 'Iron', 'Cobalt', 'Nickel',
    'Copper', 'Zinc', 'Gallium', 'Germanium', 'Arsenic', 'Selenium',
    'Bromine', 'Krypton', 'Rubidium', 'Strontium', 'Yttrium', 'Zirconium',
    'Niobium', 'Molybdenum', 'Technetium', 'Ruthenium', 'Rhodium',
    'Palladium', 'Silver', 'Cadmium', 'Indium', 'Tin', 'Antimony',
    'Tellurium', 'Iodine', 'Xenon', 'Caesium', 'Barium', 'Lanthanum',
    'Cerium', 'Praseodymium', 'Neodymium', 'Promethium', 'Samarium',
    'Europium', 'Gadolinium', 'Terbium', 'Dysprosium', 'Holmium',
    'Erbium', 'Thulium', 'Ytterbium', 'Lutetium', 'Hafnium', 'Tantalum',
    'Tungsten', 'Rhenium', 'Osmium', 'Iridium', 'Platinum', 'Gold',
    'Mercury', 'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine',
    'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium', 'Protactinium',
    'Uranium', 'Neptunium', 'Plutonium', 'Americium', 'Curium',
    'Berkelium', 'Californium', 'Einsteinium', 'Fermium', 'Mendelevium',
    'Nobelium', 'Lawrencium', 'Rutherfordium', 'Dubnium', 'Seaborgium',
    'Bohrium', 'Hassium', 'Meitnerium', 'Darmstadtium', 'Roentgenium',
    'Copernicium', '(Ununtrium)', 'Flerovium', '(Ununpentium)',
    'Livermorium', '(Ununseptium)', '(Ununoctium)'
)

# The data in '_relative_atomic_masses' is licensed under the CC-SA license
# https://en.wikipedia.org/w/index.php?title=List_of_elements&oldid=700476748
_relative_atomic_masses = (
    "1.008 4.002602(2) 6.94 9.0121831(5) 10.81 12.011 14.007 15.999"
    " 18.998403163(6) 20.1797(6) 22.98976928(2) 24.305 26.9815385(7) 28.085"
    " 30.973761998(5) 32.06 35.45 39.948(1) 39.0983(1) 40.078(4)"
    " 44.955908(5) 47.867(1) 50.9415(1) 51.9961(6) 54.938044(3) 55.845(2)"
    " 58.933194(4) 58.6934(4) 63.546(3) 65.38(2) 69.723(1) 72.630(8)"
    " 74.921595(6) 78.971(8) 79.904 83.798(2) 85.4678(3) 87.62(1)"
    " 88.90584(2) 91.224(2) 92.90637(2) 95.95(1) [98] 101.07(2) 102.90550(2)"
    " 106.42(1) 107.8682(2) 112.414(4) 114.818(1) 118.710(7) 121.760(1)"
    " 127.60(3) 126.90447(3) 131.293(6) 132.90545196(6) 137.327(7)"
    " 138.90547(7) 140.116(1) 140.90766(2) 144.242(3) [145] 150.36(2)"
    " 151.964(1) 157.25(3) 158.92535(2) 162.500(1) 164.93033(2) 167.259(3)"
    " 168.93422(2) 173.045(10) 174.9668(1) 178.49(2) 180.94788(2) 183.84(1)"
    " 186.207(1) 190.23(3) 192.217(3) 195.084(9) 196.966569(5) 200.592(3)"
    " 204.38 207.2(1) 208.98040(1) [209] [210] [222] [223] [226] [227]"
    " 232.0377(4) 231.03588(2) 238.02891(3) [237] [244] [243] [247] [247]"
    " [251] [252] [257] [258] [259] [266] [267] [268] [269] [270] [269]"
    " [278] [281] [282] [285] [286] [289] [289] [293] [294] [294]"
)


def _get_relative_atomic_masses():
    for mass in _relative_atomic_masses.split():
        if mass.startswith('[') and mass.endswith(']'):
            yield float(mass[1:-2])
        elif '(' in mass:
            yield float(mass.split('(')[0])
        else:
            yield(float(mass))

relative_atomic_masses = tuple(_get_relative_atomic_masses())


def mass_from_composition(composition):
    mass = 0.0
    for k, v in composition.items():
        if k == 0:  # electron
            mass -= v*5.489e-4
        else:
            mass += v*relative_atomic_masses[k-1]
    return mass


def _get_charge(chgstr):
    if '-' in chgstr:
        if '+' in chgstr:
            raise ValueError("Invalid charge description (+ & - present)")
        if chgstr.count('-') != 1:
            raise ValueError("Invalid charge description (multiple -)")
        s = chgstr.replace('-', '')
        return -int(1 if s == '' else s)
    elif '+' in chgstr:
        if '-' in chgstr:
            raise ValueError("Invalid charge description (+ & - present)")
        if chgstr.count('+') != 1:
            raise ValueError("Invalid charge description (multiple +)")
        s = chgstr.replace('+', '')
        return int(1 if s == '' else s)
    else:
        raise ValueError("Invalid charge description (+ or - missing)")


def _formula_to_parts(formula):
    parts = formula.split('/')

    if '+' in parts[0] or '-' in parts[0]:
        raise ValueError("Charge needs to be separated with a /")
    if len(parts) > 2:
        raise ValueError("At most one '/' allowed in formula")
    return parts


def to_composition(formula):
    """ Parse composition of formula representing a chemical formula

    Composition is represented as a dict mapping int -> int (atomic
    number -> multiplicity). "Atomic number" 0 represents net charge.

    Parameters
    ----------
    formula: str
        Chemical formula, e.g. 'H2O', 'Fe/3+', 'Cl/-'

    Examples
    --------
    >>> parse_formula('NH4/+'):
    {0: 1, 1: 4, 7: 1}
    """
    parts = _formula_to_parts(formula)
    comp = {symbols.index(k)+1: n for k, n in _parser.parseString(parts[0])}
    if len(parts) == 2:
        comp[0] = _get_charge(parts[1])
    return comp


def to_latex(formula):
    """ Convert formula string to latex representation

    Parameters
    ----------
    formula: str
        Chemical formula, e.g. 'H2O', 'Fe/3+', 'Cl/-'

    Examples
    --------
    >>> to_latex('NH4/+')
    NH_{4}^{+}
    >>> to_latex('Fe(CN)6/+2')
    Fe(CN)_{6}^{2+}

    """
    parts = _formula_to_parts(formula)

    string = re.sub(r'([0-9]+)', r'_{\1}', parts[0])
    if len(parts) == 2:
        chg = _get_charge(parts[1])
        if chg < 0:
            token = '-' if chg == -1 else '%d-' % -chg
        if chg > 0:
            token = '+' if chg == 1 else '%d+' % chg
        string += '^{%s}' % token
    return string
