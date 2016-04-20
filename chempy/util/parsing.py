# -*- coding: utf-8 -*-
""" Functions for chemical formulae and reactions """

from __future__ import (absolute_import, division, print_function)

from collections import defaultdict

import re
import warnings

from .pyutil import ChemPyDeprecationWarning, memoize

parsing_library = 'pyparsing'  # info used for selective testing.


@memoize
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
    _p = __import__(parsing_library)
    Forward, Group, OneOrMore = _p.Forward, _p.Group, _p.OneOrMore
    Optional, ParseResults, Regex = _p.Optional, _p.ParseResults, _p.Regex
    Suppress, Word, nums = _p.Suppress, _p.Word, _p.nums

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

lower_names = tuple(n.lower().lstrip('(').rstrip(')') for n in names)


def atomic_number(name):
    try:
        return symbols.index(name) + 1
    except ValueError:
        return lower_names.index(name.lower()) + 1

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
            yield float(mass[1:-1])
        elif '(' in mass:
            yield float(mass.split('(')[0])
        else:
            yield(float(mass))

relative_atomic_masses = tuple(_get_relative_atomic_masses())


def mass_from_composition(composition):
    """ Calculates molecular mass from atomic weights

    Parameters
    ----------
    composition: dict
        Dictionary mapping int (atomic number) to int (coefficient)

    Returns
    -------
    float
        molecular weight in atomic mass units


    Notes
    -----
    Atomic number 0 denotes charge or "net electron defficiency"

    Examples
    --------
    >>> '%.2f' % mass_from_composition({0: -1, 1: 1, 8: 1})
    '17.01'
    """
    mass = 0.0
    for k, v in composition.items():
        if k == 0:  # electron
            mass -= v*5.489e-4
        else:
            mass += v*relative_atomic_masses[k-1]
    return mass


def _get_charge(chgstr):

    if chgstr == '+':
        return 1
    elif chgstr == '-':
        return -1

    for token, anti, sign in zip('+-', '-+', (1, -1)):
        if token in chgstr:
            if anti in chgstr:
                raise ValueError("Invalid charge description (+ & - present)")
            before, after = chgstr.split(token)
            if len(before) > 0 and len(after) > 0:
                raise ValueError("Values both before and after charge token")
            if len(before) > 0:
                # will_be_missing_in='0.4.0'
                warnings.warn("'Fe/3+' deprecated, use e.g. 'Fe+3'",
                              ChemPyDeprecationWarning, stacklevel=3)
                return sign * int(1 if before == '' else before)
            if len(after) > 0:
                return sign * int(1 if after == '' else after)
    raise ValueError("Invalid charge description (+ or - missing)")


def _formula_to_parts(formula, prefixes, suffixes):
    # Drop prefixes and suffixes
    drop_pref, drop_suff = [], []
    for ign in prefixes:
        if formula.startswith(ign):
            drop_pref.append(ign)
            formula = formula[len(ign):]
    for ign in suffixes:
        if formula.endswith(ign):
            drop_suff.append(ign)
            formula = formula[:-len(ign)]

    # Extract charge
    if '/' in formula:
        # will_be_missing_in='0.4.0'
        warnings.warn("/ depr. (before 0.4.0): use 'Fe+3' over 'Fe/3+'",
                      ChemPyDeprecationWarning, stacklevel=3)
        parts = formula.split('/')

        if '+' in parts[0] or '-' in parts[0]:
            raise ValueError("Charge needs to be separated with a /")
        if parts[1] is not None:
            wo_pm = parts[1].replace('+', '').replace('-', '')
            if wo_pm != '' and not str.isdigit(wo_pm):
                raise ValueError("Non-digits in charge specifier")
        if len(parts) > 2:
            raise ValueError("At most one '/' allowed in formula")
    else:
        for token in '+-':
            if token in formula:
                if formula.count(token) > 1:
                    raise ValueError("Multiple tokens: %s" % token)
                parts = formula.split(token)
                parts[1] = token + parts[1]
                break
        else:
            parts = [formula, None]
    return parts + [tuple(drop_pref), tuple(drop_suff[::-1])]


def _parse_stoich(stoich):
    if stoich == 'e':  # special case, the electron is not an element
        return {}
    return {symbols.index(k)+1: n for k, n
            in _get_formula_parser().parseString(stoich)}

_greek_letters = (
    'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta',
    'iota', 'kappa', 'lambda', 'mu', 'nu', 'xi', 'omicron', 'pi', 'rho',
    'sigma', 'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega'
)
_greek_u = u'αβγδεζηθικλμνξοπρστυφχψω'

_latex_mapping = {k + '-': '\\' + k + '-' for k in _greek_letters}
_latex_mapping['epsilon-'] = '\\varepsilon-'
_latex_mapping['omicron-'] = 'o-'
_latex_mapping['.'] = '^\\bullet '

_unicode_mapping = {k + '-': v + '-' for k, v in zip(_greek_letters, _greek_u)}
_unicode_mapping['.'] = u'⋅'

_html_mapping = {k + '-': '&' + k + ';-' for k in _greek_letters}
_html_mapping['.'] = '&sdot;'


def formula_to_composition(formula, prefixes=None,
                           suffixes=('(s)', '(l)', '(g)', '(aq)')):
    """ Parse composition of formula representing a chemical formula

    Composition is represented as a dict mapping int -> int (atomic
    number -> multiplicity). "Atomic number" 0 represents net charge.

    Parameters
    ----------
    formula: str
        Chemical formula, e.g. 'H2O', 'Fe+3', 'Cl-'
    prefixes: iterable strings
        Prefixes to ignore, e.g. ('.', 'alpha-')
    suffixes: tuple of strings
        Suffixes to ignore, e.g. ('(g)', '(s)')

    Examples
    --------
    >>> formula_to_composition('NH4+') == {0: 1, 1: 4, 7: 1}
    True
    >>> formula_to_composition('.NHO-(aq)') == {0: -1, 1: 1, 7: 1, 8: 1}
    True

    """
    if prefixes is None:
        prefixes = _latex_mapping.keys()
    parts = _formula_to_parts(formula, prefixes, suffixes)
    comp = _parse_stoich(parts[0])
    if parts[1] is not None:
        comp[0] = _get_charge(parts[1])
    return comp


def _subs(string, patterns):
    for patt, repl in patterns.items():
        string = string.replace(patt, repl)
    return string


def _parse_multiplicity(strings, substance_keys=None):
    """
    Examples
    --------
    >>> _parse_multiplicity(['2 H2O2', 'O2']) == {'H2O2': 2, 'O2': 1}
    True
    >>> _parse_multiplicity(['2 * H2O2', 'O2']) == {'H2O2': 2, 'O2': 1}
    True

    """
    result = {}
    for items in [re.split(' \* | ', s) for s in strings]:
        if len(items) == 1:
            result[items[0]] = 1
        elif len(items) == 2:
            result[items[1]] = int(items[0])
        else:
            raise ValueError("To many parts in substring")
    if substance_keys is not None:
        for k in result:
            if k not in substance_keys:
                raise ValueError("Unkown substance_key: %s" % k)
    return result


def to_reaction(line, substance_keys, token, Cls, globals_=None):
    """ Parses a string into a Reaction object and substances

    Reac1 + 2 Reac2 + (2 Reac1) -> Prod1 + Prod2; 10**3.7; ref='doi:12/ab'
    Reac1 = Prod1; 2.1;

    Parameters
    ----------
    line: str
        string representation to be parsed
    substance_keys: iterable of strings
        Allowed names, e.g. ('H2O', 'H+', 'OH-')
    token : str
        delimiter token between reactant and product side
    Cls : class
        e.g. subclass of Reaction
    globals_: dict (optional)
        Globals passed on to :func:`eval`, when ``None``:
        `chempy.units.default_units` is used with 'chempy'
        and 'default_units' extra entries.

    Notes
    -----
    This function calls :func:`eval`, hence there are severe security concerns
    with running this on untrusted data.

    """
    # TODO: add handling of units.
    if globals_ is None:
        import chempy
        from chempy.units import default_units
        globals_ = {'chempy': chempy, 'default_units': default_units}
        globals_.update(default_units.as_dict())

    try:
        stoich, param, kwargs = map(str.strip, line.rstrip('\n').split(';'))
    except ValueError:
        if ';' in line:
            stoich, param = map(str.strip, line.rstrip('\n').split(';'))
        else:
            stoich, param = line.strip(), 'None'
        kwargs = {}
    else:
        kwargs = eval('dict('+kwargs+')', globals_)

    param = eval(param, globals_)

    if token not in stoich:
        raise ValueError("Missing token: %s" % token)

    reac_prod = [[y.strip() for y in x.split(' + ')] for
                 x in stoich.split(token)]

    act, inact = [], []
    for side in reac_prod:
        if side[-1].startswith('('):
            if not side[-1].endswith(')'):
                raise ValueError("Bad format (missing closing paren)")
            inact.append(_parse_multiplicity(side[-1][1:-1].split(' + '),
                                             substance_keys))
            act.append(_parse_multiplicity(side[:-1], substance_keys))
        else:
            inact.append({})
            act.append(_parse_multiplicity(side, substance_keys))

    # stoich coeff -> dict
    return Cls(act[0], act[1], param, inact_reac=inact[0],
               inact_prod=inact[1], **kwargs)


def _formula_to_format(sub, sup, formula, prefixes=None,
                       suffixes=('(s)', '(l)', '(g)', '(aq)'),
                       intercept=None):
    parts = _formula_to_parts(formula, prefixes.keys(), suffixes)
    string = re.sub(r'([0-9]+)', lambda m: sub(m.group(1)), parts[0])
    if parts[1] is not None:
        chg = _get_charge(parts[1])
        if chg < 0:
            token = '-' if chg == -1 else '%d-' % -chg
        if chg > 0:
            token = '+' if chg == 1 else '%d+' % chg
        string += sup(token)
    if len(parts) > 4:
        raise ValueError("Incorrect formula")
    pre_str = ''.join(map(lambda x: _subs(x, prefixes), parts[2]))
    return pre_str + string + ''.join(parts[3])


def formula_to_latex(formula, prefixes=None, **kwargs):
    r""" Convert formula string to latex representation

    Parameters
    ----------
    formula: str
        Chemical formula, e.g. 'H2O', 'Fe+3', 'Cl-'
    prefixes: dict
        Prefix transofmrations, default: greek letters and .
    suffixes: tuple of strings
        Suffixes to keep, e.g. ('(g)', '(s)')

    Examples
    --------
    >>> formula_to_latex('NH4+')
    'NH_{4}^{+}'
    >>> formula_to_latex('Fe(CN)6+2')
    'Fe(CN)_{6}^{2+}'
    >>> formula_to_latex('Fe(CN)6+2(aq)')
    'Fe(CN)_{6}^{2+}(aq)'
    >>> formula_to_latex('.NHO-(aq)')
    '^\\bullet NHO^{-}(aq)'
    >>> formula_to_latex('alpha-FeOOH(s)')
    '\\alpha-FeOOH(s)'

    """
    if prefixes is None:
        prefixes = _latex_mapping
    return _formula_to_format(lambda x: '_{%s}' % x, lambda x: '^{%s}' % x,
                              formula, prefixes)


def _number_to_scientific_latex(number, fmt='%.3g'):
    r"""
    Examples
    --------
    >>> _number_to_scientific_latex(3.14) == '3.14'
    True
    >>> _number_to_scientific_latex(3.14159265e-7)
    '3.14\\cdot 10^{-7}'

    """
    s = fmt % number
    if 'e' in s:
        prefix, suffix = s.split('e')
        return prefix + r'\cdot 10^{%s}' % str(int(suffix))
    else:
        return s

_unicode_sub = {}

for k, v in enumerate(u"₀₁₂₃₄₅₆₇₈₉"):
    _unicode_sub[str(k)] = v

_unicode_sup = {
    '+': u'⁺',
    '-': u'⁻',
}

for k, v in enumerate(u"⁰¹²³⁴⁵⁶⁷⁸⁹"):
    _unicode_sup[str(k)] = v


def formula_to_unicode(formula, prefixes=None, **kwargs):
    u""" Convert formula string to unicode string representation

    Parameters
    ----------
    formula: str
        Chemical formula, e.g. 'H2O', 'Fe+3', 'Cl-'
    prefixes: dict
        Prefix transofmrations, default: greek letters and .
    suffixes: tuple of strings
        Suffixes to keep, e.g. ('(g)', '(s)')

    Examples
    --------
    >>> formula_to_unicode('NH4+') == u'NH₄⁺'
    True
    >>> formula_to_unicode('Fe(CN)6+2') == u'Fe(CN)₆²⁺'
    True
    >>> formula_to_unicode('Fe(CN)6+2(aq)') == u'Fe(CN)₆²⁺(aq)'
    True
    >>> formula_to_unicode('.NHO-(aq)') == u'⋅NHO⁻(aq)'
    True
    >>> formula_to_unicode('alpha-FeOOH(s)') == u'α-FeOOH(s)'
    True

    """
    if prefixes is None:
        prefixes = _unicode_mapping

    return _formula_to_format(
        lambda x: ''.join(_unicode_sub[str(_)] for _ in x),
        lambda x: ''.join(_unicode_sup[str(_)] for _ in x),
        formula, prefixes)


def _number_to_scientific_unicode(number, fmt='%.3g'):
    u"""
    Examples
    --------
    >>> _number_to_scientific_unicode(3.14) == u'3.14'
    True
    >>> _number_to_scientific_unicode(3.14159265e-7) == u'3.14·10⁻⁷'
    True

    """
    s = fmt % number
    if 'e' in s:
        prefix, suffix = s.split('e')
        return prefix + u'·10' + u''.join(map(_unicode_sup.get, str(int(suffix))))
    else:
        return s


def formula_to_html(formula, prefixes=None, **kwargs):
    u""" Convert formula string to html string representation

    Parameters
    ----------
    formula: str
        Chemical formula, e.g. 'H2O', 'Fe+3', 'Cl-'
    prefixes: dict
        Prefix transofmrations, default: greek letters and .
    suffixes: tuple of strings
        Suffixes to keep, e.g. ('(g)', '(s)')

    Examples
    --------
    >>> formula_to_html('NH4+')
    'NH<sub>4</sub><sup>+</sup>'
    >>> formula_to_html('Fe(CN)6+2')
    'Fe(CN)<sub>6</sub><sup>2+</sup>'
    >>> formula_to_html('Fe(CN)6+2(aq)')
    'Fe(CN)<sub>6</sub><sup>2+</sup>(aq)'
    >>> formula_to_html('.NHO-(aq)')
    '&sdot;NHO<sup>-</sup>(aq)'
    >>> formula_to_html('alpha-FeOOH(s)')
    '&alpha;-FeOOH(s)'

    """
    if prefixes is None:
        prefixes = _html_mapping
    return _formula_to_format(lambda x: '<sub>%s</sub>' % x,
                              lambda x: '<sup>%s</sup>' % x,
                              formula, prefixes)


def _number_to_scientific_html(number, fmt='%.3g'):
    """
    Examples
    --------
    >>> _number_to_scientific_html(3.14) == '3.14'
    True
    >>> _number_to_scientific_html(3.14159265e-7)
    '3.14&sdot;10<sup>-7</sup>'

    """
    s = fmt % number
    if 'e' in s:
        prefix, suffix = s.split('e')
        return prefix + '&sdot;10<sup>' + str(int(suffix)) + '</sup>'
    else:
        return s
