# -*- coding: utf-8 -*-
""" Functions for chemical formulae and reactions """


from collections import defaultdict

import re

from .pyutil import memoize
from .periodic import symbols

parsing_library = "pyparsing"  # info used for selective testing.


def get_parsing_context():
    """returns the default dictionary for parsing strings in chempy"""
    import chempy
    from chempy.kinetics import rates
    from chempy.units import default_units, default_constants, to_unitless

    globals_ = dict(to_unitless=to_unitless, chempy=chempy)

    def _update(mod, keys=None):
        if keys is None:
            keys = dir(mod)
        globals_.update({k: getattr(mod, k) for k in keys if not k.startswith("_")})

    try:
        import numpy
    except ImportError:

        def _numpy_not_installed_raise(*args, **kwargs):
            raise ImportError("numpy not installed, no such method")

        class numpy:
            array = staticmethod(_numpy_not_installed_raise)
            log = staticmethod(_numpy_not_installed_raise)
            exp = staticmethod(_numpy_not_installed_raise)

    _update(numpy, keys="array log exp".split())  # could of course add more
    _update(rates)
    _update(chempy)
    for df in [default_units, default_constants]:
        if df is not None:
            globals_.update(df.as_dict())
    return globals_


@memoize()
def _get_formula_parser():
    """Create a forward pyparsing parser for chemical formulae

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

    Documentation for the desired product.  Original documentation
    above.

    Create a chemical formula parser.

    Parse a chemical formula, including elements, nested ions,
    complexes, charges (ions), hydrates, and state symbols.

    BNF for nested chemical formula with complexes

        count :: ( '1'..'9'? | '1'..'9'' '0'..'9'+ )
        element :: 'A'..'Z' 'a'..'z'*
        charge :: ( '-' | '+' ) ( '1'..'9'? | '1'..'9'' '0'..'9'+ )
        prime :: ( "*" | "'" )*
        term :: (element
                 | '(' formula ')'
                 | '{' formula '}'
                 | '[' formula ']' ) count prime charge?
        formula :: term+
        hydrate :: ( '..' | '\u00B7' ) count? formula
        state :: '(' ( 's' | 'l' | 'g' | 'aq' | 'cr' ) ')'
        compound :: count formula hydrate? state?

    Parse a chemical formula, including elements, non-integer
    subscripts, nested ions, complexes, charges (ions), hydrates, and
    state symbols.

    BNF for nested chemical formula with complexes

        count :: ( '1'..'9'? | '1'..'9'' '0'..'9'+ | '0'..'9'+ '.' '0'..'9'+ )
        element :: 'A'..'Z' 'a'..'z'*
        charge :: ( '-' | '+' ) ( '1'..'9'? | '1'..'9'' '0'..'9'+ )
        prime :: ( "*" | "'" )*
        term :: (element
                 | '(' formula ')'
                 | '{' formula '}'
                 | '[' formula ']' ) count prime charge?
        formula :: term+
        hydrate :: ( '..' | '\u00B7' ) count? formula
        state :: '(' ( 's' | 'l' | 'g' | 'aq' | 'cr' ) ')'
        compound :: count formula hydrate? state?
    """
    _p = __import__(parsing_library)
    Forward, Group, OneOrMore = _p.Forward, _p.Group, _p.OneOrMore
    Optional, ParseResults, Regex = _p.Optional, _p.ParseResults, _p.Regex
    Suppress = _p.Suppress

    # Define and suppress the grouping symbols.
    LCB = Suppress(Regex(r"\{"))
    RCB = Suppress(Regex(r"\}"))
    LSB = Suppress(Regex(r"\["))
    RSB = Suppress(Regex(r"\]"))
    LP = Suppress(Regex(r"\("))
    RP = Suppress(Regex(r"\)"))

    # Define and suppress the caged symbol.
    caged = Suppress(Regex(r"\@"))

    # Primes/stars for marking special species in reactions.
    primes = Suppress(Regex(r"[*']+"))

    # Parse counts (subscripts and coefficients).
    count = Regex(r"(\d+\.\d+|\d*)")
    count.setParseAction(lambda t: 1 if t[0] == "" else float(t[0]))

    # Parse states.
    state = Suppress(Regex(r"\((s|l|g|aq|cr)\)"))

    # Elements, 1-118, official symbols.
    element = Regex(
        r"A[cglmrstu]"
        "|B[aehikr]?"
        "|C[adeflmnorsu]?"
        "|D[bsy]"
        "|E[rsu]"
        "|F[elmr]?"
        "|G[ade]"
        "|H[efgos]?"
        "|I[nr]?"
        "|Kr?"
        "|L[airuv]"
        "|M[cdgnot]"
        "|N[abdehiop]?"
        "|O[gs]?"
        "|P[abdmortu]?"
        "|R[abefghnu]"
        "|S[bcegimnr]?"
        "|T[abcehilms]"
        "|U"
        "|V"
        "|W"
        "|Xe"
        "|Yb?"
        "|Z[nr]"
    ).setResultsName("element", listAllMatches=True)

    # forward declare 'formula' so it can be used in definition of 'term'
    formula = Forward()

    term = Group(
        (
            element
            | Group(LP + formula + RP)("subgroup")
            | Group(LSB + formula + RSB)("subgroup")
            | Group(LCB + formula + RCB)("subgroup")
            | Group(caged + formula)("subgroup")
        )
        + Optional(count, default=1)("mult")
        + Optional(state)("state")
        + Optional(primes)("primes")
    )

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

    # define contents of a formula as one or more terms
    formula << OneOrMore(term)
    formula.setParseAction(sumByElement)

    return formula


def _get_charge(chgstr):

    if chgstr == "+":
        return 1
    elif chgstr == "-":
        return -1

    for token, anti, sign in zip("+-", "-+", (1, -1)):
        if token in chgstr:
            if anti in chgstr:
                raise ValueError("Invalid charge description (+ & - present)")

            before, after = chgstr.split(token)

            if len(before) > 0 and len(after) > 0:
                raise ValueError("Values both before and after charge token")

            if len(after) > 0:
                return sign * int(1 if after == "" else after)

    raise ValueError("Invalid charge description (+ or - missing)")


def _formula_to_parts(formula, prefixes, suffixes):
    # Drop prefixes and suffixes.
    drop_pref, drop_suff = [], []
    for ign in prefixes:
        if formula.startswith(ign):
            drop_pref.append(ign)
            formula = formula[len(ign) :]
    for ign in suffixes:
        if formula.endswith(ign):
            drop_suff.append(ign)
            formula = formula[: -len(ign)]

    # Extract charge.
    if "/" in formula:
        raise ValueError(
            "Slashes ('/') in charge strings are deprecated."
            "  Use `Fe+3` instead of `Fe/3+`."
        )
    else:
        for token in "+-":
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
    # Special case:  the electron is not an element.
    if stoich == "e":
        return {}

    comp = {}
    for k, n in _get_formula_parser().parseString(stoich, parseAll=True):
        # Only use rational subscripts if necessary as
        # ``sympy.linsolve()`` does not like non-integers when
        # balancing reactions.
        if n == int(n):
            comp[symbols.index(k) + 1] = int(n)
        else:
            comp[symbols.index(k) + 1] = n

    return comp


_greek_letters = (
    "alpha",
    "beta",
    "gamma",
    "delta",
    "epsilon",
    "zeta",
    "eta",
    "theta",
    "iota",
    "kappa",
    "lambda",
    "mu",
    "nu",
    "xi",
    "omicron",
    "pi",
    "rho",
    "sigma",
    "tau",
    "upsilon",
    "phi",
    "chi",
    "psi",
    "omega",
)
_greek_u = "αβγδεζηθικλμνξοπρστυφχψω"

_latex_mapping = {k + "-": "\\" + k + "-" for k in _greek_letters}
_latex_mapping["epsilon-"] = "\\varepsilon-"
_latex_mapping["omicron-"] = "o-"
_latex_mapping["."] = "^\\bullet "
_latex_infix_mapping = {"..": "\\cdot "}

_unicode_mapping = {k + "-": v + "-" for k, v in zip(_greek_letters, _greek_u)}
_unicode_mapping["."] = "⋅"
_unicode_infix_mapping = {"..": "\u00b7"}  # 0x00b7: '·'

_html_mapping = {k + "-": "&" + k + ";-" for k in _greek_letters}
_html_mapping["."] = "&sdot;"
# _html_infix_mapping = _html_mapping
_html_infix_mapping = {"..": "&sdot;"}


def _get_leading_integer(s):
    m = re.findall(r"^\d+", s)
    if len(m) == 0:
        m = 1
    elif len(m) == 1:
        s = s[len(m[0]) :]
        m = int(m[0])
    else:
        raise ValueError("Failed to parse: %s" % s)
    return m, s


def formula_to_composition(
    formula, prefixes=None, suffixes=("(s)", "(l)", "(g)", "(aq)")
):
    """Parse composition of formula representing a chemical formula

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
    >>> formula_to_composition('Na2CO3..7H2O') == {11: 2, 6: 1, 8: 10, 1: 14}
    True
    >>> formula_to_composition('UO2.3') == {92: 1, 8: 2.3}
    True

    """
    if prefixes is None:
        prefixes = _latex_mapping.keys()

    stoich_tok, chg_tok = _formula_to_parts(formula, prefixes, suffixes)[:2]
    tot_comp = {}
    if '\u00b7' in stoich_tok:
        parts = stoich_tok.split('\u00b7')
    else:
        parts = stoich_tok.split("..")

    for idx, stoich in enumerate(parts):
        if idx == 0:
            m = 1
        else:
            m, stoich = _get_leading_integer(stoich)
        comp = _parse_stoich(stoich)
        for k, v in comp.items():
            if k not in tot_comp:
                tot_comp[k] = m * v
            else:
                tot_comp[k] += m * v

    if chg_tok is not None:
        tot_comp[0] = _get_charge(chg_tok)

    return tot_comp


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
    >>> _parse_multiplicity(['']) == {}
    True
    >>> _parse_multiplicity(['H2O', 'H2O']) == {'H2O': 2}
    True

    """
    result = {}
    for items in [re.split(" \\* | ", s) for s in strings]:
        items = [x for x in items if x != ""]
        if len(items) == 0:
            continue
        elif len(items) == 1:
            if items[0] not in result:
                result[items[0]] = 0
            result[items[0]] += 1
        elif len(items) == 2:
            if items[1] not in result:
                result[items[1]] = 0
            result[items[1]] += (
                float(items[0]) if "." in items[0] or "e" in items[0] else int(items[0])
            )
        else:
            raise ValueError("To many parts in substring")
    if substance_keys is not None:
        for k in result:
            if k not in substance_keys:
                raise ValueError("Unknown substance_key: %s" % k)
    return result


def to_reaction(line, substance_keys, token, Cls, globals_=None, **kwargs):
    """Parses a string into a Reaction object and substances

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
    if globals_ is None:
        globals_ = get_parsing_context()
    parts = line.rstrip("\n").split(";")
    stoich = parts[0].strip()
    if len(parts) > 2:
        kwargs.update(eval("dict(" + ";".join(parts[2:]) + "\n)", globals_ or {}))
    if len(parts) > 1:
        param = parts[1].strip()
    else:
        param = kwargs.pop("param", "None")

    if isinstance(param, str):
        if param.startswith("'") and param.endswith("'") and "'" not in param[1:-1]:
            from ..kinetics.rates import MassAction
            from ._expr import Symbol

            param = MassAction(Symbol(unique_keys=(param[1:-1],)))
        else:
            param = None if globals_ is False else eval(param, globals_)

    if token not in stoich:
        raise ValueError("Missing token: %s" % token)

    reac_prod = [[y.strip() for y in x.split(" + ")] for x in stoich.split(token)]

    act, inact = [], []
    for elements in reac_prod:
        act.append(
            _parse_multiplicity(
                [x for x in elements if not x.startswith("(")], substance_keys
            )
        )
        inact.append(
            _parse_multiplicity(
                [x[1:-1] for x in elements if x.startswith("(") and x.endswith(")")],
                substance_keys,
            )
        )

    # stoich coeff -> dict
    return Cls(
        act[0], act[1], param, inact_reac=inact[0], inact_prod=inact[1], **kwargs
    )


def _formula_to_format(
    sub,
    sup,
    formula,
    prefixes=None,
    infixes=None,
    suffixes=("(s)", "(l)", "(g)", "(aq)"),
):
    parts = _formula_to_parts(formula, prefixes.keys(), suffixes)
    if '\u00b7' in parts[0]:
        stoichs = parts[0].split('\u00b7')
    else:
        stoichs = parts[0].split("..")
    string = ""
    for idx, stoich in enumerate(stoichs):
        if idx == 0:
            m = 1
        else:
            m, stoich = _get_leading_integer(stoich)
            string += _subs("..", infixes)
        if m != 1:
            string += str(m)
        string += re.sub(r"([0-9]+\.[0-9]+|[0-9]+)", lambda m: sub(m.group(1)), stoich)

    if parts[1] is not None:
        chg = _get_charge(parts[1])
        if chg < 0:
            token = "-" if chg == -1 else "%d-" % -chg
        if chg > 0:
            token = "+" if chg == 1 else "%d+" % chg
        string += sup(token)
    if len(parts) > 4:
        raise ValueError("Incorrect formula")
    pre_str = "".join(map(lambda x: _subs(x, prefixes), parts[2]))
    return pre_str + string + "".join(parts[3])


def formula_to_latex(formula, prefixes=None, infixes=None, **kwargs):
    r"""Convert formula string to latex representation

    Parameters
    ----------
    formula: str
        Chemical formula, e.g. 'H2O', 'Fe+3', 'Cl-'
    prefixes: dict
        Prefix transformations, default: greek letters and .
    infixes: dict
        Infix transformations, default: .
    suffixes: iterable of str
        What suffixes not to interpret, default: (s), (l), (g), (aq)

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
    if infixes is None:
        infixes = _latex_infix_mapping
    return _formula_to_format(
        lambda x: "_{%s}" % x,
        lambda x: "^{%s}" % x,
        # formula,
        re.sub(r"([{}])", r"\\\1", formula) if re.search(r"[{}]", formula) else formula,
        prefixes,
        infixes,
        **kwargs
    )


_unicode_sub = {
    ".": ".",
}

for k, v in enumerate("₀₁₂₃₄₅₆₇₈₉."):
    _unicode_sub[str(k)] = v

_unicode_sup = {
    "+": "⁺",
    "-": "⁻",
}

for k, v in enumerate("⁰¹²³⁴⁵⁶⁷⁸⁹"):
    _unicode_sup[str(k)] = v


def formula_to_unicode(formula, prefixes=None, infixes=None, **kwargs):
    """Convert formula string to unicode string representation

    Parameters
    ----------
    formula : str
        Chemical formula, e.g. 'H2O', 'Fe+3', 'Cl-'
    prefixes : dict
        Prefix transofmrations, default: greek letters and .
    infixes : dict
        Infix transofmrations, default: .
    suffixes : tuple of strings
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
    if infixes is None:
        infixes = _unicode_infix_mapping
    return _formula_to_format(
        lambda x: "".join(_unicode_sub[str(_)] for _ in x),
        lambda x: "".join(_unicode_sup[str(_)] for _ in x),
        formula,
        prefixes,
        infixes,
        **kwargs
    )


def formula_to_html(formula, prefixes=None, infixes=None, **kwargs):
    """Convert formula string to html string representation

    Parameters
    ----------
    formula : str
        Chemical formula, e.g. 'H2O', 'Fe+3', 'Cl-'
    prefixes : dict
        Prefix transformations, default: greek letters and .
    infixes : dict
        Infix transformations, default: .
    suffixes : tuple of strings
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
    if infixes is None:
        infixes = _html_infix_mapping
    return _formula_to_format(
        lambda x: "<sub>%s</sub>" % x,
        lambda x: "<sup>%s</sup>" % x,
        formula,
        prefixes,
        infixes,
        **kwargs
    )
