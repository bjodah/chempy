# -*- coding: utf-8 -*-
import sys
from .string import StrPrinter
from .numbers import number_to_scientific_unicode


def _formula_fmt(s):
    def _script(s, chs):
        for i in range(10):
            s = s.replace(chr(ord("0") + i), chs[i])
        return s

    e = ""
    if "-" in s:
        s, e = s.split("-")
        e += u"⁻"
    elif "+" in s:
        s, e = s.split("+")
        e += u"⁺"

    return _script(s, u"₀₁₂₃₄₅₆₇₈₉") + _script(e, u"⁰¹²³⁴⁵⁶⁷⁸⁹")


class UnicodePrinter(StrPrinter):

    _default_settings = dict(
        StrPrinter._default_settings,
        repr_name="unicode",
        Equilibrium_arrow=u"⇌",
        Reaction_arrow=u"→",
        Reaction_formula_fmt=_formula_fmt,
        magnitude_fmt=number_to_scientific_unicode,
        unit_fmt=lambda dim: (
            dim.unicode
            if sys.version_info[0] > 2
            else dim.unicode.decode(encoding="utf-8")
        ),
    )
    _str = str if sys.version_info[0] > 2 else unicode  # noqa

    def _print_Substance(self, s, **kwargs):
        return s.unicode_name or s.name


def unicode_(obj, **settings):  # Python 2 keyword, hence the trailing '_'
    return UnicodePrinter(settings).doprint(obj)
