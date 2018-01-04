# -*- coding: utf-8 -*-
import sys
from .string import StrPrinter
from .numbers import number_to_scientific_unicode


class UnicodePrinter(StrPrinter):

    _default_settings = dict(
        StrPrinter._default_settings,
        repr_name='unicode',
        Equilibrium_arrow=u'⇌',
        Reaction_arrow=u'→',
        magnitude_fmt=number_to_scientific_unicode,
        unit_fmt=lambda dim: (
            dim.unicode if sys.version_info[0] > 2
            else dim.unicode.decode(encoding='utf-8')
        )
    )
    _str = str if sys.version_info[0] > 2 else unicode  # noqa

    def _print_Substance(self, s, **kwargs):
        return s.unicode_name or s.name


def unicode_(obj, **settings):  # Python 2 keyword, hence the trailing '_'
    return UnicodePrinter(settings).doprint(obj)
