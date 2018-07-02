from .numbers import number_to_scientific_latex
from .string import StrPrinter
from ..units import _latex_from_dimensionality


class LatexPrinter(StrPrinter):

    _default_settings = dict(
        StrPrinter._default_settings,
        repr_name='latex',
        Equilibrium_arrow=r'\rightleftharpoons',
        Reaction_arrow=r'\rightarrow',
        magnitude_fmt=number_to_scientific_latex,
        unit_fmt=_latex_from_dimensionality
    )

    def _print_Substance(self, substance, **kwargs):
        return substance.latex_name or substance.name


def latex(obj, **settings):
    return LatexPrinter(settings).doprint(obj)
