from .numbers import number_to_scientific_html
from .string import StrPrinter


def _html_clsname(key):
    return "chempy_" + key.replace(
        '+', 'plus').replace(
        '-', 'minus').replace(
        '(', 'leftparen').replace(
        ')', 'rightparen')

_html_semicolon = '&#59; '

class HTMLPrinter(StrPrinter):

    printmethod_attr = '_html'
    _default_settings = dict(
        StrPrinter._default_settings,
        repr_name='html',
        Equilibrium_arrow='&harr;',
        Reaction_arrow='&rarr;',
        Reaction_param_separator=_html_semicolon,
        magnitude_fmt=number_to_scientific_html
    )

    def _print_Substance(self, s, **kwargs):
        return s.html_name or s.name

    def _print_ReactionSystem(self, rsys, **kwargs):
        return super(HTMLPrinter, self)._print_ReactionSystem(rsys, **kwargs).replace('\n', '<br>\n')

def html(obj, **settings):
    return HTMLPrinter(settings).doprint(obj)


class CSSPrinter(HTMLPrinter):
    def _print_Substance(self, s, **kwargs):
        key = s.name
        name = s.html_name or key
        common_sty = 'border-radius: 5pt; padding: 0pt 3pt 0pt 3pt;'
        colors = self._get('colors', **kwargs)
        if key in colors:
            style = 'background-color:#%s; border: 1px solid #%s; %s' % (colors[key] + (common_sty,))
        else:
            style = 'style="%s"' % common_sty
        fmt = '<span class="%s" style="%s">%s</span>'
        return fmt % (_html_clsname(key), style, name)

    def _print_ReactionSystem(self, rsys, **kwargs):
        sep = '</td><td style="text-align:left;">&nbsp;'
        around = '</td><td style="text-align:center;">', '</td><td style="text-align:left;">'
        # cf. https://github.com/jupyter/notebook/issues/2160#issuecomment-352216152
        rxns = '</td></tr>\n<tr><td style="text-align:right;">'.join(map(lambda r: self._print(
            r, Reaction_param_separator=sep, Reaction_around_arrow=around), rsys.rxns))
        tab = '<table class="chempy_table"><tr><td style="text-align:right;">%s</td></tr></table>'
        return tab % rxns


def css(obj, **settings):
    return CSSPrinter(settings).doprint(obj)


def bimolecular_table(rsys, sinks_sources_disjoint=None, cell_label=None):
    if sinks_sources_disjoint is True:
        sinks_sources_disjoint = rsys.sinks_sources_disjoint()

    if sinks_sources_disjoint:
        _str_formula = _str_formula_factory(*sinks_sources_disjoint[:2])
    else:
        def _str_formula(s, k):
            return s

def unimolecular_table(rsys, sinks_sources_disjoint=None, cell_label=None):
    if sinks_sources_disjoint is True:
        sinks_sources_disjoint = rsys.sinks_sources_disjoint()

    if sinks_sources_disjoint:
        _str_formula = _str_formula_factory(*sinks_sources_disjoint[:2])
    else:
        def _str_formula(s, k):
            return s

    A, unconsidered = rsys._unimolecular_reactions()
