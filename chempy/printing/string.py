from operator import itemgetter
from .printer import Printer


class StrPrinter(Printer):
    _default_settings = dict(
        Printer._default_settings,
        repr_name='string',
        Equilibrium_arrow='=',
        Reaction_arrow='->',
    )

    def _Reaction_parts(self, rxn, **kwargs):
        str_ = self._str
        coeff_fmt = self._get('Reaction_coeff_fmt', **kwargs)
        formula_fmt = self._get('Reaction_formula_fmt', **kwargs)
        substances = self._get('substances', **kwargs) or {}
        nullstr, space = str_(''), str_(' ')
        reac, prod, i_reac, i_prod = [[
            (
                ((coeff_fmt(v)+space) if v != 1 else nullstr) +
                formula_fmt(self._print(substances.get(k, k)))
            ) for k, v in filter(itemgetter(1), d.items())
        ] for d in (rxn.reac, rxn.prod, rxn.inact_reac, rxn.inact_prod)]
        r_str = str_(" + ").join(reac)
        ir_str = (str_(' + ( ') + str_(" + ").join(i_reac) + str_(')')
                  if len(i_reac) > 0 else nullstr)
        arrow_str = self._get('%s_arrow' % rxn.__class__.__name__, **kwargs)
        p_str = str_(" + ").join(prod)
        ip_str = (str_(' + ( ') + str_(" + ").join(i_prod) + str_(')')
                  if len(i_prod) > 0 else nullstr)
        return r_str, ir_str, arrow_str, p_str, ip_str

    def _Reaction_str(self, rxn, **kwargs):
        fmtstr = self._str("{}{}%s{}%s{}{}") % self._get('Reaction_around_arrow', **kwargs)
        return fmtstr.format(*self._Reaction_parts(rxn, **kwargs))

    def _Reaction_param_str(self, rxn, **kwargs):
        mag_fmt = self._get('magnitude_fmt', **kwargs)
        unit_fmt = self._get('unit_fmt', **kwargs)
        try:
            magnitude_str = mag_fmt(rxn.param.magnitude)
            unit_str = unit_fmt(rxn.param.dimensionality)
        except AttributeError:
            try:
                return mag_fmt(rxn.param)
            except TypeError:
                return str(rxn.param)
        else:
            return magnitude_str + self._str(' ') + unit_str

    def _print_Reaction(self, rxn, **kwargs):
        res = self._Reaction_str(rxn, **kwargs)
        if self._get('with_param', **kwargs) and rxn.param is not None:
            res += self._get('Reaction_param_separator', **kwargs)
            try:
                res += getattr(rxn.param, self._get('repr_name', **kwargs))(
                    self._get('magnitude_fmt', **kwargs))
            except AttributeError:
                res += self._Reaction_param_str(rxn, **kwargs)
        return res

    def _print_ReactionSystem(self, rsys, **kwargs):
        header = (rsys.name + '\n') if rsys.name else ''
        return header + '\n'.join(map(self._print, rsys.rxns)) + '\n'


def str_(obj, **settings):  # Python keyword, hence the trailing '_'
    return StrPrinter(settings).doprint(obj)
