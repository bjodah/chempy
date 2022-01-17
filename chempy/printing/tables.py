
class _RxnTable(object):

    _rsys_meth = None

    def __init__(self, idx_rxn_pairs, substances, colors=None, missing=None, missing_color='eee8aa', prefer_name=True):
        self.idx_rxn_pairs = idx_rxn_pairs
        self.substances = substances
        self.colors = colors or {}
        self.missing = missing or []
        self.missing_color = missing_color
        self.prefer_name = prefer_name

    @classmethod
    def from_ReactionSystem(cls, rsys, color_categories=True, **kwargs):
        idx_rxn_pairs, unconsidered_ri = getattr(rsys, cls._rsys_meth)()
        colors = rsys._category_colors() if color_categories else {}
        missing = [not rsys.substance_participation(sk) for sk in rsys.substances]
        return cls(idx_rxn_pairs, rsys.substances, colors=colors, missing=missing, **kwargs), unconsidered_ri

    def _repr_html_(self):
        from .web import css
        return css(self, substances=self.substances, colors=self.colors)

    def _cell_label_html(self, printer, ori_idx, rxn):
        """ Reaction formatting callback. (reaction index -> string) """
        pretty = rxn.unicode(self.substances, with_param=True, with_name=False)
        return '<a title="%d: %s">%s</a>' % (ori_idx, pretty, printer._print(
            (rxn.name or rxn.param) if self.prefer_name else rxn.param))

    @staticmethod
    def _cell_color(ri, ci, scale=1.0):
        if ci is None:
            return []
        else:
            colors = {
                (0, 0): (180,221,194),
                (0, 1): (180,214,221),
                (1, 0): (165,211,183),
                (1, 1): (165,193,201),
            }
            return ['style="background-color: #%x%x%x;"' % tuple(map(lambda c: min(255, int(c*scale)), colors[ri%2, ci%2]))]

    def _cell_html(self, printer, A, ri, ci=None):
        args = []
        if ci is not None and ri > ci:
            r = '-'
            args.extend(self._cell_color(ri, ci, scale=0.9))
        else:
            if ci is None:  # A is a vector
                c = A[ri]
                is_missing = self.missing[ri]
            else:  # A is a matrix
                c = A[ri][ci]
                is_missing = self.missing[ri] or self.missing[ci]

            if c is None:
                r = ''
            else:
                r = ', '.join(self._cell_label_html(printer, *r) for r in c)

            if is_missing:
                args.append('style="background-color: #%s;"' % self.missing_color)
            else:
                args.extend(self._cell_color(ri, ci, scale=1.1))

        return '<td %s>%s</td>' % (' '.join(args), r)


class UnimolecularTable(_RxnTable):
    """ Table of unimolecular reactions in a ReactionSystem

    Parameters
    ----------
    rsys : ReactionSystem
    sinks_sources_disjoint : tuple, None or True
        Colors sinks & sources. When ``True`` :meth:`sinks_sources_disjoint` is called.
    html_cell_label : Reaction formatting callback
        The function takes an integer, a Reaction instance and a dict of Substances as
        parameters and return a string.

    Returns
    -------
    string: html representation
    list: reactions not considered
    """

    _rsys_meth = '_unimolecular_reactions'

    def _html(self, printer, **kwargs):
        if 'substances' not in kwargs:
            kwargs['substances'] = self.substances
        ss = printer._get('substances', **kwargs)
        rows = '\n'.join('<tr><td>%s</td>%s</tr>' % (
            printer._print(s), self._cell_html(printer, self.idx_rxn_pairs, rowi)
        ) for rowi, s in enumerate(ss.values()))
        return '<table>%s</table>' % rows


class BimolecularTable(_RxnTable):
    """ Table of bimolecular reactions

    Parameters
    ----------
    idx_rxn_pairs : iterable of (int, Reaction) pairs
    substances : dict
        Mapping substance key to Substance instance.
    sinks_sources_disjoint : tuple, None or True
        Colors sinks & sources. When ``True`` :meth:`sinks_sources_disjoint` is called.

    Returns
    -------
    string: html representation
    list: reactions not considered
    """

    _rsys_meth = '_bimolecular_reactions'

    def _html(self, printer, **kwargs):
        if 'substances' not in kwargs:
            kwargs['substances'] = self.substances
        ss = printer._get('substances', **kwargs)
        header = '<th></th>' + ''.join('<th>%s</th>' % printer._print(s) for  # style="text-align: left; transform: rotate(-90deg);"
                                       s in ss.values())
        rows = ['<tr><td>%s</td>%s</tr>' % (
            printer._print(s), ''.join(self._cell_html(printer, self.idx_rxn_pairs, rowi, ci)
                                       for ci in range(len(ss)))
        ) for rowi, s in enumerate(ss.values())]
        return '<table>%s</table>' % '\n'.join([header, '\n'.join(rows)])  # cellpadding="0"
