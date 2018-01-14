
class _RxnTable(object):

    _rsys_meth = None

    def __init__(self, idx_rxn_pairs, substances, colors=None, missing=None, missing_color='eee8aa'):
        self.idx_rxn_pairs = idx_rxn_pairs
        self.substances = substances
        self.colors = colors or {}
        self.missing = missing or []
        self.missing_color = missing_color

    @classmethod
    def from_ReactionSystem(cls, rsys, color_categories=True):
        idx_rxn_pairs, unconsidered_ri = getattr(rsys, cls._rsys_meth)()
        colors = rsys._category_colors() if color_categories else {}
        missing = [not rsys.substance_participation(sk) for sk in rsys.substances]
        return cls(idx_rxn_pairs, rsys.substances, colors=colors, missing=missing), unconsidered_ri

    def _repr_html_(self):
        from .web import css
        return css(self, substances=self.substances, colors=self.colors)

    def _cell_label_html(self, ori_idx, rxn):
        """ Reaction formatting callback. (reaction index -> string) """
        pretty = rxn.unicode(self.substances, with_param=True)
        return '<a title="%d: %s">%s</a>' % (ori_idx, pretty, rxn.name or rxn.param)

    def _cell_html(self, A, ri, ci=None):
        args = []
        if ci is not None and ri > ci:
            r = '-'
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
                r = ', '.join(self._cell_label_html(*r) for r in c)

            if is_missing:
                args.append('style="background-color: #%s;"' % self.missing_color)

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
            printer._print(s), self._cell_html(self.idx_rxn_pairs, rowi)
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
        header = '<th></th>' + ''.join('<th>%s</th>' % printer._print(s) for
                                       s in ss.values())
        rows = ['<tr><td>%s</td>%s</tr>' % (
            printer._print(s), ''.join(self._cell_html(self.idx_rxn_pairs, rowi, ci)
                                       for ci in range(len(ss)))
        ) for rowi, s in enumerate(ss.values())]
        return '<table>%s</table>' % '\n'.join([header, '\n'.join(rows)])
