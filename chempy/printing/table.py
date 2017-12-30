# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict
from ..chemistry import Substance
from .numbers import number_to_scientific_html


class Table(object):
    def __init__(self, rows, headers=None):
        self.rows, self.headers = rows, headers

    def _html(self, printer, **kwargs):
        def map_fmt(cont, fmt, joiner='\n'):
            return joiner.join(map(lambda x: fmt % printer._print(x, **kwargs), cont))
        rows = (
            [map_fmt(self.headers, '<th>%s</th>')] +
            [map_fmt(row, '<td>%s</td>') for row in self.rows]
        )
        return '<table>%s</table>' % map_fmt(rows, '<tr>%s</tr>')

    def _repr_html_(self):
        from .web import html
        return html(self)


def as_per_substance_html_table(cont, substances=None, header=None,
                                substance_factory=Substance.from_formula):
    """ """
    if substances is None:
        substances = OrderedDict([(k, substance_factory(k)) for k in cont])

    def _elem(k):
        try:
            return cont[k]
        except (IndexError, TypeError):
            return cont[list(substances.keys()).index(k)]
    rows = [(v.html_name, number_to_scientific_html(_elem(k))) for k, v in substances.items()]
    return Table(rows, ['Substance', header or ''])
