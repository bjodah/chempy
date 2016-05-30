# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from .util.parsing import number_to_scientific_html


class Table(object):
    def __init__(self, rows, headers=None):
        self.rows, self.headers = rows, headers

    def html(self):
        rows = ['\n'.join('<th>%s</th>' % _ for _ in self.headers)]
        rows += ['\n'.join('<td>%s</td>' % _ for _ in r) for r in self.rows]
        return '<table>%s</table>' % '\n'.join(['<tr>%s</tr>' % r for r in rows])

    def _repr_html_(self):
        return self.html()


def as_per_substance_html_table(cont, substances, header='Concentration'):
    rows = [(v.html_name, number_to_scientific_html(cont[k])) for k, v in substances.items()]
    return Table(rows, ['Substance', header])
