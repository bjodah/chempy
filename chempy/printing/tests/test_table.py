# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict

from chempy import Substance
from chempy.util.testing import requires
from chempy.util.parsing import parsing_library
from ..table import as_per_substance_html_table
from ..web import html


@requires(parsing_library)
def test_as_per_substance_html_table():
    substances = OrderedDict([(k, Substance.from_formula(k)) for k in 'H OH'.split()])
    assert html(as_per_substance_html_table([2, 3], substances)).count('<tr>') == 3
    assert html(as_per_substance_html_table({'H': 2, 'OH': 3}, substances)).count('<tr>') == 3
