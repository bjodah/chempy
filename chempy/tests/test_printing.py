# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict
from .. import Substance
from ..printing import as_per_substance_html_table
from ..util.testing import requires
from ..util.parsing import parsing_library


@requires(parsing_library)
def test_as_per_substance_html_table():
    substances = OrderedDict([(k, Substance.from_formula(k)) for k in 'H OH'.split()])
    assert as_per_substance_html_table([2, 3], substances).html().count('<tr>') == 3
    assert as_per_substance_html_table({'H': 2, 'OH': 3}, substances).html().count('<tr>') == 3
