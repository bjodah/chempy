# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from .util.pyutil import ChemPyDeprecationWarning

import warnings

from .electrolytes import (
    A, B, limiting_log_gamma, extended_log_gamma, davies_log_gamma,
    limiting_activity_product, extended_activity_product,
    davies_activity_product
)


warnings.warn("use .electrolytes instead of .debye_huckel", ChemPyDeprecationWarning)
