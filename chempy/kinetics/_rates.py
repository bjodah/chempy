# -*- coding: utf-8 -*-
"""
Non-public API (classes in this module may change without notice).

The purpose here is to define conventions, e.g. lower-case string
 'temperature' is used, opposed to e.g. 'T', 'Temperature', etc.
"""
from __future__ import (absolute_import, division, print_function)

from ..util._expr import create_Poly, create_Piecewise

TPoly = create_Poly('temperature')
RTPoly = create_Poly('temperature', reciprocal=True)
Log10TPoly = create_Poly('log10_temperature')
ShiftedTPoly = create_Poly('temperature', shift='Tref', name='ShiftedTPoly')
ShiftedLog10TPoly = create_Poly('log10_temperature', shift='log10_Tref')
ShiftedRTPoly = create_Poly('temperature', shift='Tref', reciprocal=True)
TPiecewise = create_Piecewise('temperature', nan_fallback=True)
