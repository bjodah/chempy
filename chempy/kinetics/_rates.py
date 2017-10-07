# -*- coding: utf-8 -*-
"""
Non-public API (classes in this module may change without notice).

The purpose here is to define conventions, e.g. lower-case string
 'temperature' is used, opposed to e.g. 'T', 'Temperature', etc.
"""
from __future__ import (absolute_import, division, print_function)

from ..util._expr import create_Poly, create_Piecewise

TPoly = create_Poly('temperature')
RTPoly = create_Poly('temperature', reciprocal=True, args_dimensionality=lambda self, **kwargs: tuple(
    [{'temperature': i} for i in range(len(self.args))]))
Log10TPoly = create_Poly('log10_temperature', args_dimensionality=lambda self, **kwargs: tuple(
    [{'temperature': -i} for i in range(self.nargs)]))
ShiftedTPoly = create_Poly('temperature', shift='Tref', name='ShiftedTPoly',
                           args_dimensionality=lambda self, **kwargs: (({'temperature': 1}), tuple(
                               [{'temperature': -i} for i in range(self.nargs)])))
ShiftedLog10TPoly = create_Poly('log10_temperature', shift='log10_Tref',
                                args_dimensionality=lambda self, **kwargs: (({'temperature': 1}), tuple(
                                    [{'temperature': -i} for i in range(self.nargs)])))
ShiftedRTPoly = create_Poly('temperature', shift='Tref', reciprocal=True,
                            args_dimensionality=lambda self, **kwargs: (({'temperature': 1}), tuple(
                                [{'temperature': i} for i in range(self.nargs)])))
TPiecewise = create_Piecewise(
    'temperature', nan_fallback=True, args_dimensionality=lambda self, **kwargs: tuple([
        {'temperature': 1} if idx % 2 == 0 else self.args[idx].args_dimensionality()
        for idx in range(self.nargs)]))
