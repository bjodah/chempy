# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from itertools import chain
from .periodic import groups, symbols, names

_anions = {
    'F-': 'fluoride',
    'Cl-': 'chloride',
    'Br-': 'bromide',
    'I-': 'iodide',
    'OH-': 'hydroxide',
    'CN-': 'cyanide',
    'SCN-': 'thiocyanate',
    'SO4-2': 'sulphate',
    'HSO4-': 'hydrogensulphate',
    'PO4-3': 'phospahte',
    'HPO4-2': 'hydrogenphospahte',
    'H2PO4-': 'dihydrogenphospahte',
    'NO3-': 'nitrate',
    'NO2-': 'nitrite',
    'ClO-': 'hypochlorite',
    'ClO2-': 'chlorite',
    'ClO3-': 'chlorate',
    'ClO4-': 'perchlorate',
}
_alkali = [
    (symbols[n]+'+', names[n].lower()) for n in groups[1]
]
_alkaline_earth = [
    (symbols[n]+'+2', names[n].lower()) for n in groups[2]
]
_all_names = dict(chain(_alkali, _alkaline_earth, _anions.items()))


def name(ion):
    return _all_names[ion]
