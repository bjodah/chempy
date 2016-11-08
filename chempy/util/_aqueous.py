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
    'CO3-2': 'carbonate',
    'C2O4-2': 'oxalate',
    'HCO3-': 'hydrogencarbonate',
    'NO3-': 'nitrate',
    'NO2-': 'nitrite',
    'PO4-3': 'phospahte',
    'HPO4-2': 'hydrogenphospahte',
    'H2PO4-': 'dihydrogenphospahte',
    'P-3': 'phosphide',
    'SO4-2': 'sulphate',
    'HSO4-': 'hydrogensulphate',
    'SO3-2': 'sulphite',
    'HSO3-': 'hydrogensulphite',
    'S-2': 'sulfide',
    'ClO-': 'hypochlorite',
    'ClO2-': 'chlorite',
    'ClO3-': 'chlorate',
    'ClO4-': 'perchlorate',
    'CrO4-2': 'chromate(VI)',
    'Cr2O7-2': 'dichromate(VI)',
    'MnO4-2': 'manganate(VI)',
    'MnO4-': 'permanganate(VII)',
    'FeO4-2': 'ferrate(VI)',
    'OsO4-2': 'osmate(VI)',
    'Bo3-3': 'borate',
    'BiO3-': 'bismuthate(V)',
}
_cations = {
    'H3O+': 'hydronium',
}
_cation_oxidation_states = {  # This needs to be reviewed, just from the top of my head
    'Cr': (2, 3),
    'Fe': (2, 3),
    'Mn': (2,),
    'Co': (2, 3),
    'Ni': (2, 3),
    'Cu': (1, 2, 3),
    'Ag': (1, 2),
    'Au': (3,),
    'Zn': (2,),
    'Cd': (2,),
    'Hg': (1, 2),  # Tricky: Hg2+2
    'Al': (3,),
    'Ga': (3,),
    'In': (3,),
    'Tl': (1, 3),
    'Sn': (2, 4),
    'Pb': (2, 4),
    'Bi': (3,),
    'Sb': (3,),
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


def ions_from_formula(formula):
    """
    This will be working examples eventually:

    #>>> ions_from_formula('NaCl') == {'Na+': 1, 'Cl-': 1}
    #True
    #>>> ions_from_formula('Fe(NO3)3') == {'Fe+3': 1, 'NO3-': 3}
    #True
    #>>> ions_from_formula('FeSO4') == {'Fe+2': 1, 'SO4-2': 1}
    #True
    #>>> ions_from_formula('(NH4)3PO4') == {'NH4+': 3, 'PO4-3': 1}
    #True
    #>>> ions_from_formula('KAl(SO4)2.11H2O') == {'K+': 1, 'Al+3': 1, 'SO4-2': 2}
    #True

    """
    pass
