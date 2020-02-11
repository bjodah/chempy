# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from itertools import chain
from .periodic import groups, symbols, names

_anions = {
    # Atomic anions.
    'N-3': 'nitride',
    'P-3': 'phosphide',
    'As-3': 'arsenide',
    'O-2': 'oxide',
    'S-2': 'sulfide',
    'Se-2': 'selenide',
    'Te-2': 'telluride',
    'F-': 'fluoride',
    'Cl-': 'chloride',
    'Br-': 'bromide',
    'I-': 'iodide',
    'BO3-3': 'borate',
    'BiO3-': 'bismuthate(V)',
    'C2O4-2': 'oxalate',
    'CN-': 'cyanide',
    'CO3-2': 'carbonate',
    'HCO3-': 'hydrogen carbonate',
    'ClO-': 'hypochlorite',
    'ClO2-': 'chlorite',
    'ClO3-': 'chlorate',
    'ClO4-': 'perchlorate',
    'Cr2O7-2': 'dichromate(VI)',
    'CrO4-2': 'chromate(VI)',
    'FeO4-2': 'ferrate(VI)',
    'MnO4-': 'permanganate(VII)',
    'MnO4-2': 'manganate(VI)',
    'NO2-': 'nitrite',
    'NO3-': 'nitrate',
    'OH-': 'hydroxide',
    'OsO4-2': 'osmate(VI)',
    'PO4-3': 'phosphate',
    'HPO4-2': 'hydrogen phosphate',
    'H2PO4-': 'dihydrogen phosphate',
    'SCN-': 'thiocyanate',
    'SO3-2': 'sulphite',
    'HSO3-': 'hydrogen sulphite',
    'SO4-2': 'sulphate',
    'HSO4-': 'hydrogen sulphate',
}

_cations = {
    'H3O+': 'hydronium',
    'NH4+': 'ammonium',
    'C7H7+': 'tropylium',
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
_all_names = dict(chain(_alkali, _alkaline_earth,
                        _cations.items(), _anions.items()))


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
