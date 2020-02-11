# -*- coding: utf-8 -*-
"""Provides function and data for ionic species."""
from __future__ import (absolute_import, division, print_function)

from itertools import chain
from .periodic import groups, symbols, names

# Should be able to infer transition metal oxidation state from
# charge.
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
    'CH3COO-': 'acetate',
    'C6H5COO-': 'benzoate',
    'C4H406-2': 'tartrate',
    'BF4-': 'tetrafluoroborate',
    'S4O6-2': 'tetrathionate',
    'AsO4-3': 'arsenate',
    'AsO3-3': 'arsenite',
    'N3-': 'azide',
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
    'BrO3-': 'bromate',
    'IO3-': 'iodate',
    'SiO3-2': 'silicate',
    'Cr2O7-2': 'dichromate(VI)',
    'CrO4-2': 'chromate(VI)',
    'FeO4-2': 'ferrate(VI)',
    'MnO4-': 'permanganate(VII)',
    'MnO4-2': 'manganate(VI)',
    'NO2-': 'nitrite',
    'NO3-': 'nitrate',
    'OH-': 'hydroxide',
    'O2-': 'peroxide',
    'OsO4-2': 'osmate(VI)',
    'PO4-3': 'phosphate',
    'HPO4-2': 'hydrogen phosphate',
    'H2PO4-': 'dihydrogen phosphate',
    'PO3-3': 'phosphite',
    'HPO3-2': 'hydrogen phosphite',
    'H2PO3-': 'dihydrogen phosphite',
    'SCN-': 'thiocyanate',
    'OCN-': 'cyanate',
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

# Fixed.  Common oxidation states for the transition metals.  This is
# not the complete list of oxidation states.
#
# From: Chemistry of the Elements; Greenwood, NN, Earnshaw, A; second
# edition; Butterworth-Heinemann, 1997, page 28.
_cation_oxidation_states = {
    'Al': (3,),
    'Sc': (3,),
    'Ti': (4,),
    'V': (5,),
    'Mn': (2, 4, 7),
    'Cr': (3, 6),
    'Fe': (2, 3),
    'Co': (2, 3),
    'Ni': (2,),
    'Cu': (2,),
    'Zn': (2,),
    'Ga': (1, 2),
    'Ge': (1, 3),
    'Y': (3,),
    'Zr': (4,),
    'Nb': (5,),
    'Mo': (3, 6),
    'Tc': (4, 7),
    'Ru': (3, 4),
    'Rh': (3,),
    'Pd': (2, 4),
    'Ag': (1,),
    'Cd': (2,),
    'In': (3,),
    'Sn': (2, 4),
    'Sb': (3, 5),
    'Hf': (4,),
    'Ta': (5,),
    'W': (4, 6),
    'Re': (4,),
    'Os': (4,),
    'Ir': (3, 4),
    'Pt': (2, 4),
    'Au': (3,),
    'Hg': (1, 2),  # Tricky: Hg2+2
    'Tl': (1, 3),
    'Pb': (2, 4),
    'Bi': (3,),
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
    """Name all the available ions."""
    return _all_names[ion]


def ions_from_formula(formula):
    """Parse ions from formula.

    #>>> ions_from_formula('NaCl') == {'Na+': 1, 'Cl-': 1}
    #True
    #>>> ions_from_formula('Fe(NO3)3') == {'Fe+3': 1, 'NO3-': 3}
    #True
    #>>> ions_from_formula('FeSO4') == {'Fe+2': 1, 'SO4-2': 1}
    #True
    #>>> ions_from_formula('(NH4)3PO4') == {'NH4+': 3, 'PO4-3': 1}
    #True
    #>>> ions_from_formula('KAl(SO4)2.11H2O')
    #...     == {'K+': 1, 'Al+3': 1, 'SO4-2': 2}
    #True

    """
    pass
