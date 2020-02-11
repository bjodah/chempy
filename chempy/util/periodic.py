# -*- coding: utf-8 -*-
"""Contains atomic data from the periodic table.

Contains atomic data from the periodic table and from other sources,
such as the aufbau rules.  Current data are either from the original
data (version 0.7.9, marked '*_old') or from NIST ('*_NIST').  The old
atomic mass data were licensed under the CC-SA license and were found
at
https://en.wikipedia.org/w/index.php?title=List_of_elements&oldid=700476748.
All updated periodic table data marked NIST is from the NIST Periodic
Table, July 2019 version from http://www.nist.gov/.
"""
import importlib.resources
import re
from . import data

# Some preliminary work on useful(?) classes.


class Atom:
    """Manipulate data concerning an atom."""

    def __init__(self, atomic_number, NIST=False):
        """Initialize new Atom() objects.

        Parameters
        ----------
        atomic_number: int
            The atomic number of the atom/ion.
        NIST: bool, optional
            Use NIST data, or not (default is False, use original data).
        """
        self.atomic_number = atomic_number
        self.electrons = atomic_number
        self.charge = 0
        self.symbol = symbols[atomic_number - 1]
        self.name = names[atomic_number - 1]
        self.atomic_mass = relative_atomic_masses[atomic_number - 1]
        self.ionization_energy = ''
        self.ec_naive = naive_electron_configuration_builder(atomic_number, 0)
        self.ec = ''

        # Option to use old data or current NIST data.
        self.NIST = NIST

        return


class AtomicIon(Atom):
    """Represents data concerning an atomic ion.  Subclasses Atom()."""

    def __init__(self, atomic_number, charge):
        """Initialize new AtomicIon() objects."""
        super().__init__(atomic_number)
        self.charge = charge
        self.electrons = atomic_number - charge

        return


class Isotope(Atom):
    """Represents data concerning an isotope.  Subclasses Atom()."""

    def __init__(self, atomic_number, neutrons):
        """Initialize new Isotope() objects."""
        super().__init__(atomic_number)
        self.neutrons = neutrons

        return

    def getMassNumber(self):
        """Calculate the mass number of the isotope."""
        return self.atomic_number + self.neutrons


def naive_electron_configuration_builder(
        atomic_number,
        charge,
        noble_gas_shorthand=False):
    """Build the naive ground-state electron configuration.

    Build the naive (aufbau rules) ground-state electron configuration
    of a given atom or ion.

    Parameters
    ----------
    atomic_number: int
        The atomic number of the atom/ion.
    charge: int
        The charge of the ion, or zero for an atom.
    noble_gas_shorthand: bool
        Output the Noble gas shorthand.  Default is False.

    Returns
    -------
    ec: string
        A naive, aufbau electron configuration like you do in high school.

    Notes
    -----
    This result often disagrees with experimental data, especially of d-
    and f-block elements.

    Examples
    --------
    >>> naive_electron_configuration_builder(2, 0)
    '1s2'
    >>> naive_electron_configuration_builder(9, 0)
    '1s2 2s2 2p5'
    >>> naive_electron_configuration_builder(9, -1)
    '1s2 2s2 2p6'
    >>> naive_electron_configuration_builder(9, -1, True)
    '[He] 2s2 2p6'
    """
    electrons = atomic_number - charge

    return ec


# Load atomic data.
#
# All data is stored in 118 line files in ./data and imported with
# importlib.resources.  Original data is in '-old' files and recent
# NIST data is in '-NIST' files.  Contents:
#
# symbols-*:  element symbols
# names-*:  element names
# masses-*:  atomic masses
# ionization-energies-*:  ionization energies (NIST only)
#
# Currently, the arrays are indexed from zero and not atomic number.
# Using None as the zeroth element to allow atomic number indexing may
# be beneficial.

atomic_data = {}
sources = (
    'symbols_old',
    'symbols_NIST',
    'names_old',
    'names_NIST',
    'masses_old',
    'masses_NIST',
    'ionization_energies_NIST',
    )

for source in sources:
    atomic_data[source] = []
    with importlib.resources.open_text(data, source) as file:
        if re.match(r'^masses_', source):
            for line in file:
                line = line.strip()
                if line.startswith('[') and line.endswith(']'):
                    atomic_data[source].append(float(line[1:-1]))
                elif '(' in line:
                    atomic_data[source].append(float(line.split('(')[0]))
                else:
                    atomic_data[source].append(float(line))
        else:
            for line in file:
                atomic_data[source].append(line.strip())

symbols_NIST = atomic_data['symbols_NIST']
symbols_old = atomic_data['symbols_old']
names_NIST = atomic_data['names_NIST']
names_old = atomic_data['names_old']
relative_atomic_masses_NIST = atomic_data['masses_NIST']
relative_atomic_masses_old = atomic_data['masses_old']
ionization_energies = atomic_data['ionization_energies_NIST']

period_lengths = (2, 8, 8, 18, 18, 32, 32)
accum_period_lengths = (2, 10, 18, 36, 54, 86, 118)

# groups: dict with the group number as key and an array of atomic
# number of the group memmbers as value.  Includes the s- and p-blocks
# only.
# Icosagens, crystallogens, pnictogens, chalcogens, halogens.
groups = {g: tuple(x - 18 + g for x in accum_period_lengths[1:])
          for g in range(13, 18)}
# Alkali metals.
groups[1] = (1,) + tuple(x + 1 for x in accum_period_lengths[:-1])
# Alkaline earth metals.
groups[2] = tuple(x + 2 for x in accum_period_lengths[:-1])
# Noble gases.
groups[18] = accum_period_lengths

# Boolean to control use of NIST data or old data.
# Default is False, use old data.
use_NIST = False


def set_NIST(value):
    """Set the data source of chempy.util.periodic.

    Set the behavior of chempy.util.periodic to use either the
    original data (False) or updated NIST data (True).

    Parameters
    ----------
    value: bool
        Use NIST (True) or not (False).

    Returns
    -------
    use_NIST: bool
        The current value of use_NIST.

    Notes
    -----
    The function coerces whatever you send it to type bool.


    Examples
    --------
    >>> set_NIST(True)
    'True'
    >>> set_NIST(False)
    'False'
    """
    use_NIST = bool(value)

    return use_NIST


# Boolean to control whether we have initialized or not.
_initialized = False

symbols = ''
names = ''
lower_names = ''
relative_atomic_masses = ''
electron_mass = ''


def _initialize():

    if not _initialized:
        if use_NIST:
            toNIST()
        else:
            toOld()

    return _initialized


def toNIST():

    global _initialized
    global symbols
    global names
    global lower_names
    global relative_atomic_masses
    global electron_mass

    symbols = symbols_NIST
    names = names_NIST
    electron_mass = 6.02214076e23 * 9.10938370e-31 * 1000
    lower_names = tuple(n.lower().lstrip('(').rstrip(')') for n in names)
    relative_atomic_masses = relative_atomic_masses_NIST
    _initialized = True

    return _initialized


def toOld():

    global _initialized
    global symbols
    global names
    global lower_names
    global relative_atomic_masses
    global electron_mass

    symbols = symbols_old
    names = names_old
    electron_mass = 5.489e-4
    lower_names = tuple(n.lower().lstrip('(').rstrip(')') for n in names)
    relative_atomic_masses = relative_atomic_masses_old
    _initialized = True

    return _initialized


def atomic_number(name):
    """Find the atomic number of the given symbol or name.

    Parameters
    ----------
    name: string
        The name or symbol of the element.

    Returns
    -------
    int
        The atomic number of the element.

    Examples
    --------
    >>> atomic_number('He')
    '2'
    >>> atomic_number('Lithium')
    '3'
    """
    if not _initialized:
        _initialize()
    # If name is actually the atomic number.
    if isinstance(name, int):
        if 1 <= name <= 118:
            return name
    # Try using the symbol.
    try:
        return symbols.index(name) + 1
    # If not the symbol, try the name.
    except ValueError:
        return lower_names.index(name.lower()) + 1


def mass_from_composition(composition):
    """Calculate the molecular mass from atomic masses.

    Parameters
    ----------
    composition: dict
        Dictionary mapping int (atomic number) to int (coefficient)

    Returns
    -------
    mass: float
        molecular mass in atomic mass units (g/mol)

    Notes
    -----
    Atomic number 0 denotes charge or "net electron deficiency"

    Examples
    --------
    >>> '%.2f' % mass_from_composition({0: -1, 1: 1, 8: 1})
    '17.01'
    """
    if not _initialized:
        _initialize()
    mass = 0.0
    for k, v in composition.items():
        # For electrons.
        if k == 0:
            mass -= v * electron_mass
        # For elements.
        else:
            mass += v * relative_atomic_masses[k - 1]
    return mass


_initialize()
