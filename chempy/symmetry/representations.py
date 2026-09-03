# -*- coding: utf-8 -*-

"""
Contains chemical group theory functions for calculating numbers and type of
irreducible representations in a reducible representation and number of IR
and Raman active modes. Below uses the inverse matrix method of solving for
irreducible representations described in J. Chem. Educ. 2009, 86, 251-253.
https://doi.org/10.1021/ed086p251.
"""

import numpy as np
from functools import wraps
from tabulate import tabulate

# data tables in character_tables.py file
from .tables import (
    tables,
    rot_trans_modes,
    IR_active,
    Raman_active,
    masks,
    atom_contribution,
    mulliken,
    headers,
    row_coeffs,
    column_coeffs
)


def print_header(group):
    """
    Print the header for a character table given the point group.

    Parameters
    ----------
    group : str
        Point group Schoenflies notation (e.g., 'C2v').  This is
        case-insensitive.

    Returns
    -------
    Prints character table header indicating order of symmetry operations.
    """
    symbols = headers[group.lower()]
    numbers = column_coeffs[group.lower()]

    header = []
    for i in range(len(numbers)):
        if numbers[i] != 1:
            header.append(str(numbers[i]) + symbols[i])
        else:
            header.append(symbols[i])

    print(*header)


def print_point_groups():
    """
    Print supported point group Schoenflies notations.

    Parameters
    ----------
    None

    Returns
    -------
    Prints Schoenflies notations for character tables available.

    """

    print(*[x.capitalize() for x in tables.keys()])


def print_mulliken(group):
    """
    Print Mulliken symbols of irreducible representation in order.

    Parameters
    ----------
    group : str
        Point group Schoenflies notation (e.g., 'C2v').  This is
        case-insensitive.

    Returns
    -------
    Prints Mulliken notations for irreducible representations

    """
    print(*mulliken[group.lower()])


def print_table(group):
    """
    Print character table for given point group.

    Parameters
    ----------
    group : str
        Point group Schoenflies notation (e.g., 'C2v').  This is
        case-insensitive.

    Returns
    -------
    Prints the character table for the point group.

    """
    group = group.lower()
    symbols = headers[group]
    mull_symbols = mulliken[group]
    numbers = column_coeffs[group.lower()]

    # generate list of symmetry symbols with coefficients
    header = []
    for i in range(len(numbers)):
        if numbers[i] != 1:
            header.append(str(numbers[i]) + symbols[i])
        else:
            header.append(symbols[i])

    # generate list of Mulliken symbols including duplicates if imaginary rep
    if group.lower() in row_coeffs.keys():
        rows = []
        for i in range(len(mull_symbols)):
            rows.extend([mull_symbols[i]] * row_coeffs[group][i])
    else:
        rows = mull_symbols

    table = [[group.capitalize(), *header]]
    if group == 'c1':
        table.append([rows[0], 1])
    else:
        for i_row in range(len(rows)):
            table.append([rows[i_row], *tables[group][i_row]])

    print(tabulate(table, tablefmt='rounded_grid'))


@np.vectorize
def _sympy_to_num(sp):
    """
    Convert array with sympy objects to array of floats and ints.

    Convert numpy array of sympy objects from table dictionary in tables.py
    to numpy array of floats and ints for use in calculations. Three types of
    values occur in the tables: ints, real sympy values, and imaginary sympy
    values.

    Parameters
    ----------
    sp : numpy array with sympy objects
        numpy array with sympy objects from table dictionary in tables.py

    Returns
    -------
    numpy array of floats and ints.

    """
    return complex(sp)


class Reducible:
    """Reducible representation object.

    Reducible representation object for calculating number and type of
    irreducible representations and calculating IR and Raman active vibrational
    modes.
    """

    def __init__(self, gamma, group, all_motion=False):
        """
        Initialize Reducible representation object.

        Parameters
        ----------
        gamma : List, tuple, or array
            Reducible representation.
        group : str
            Point group Schoenflies notation (e.g., 'C2v').  This is
            case-insensitive.
        all_motion : bool, optional
            False if reducible representation is only vibrational modes and
            True if the reducible describes all motions (rotational,
            translational, and vibrational). The default is False.

        Returns
        -------
        None.
        """
        if np.any(np.mod(gamma, 1) != 0):
            raise ValueError('Invalid representation - must be whole numbers.')
        elif group.lower() not in tables.keys():
            raise ValueError('Invalid point group.')
        elif np.array(gamma).size != tables[group.lower()].shape[0]:
            raise ValueError(f'Invalid representation size for {group}'
                             ' point group.')

        self.group = group.lower()
        self.gamma = np.atleast_1d(gamma)
        self.all_motion = all_motion

    def _return_dict(func):
        """
        Return results as a dictionary.

        Return a list or array as a dictionary with Mulliken symbols as
        the keys.

        Returns
        -------
        Dictionary.

        """
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            if kwargs.get('to_dict'):
                keys = mulliken[self.group.lower()]
                values = func(self, *args, **kwargs)
                return dict(zip(keys, values.tolist(), strict=True))
            else:
                return func(self, *args, **kwargs)
        return wrapper

    @_return_dict
    def decomp(self, to_dict=False):
        """
        Decompose reducible representation into number of irreducibles.

        Decompose a reducible representation for a specific point group and
        return the number of each irreducible representation in the reducible.

        The order of irreducibles can be looked up by print_mulliken()
        or print_table().

        Parameters
        ----------
        to_dict : bool
            True causes the function to return a dictionary with Mulliken
            symbols as the keys.

        Returns
        -------
        NumPy array with the number of each irreducible representation in the
        provided reducible representation. Use print_mulliken() for listing of
        irreducible representations and order. If to_dict=True, the function
        returns a dictionary with the Mulliken symbols as keys.

        Examples
        --------
        >>> rep = Reducible([9, -1, 3, 1], 'c2v', all_motion=True)
        >>> rep.decomp()
        array([3, 1, 3, 2])
        >>> rep = Reducible([15, 0, 0, 7, -2, -2], 'C3h', all_motion=False)
        >>> rep.decomp()
        array([3, 4, 2, 1])
        >>> rep.decomp(to_dict=True)
        {"A'": 3, "E'": 4, 'A"': 2, 'E"': 1}
        """
        table = _sympy_to_num(tables[self.group])
        gamma = np.array(self.gamma)

        if self.group == 'c1':
            n_i = self.gamma
        else:
            mask = np.array(masks[self.group], dtype=bool)
            # mask removes complex conjugate to avoid "doubling problem"
            n_i = gamma.dot(np.linalg.inv(table)).real[mask]

        if np.any(np.abs(n_i - np.rint(n_i)) > 0.02):
            raise ValueError('Invalid reducible representation. Does not'
                             ' decompose into irreducibles for this group.')

        return np.rint(n_i).astype(int)

    @_return_dict
    def vibe_modes(self, to_dict=False):
        """Return vibrational modes.

        Return the number of vibrational modes after rotation and translation
        are subtracted out.

        Parameters
        ----------
        to_dict : bool
            True causes the function to return a dictionary with Mulliken
            symbols as the keys.

        Returns
        -------
        Numpy array or dictionary

        Examples
        --------
        >>> rep = Reducible([9, -1, 3, 1], 'c2v', all_motion=True)
        >>> rep.vibe_modes()
        array([2, 0, 1, 0])
        >>> rep.vibe_modes(to_dict=True)
        {'A1': 2, 'A2': 0, 'B1': 1, 'B2': 0}
        """
        if self.all_motion is False:
            return self.decomp()
        else:
            rot_trans = rot_trans_modes[self.group]
            irreducibles = self.decomp()

        return np.array(irreducibles) - np.array(rot_trans)

    @_return_dict
    def ir_active(self, to_dict=False):
        """Return IR active vibrational modes.

        Return the number of each irreducible representation that is IR active
        in the given reducible representation. If all_motion=True for the
        Reducible, the rotational and translational modes are automatically
        subtracted out.

        Parameters
        ----------
        to_dict : bool
            True causes function to return a dictionary with Mulliken symbols
            as the keys.

        Returns
        -------
        Numpy array or dictionary

        Examples
        --------
        >>> rep = Reducible([9, -1, 3, 1], 'c2v', all_motion=True)
        >>> rep.ir_active()
        array([2, 0, 1, 0])
        >>> rep = Reducible([5, 2, 1, 3, 0, 3], 'd3h', all_motion=False)
        >>> rep.ir_active()
        array([0, 0, 1, 0, 1, 0])
        >>> rep.ir_active(to_dict=True)
        {"A'1": 0, "A'2": 0, "E'": 1, 'A"1': 0, 'A"2': 1, 'E"': 0}

        """
        return self.vibe_modes() * np.array(IR_active[self.group])

    @_return_dict
    def raman_active(self, to_dict=False):
        """Return Raman active vibrational modes.

        Return the number of each irreducible representation that is Raman
        active in the given reducible representation. If all_motion=True for
        the Reducible, the rotational and translational modes are automatically
        subtracted out.

        Parameters
        ----------
        to_dict : bool
            True causes the function to return a dictionary with Mulliken
            symbols as the keys.

        Returns
        -------
        Numpy array or dictionary

        Examples
        --------
        >>> rep = Reducible([9, -1, 3, 1], 'c2v', all_motion=True)
        >>> rep.raman_active()
        array([2, 0, 1, 0])
        >>> rep = Reducible([5, 2, 1, 3, 0, 3], 'd3h', all_motion=False)
        >>> rep.raman_active()
        array([2, 0, 1, 0, 0, 0])
        >>> rep.raman_active(to_dict=True)
        {"A'1": 2, "A'2": 0, "E'": 1, 'A"1': 0, 'A"2': 0, 'E"': 0}

        """
        return self.vibe_modes() * np.array(Raman_active[self.group])

    @classmethod
    def from_irred(cls, n_irred, group, all_motion=False):
        """Create reducible from number of irreducible representations.

        Alternative constructor that returns a Reducible representation object
        given the number of each irreducible representation that comprises
        the reducible representation in the point group.

        Parameters
        ----------
        n_irred: List, tuple, or array
            Number of each irreducible representation in the returned
            reducible representation.
        group: str
            Point group in Schoenflies notation (e.g., 'C2v').  This is
            case-insensitive.
        all_motion: bool
            Whether the resulting reducible representation represents all
            motions (rotation, vibration, and translation).

        Returns
        -------
        Reducible object

        Examples
        --------
        >>> rep = Reducible.from_irred([1, 0, 1, 0], 'c2v')
        >>> rep.gamma
        array([2, 0, 2, 0])
        >>> rep = Reducible.from_irred([3, 1, 1], 'C3v')
        >>> rep.gamma
        array([6, 3, 2])
        """
        # double each irreducible for groups with complex conjugates
        n_irred = np.asarray(n_irred)[np.cumsum(masks[group.lower()]) - 1]
        irred_sum = np.sum((tables[group.lower()].T * n_irred).T, axis=0)

        return cls(np.rint(_sympy_to_num(irred_sum).real).astype(int), group,
                   all_motion=all_motion)

    @classmethod
    def from_atoms(cls, n_atoms, group):
        """Create a representation based on number of stationary atoms.

        Alternative constructor that returns a Reducible representation object
        with all motions (rotation, translation, and vibration) given the
        number of atoms that do NOT move (i.e., translate) when carrying out
        each symmetry operation in the point group.

        Note: all_motion parameter in the returned Reducible object is set
        to True.

        Parameters
        ----------
        n_atoms: List, tuple, or array
            Number of atoms that remain stationary during each operation in
            the point group - see print_header() for the group operations.
        group: str
            Point group in Schoenflies notation (e.g., 'C2v').  This is
            case-insensitive.

        Returns
        -------
        Reducible object

        Examples
        --------
        >>> rep = Reducible.from_atoms([4, 2, 4, 2], 'c2v')
        >>> rep.gamma  # doctest: +NORMALIZE_WHITESPACE
        array([12, -2, 4, 2])

        """
        n_atoms = np.array(n_atoms)

        if np.any(np.mod(n_atoms, 1) != 0):
            raise ValueError('Number of stationary atoms (n_atoms) must be'
                             ' an integer value.')

        gamma = np.rint(
            n_atoms * _sympy_to_num(
                np.array(atom_contribution[group.lower()])).real)
        return cls(gamma.astype(int), group, all_motion=True)
