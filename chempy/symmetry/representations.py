# -*- coding: utf-8 -*-

"""
Contains chemical group theory functions for calculating numbers and type of
irreducible represetantions in a reducible representation and number of IR
and Raman active modes. Below uses the inverse matrix method of solving for
irreducible representations described in J. Chem. Educ. 2009, 86, 251-253
https://doi.org/10.1021/ed086p251.
"""

import numpy as np

# data tables in character_tables.py file
from .tables import (
    tables,
    rot_trans_modes,
    IR_active,
    Raman_active,
    masks,
    atom_contribution,
    mulliken,
    headers
)


def get_header(group):
    """
    Print the header for a character table given the point group.

    Parameters
    ----------
    group : string
        Point group Schoelflies notation (e.g., 'C2v').  This is
        case-insentive.

    Returns
    -------
    Prints character table header indicating order of symmetry operations
    """
    print(*headers[group.lower()])


def get_point_groups():
    """
    Print supported point group Schoenflies notations.

    Parameters
    ----------
    None

    Returns
    -------
    Prints Schoelflies notations for character tables available

    """
    pg = ('C1', 'Cs', 'Ci', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'D2',
          'D3', 'D4', 'D5', 'D6', 'C2v', 'C3v', 'C4v', 'C5v', 'C6v', 'C2h',
          'C3h', 'C4h', 'C5h', 'C6h', 'D2h', 'D3h', 'D4h', 'D5h', 'D6h', 'D8h',
          'D2d', 'D3d', 'D4d', 'D5d', 'D6d', 'S4', 'S6', 'S8', 'T', 'Th', 'Td',
          'O', 'Oh', 'I', 'Ih')

    print(*pg)


def get_mulliken(group):
    """
    Print mulliken symbols of irreducible reprensentation in order.

    Parameters
    ----------
    group : string
        Point group Schoelflies notation (e.g., 'C2v').  This is
        case-insentive.

    Returns
    -------
    Prints mulliken notations for irreducible representations

    """
    print(*mulliken[group.lower()])


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
        gamma : array_like
            Reducible representation.
        group : str
            Point group Schoelflies notation (e.g., 'C2v').  This is
            case-insentive.
        all_motion : bool, optional
            False if reducible representation is only vibrational modes and
            True if the reducible describes all motions (rotational,
            translational, and vibrational). The default is False.

        Returns
        -------
        None.
        """
        if np.all(np.mod(gamma, 1) != 0):
            print('Invalid representation - must be whole numbers.')
        elif group.lower() not in tables.keys():
            print('Invalid point group.')
        elif np.array(gamma).size != tables[group.lower()].shape[0]:
            print(f'Invalid representation size for {group} point group.')
        else:
            self.group = group.lower()
            self.gamma = gamma
            self.all_motion = all_motion

    def return_dict_method(func):
        """
        Return results as a dictionary.

        Return a list or array as a dictionary with with mulliken symbols as
        the keys.

        Parameters
        ----------
        arr : array_like
            List or array containing results corresponding to irreducible
            representations.
        group : string
            Point group Schoelflies notation (e.g., 'C2v').  This is
            case-insentive.

        Returns
        -------
        Dictionary.

        """
        def wrapper(self, *args, **kwargs):
            if kwargs.get('to_dict'):
                keys = mulliken[self.group.lower()]
                values = func(self, *args, **kwargs)
                return dict(zip(keys, values.tolist()))
            else:
                return func(self, *args, **kwargs)
        return wrapper

    @return_dict_method
    def decomp(self, to_dict=False):
        """
        Decompose reducible representation into number of irreducbiles.

        Decompose a reducible representation for a specific point group and
        returns the number of each irreducible representation in the reducible.

        The order of irreducibles can be determined by get_mulliken(group).

        Parameters
        ----------
        gamma : array_like
               Reducible representation with symmetry operations in the order
               provided by get_header(group).
        group : str
            Point group of representation in Schoelflies notation
            (e.g., 'C2v'). This is case-insentive.
        to_dict : bool
            True causes function to return a dictionary with mulliken symbols
            as the keys.

        Returns
        -------
        np.array with  number of each irreducible representation in the
        provided reducible representation. Use get_mulliken() for listing of
        irreducible representations and order.

        Examples
        --------
        >>> rep = Reducible([9, -1, 3, 1], 'c2v', all_motion=True)
        >>> rep.decomp()
        array([3, 1, 3, 2])

        >>> rep = Reducible([15, 0, 0, 7, -2, -2], 'C3h', all_motion=False)
        >>> rep.decomp()
        array([3, 4, 1, 1])
        >>> rep.decomp(to_dict=True)
        >>> {"A'": 3, 'A"': 4, "E'": 2, 'E"': 1}
        """
        table = tables[self.group]
        gamma = np.array(self.gamma)

        if self.group == 'c1':
            n_i = self.gamma
        else:
            mask = np.array(masks[self.group], dtype=bool)
            n_i = gamma.dot(np.linalg.inv(table)).real[mask]

        return np.rint(n_i).astype(int)

    @return_dict_method
    def vibe_modes(self, to_dict=False):
        """Return vibrational modes.

        Return the number of vibrational modes after rotation and translation
        are subtracted out.

        Parameters
        ----------
        to_dict : bool
            True causes function to return a dictionary with mulliken symbols
            as the keys.

        Returns
        -------
        Numpy array or dictionary

        Examples
        --------
        >>> rep = Reducible([9, -1, 3, 1], 'c2v', all_motion=True)
        >>> rep.vide_modes()
        array([2, 0, 1, 0])
        >>> rep.vibe_modes(to_dict)
        >>> {'A1': 2, 'A2': 0, 'B1': 1, 'B2': 0}
        """
        if self.all_motion is False:
            return self.decomp()
        else:
            rot_trans = rot_trans_modes[self.group]
            irreducibles = self.decomp()

        return np.array(irreducibles) - np.array(rot_trans)

    @return_dict_method
    def ir_active(self, to_dict=False):
        """Return IR active vibrational modes.

        Return the number of each irreducible representation that are IR active
        in the given reducible representation. If vibe_only=False for the
        reducible representation, the rotational and translational modes are
        automatically subtracted out.

        Parameters
        ----------
        to_dict : bool
            True causes function to return a dictionary with mulliken symbols
            as the keys.

        Returns
        -------
        Numpy array or dictionary

        Examples
        --------
        >>> rep = Reducible([9, -1, 3, 1], 'c2v', all_motion=True)
        >>> rep.ir_active([3, 1, 3, 2, 'C2v')
        array([2, 0, 1, 0])
        >>> rep = Reducible([5, 2, 1, 3, 0, 3], 'd3h', all_motion=False)
        >>> rep.ir_active()
        >>> array([0, 0, 1, 0, 1, 0])
        >>> rep.ir_active(to_dict=True)
        >>>  {"A1'": 0, "A2'": 0, "E'": 1, 'A1"': 0, 'A2"': 1, 'E"': 0}

        """
        return self.vibe_modes() * np.array(IR_active[self.group])

    @return_dict_method
    def raman_active(self, to_dict=False):
        """Return Raman active vibrational modes.

        Return the number of each irreducible representation that are Raman
        active in the given reducible representation. If vibe_only=False for
        the reducible representation, the rotational and translational modes
        are automatically subtracted out.

        Parameters
        ----------
        to_dict : bool
            True causes function to return a dictionary with mulliken symbols
            as the keys.

        Returns
        -------
        Numpy array or dictionary

        Examples
        --------
        >>> rep = Reducible([9, -1, 3, 1], 'c2v', all_motion=True)
        >>> rep.raman_active([3, 1, 3, 2, 'C2v')
        array([2, 0, 1, 0])
        >>> rep = Reducible([5, 2, 1, 3, 0, 3], 'd3h', all_motion=False)
        >>> rep.raman_active()
        >>> array([2, 0, 1, 0, 0, 0])
        >>> rep.raman_active(to_dict=True)
        >>>  {"A1'": 2, "A2'": 0, "E'": 1, 'A1"': 0, 'A2"': 0, 'E"': 0}

        """
        return self.vibe_modes() * np.array(Raman_active[self.group])

    @classmethod
    def from_irred(cls, n_irred, group, all_motion=False):
        """Create reducible from number of irreducibiel representations.

        Alternative constructor that returns a reducible representation
        given the number of each irreducible representations that comprise
        the reducible representation and the point group.

        Parameters
        ----------
        n_irred: array_like
            Number of each irreducible representations in the returned
            reducible representation.
        group: str
            Point group in Schoelflies notation (e.g., 'C2v').  This is
            case-insentive.
        all_motion: bool
            Whether the resulting reducible representation represents all
            motions (rotration, vibration, and translation).

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
        irred_sum = np.sum((tables[group.lower()].T * n_irred).T, axis=0)

        return cls(np.rint(irred_sum).astype(int), group, all_motion=False)

    @classmethod
    def from_atoms(cls, n_atoms, group):
        """Create a representation based on number of stationary atoms.

        Alternative constructor that returns a reducible representation
        with all motions (rotation, translation, and vibration) given the
        number of atoms that do NOT move (i.e., translate) when carrying out
        each symmetry operation in the point group. Note: vibe_only parameter
        in the returned Reducible object is set to False.

        Parameters
        ----------
        n_atoms: array_like
            Number of atoms that remain stationary during each operation in
            the point group - see get_header() for list of operations.
        group: str
            Point group in Schoelflies notation (e.g., 'C2v').  This is
            case-insentive.


        Returns
        -------
        Reducible object

        Examples
        --------
        >>> rep = Reducible.from_atoms([4, 2, 4, 2], 'c2v')
        >>> rep.gamma
        >>> array([12, -2,  4,  2])
        """
        n_atoms = np.array(n_atoms)

        if np.all(np.mod(n_atoms, 1) == 0):
            gamma = np.rint(n_atoms * np.array(atom_contribution[group]))
            return cls(gamma.astype(int), group, all_motion=True)
        else:
            print("""Number of stationary atoms (n_atoms) must be an integer
                  value.""")
