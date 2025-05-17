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
    headers,
    mulliken,
    rot_trans_modes,
    IR_active,
    Raman_active,
    masks,
    atom_contribution,
)


def get_header(group):
    """
    Print the header for a character table given the point group.

    Parameters
    ----------
    group : string
        point group Schoelflies notation (e.g., 'C2v')

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
        point group Schoelflies notation (e.g., 'C2v')

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

    def __init__(self, gamma, group, vibe_only=True):
        """

        Parameters
        ----------
        gamma : array_like
            Reducible representation.
        group : str
            Point group Schoelflies notation (e.g., 'C2v')
        vibe_only : bool, optional
            True if reducible representation is only vibrational modes and
            False if the reducible describes all motions. The default is True.

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
            self.vibe_only = vibe_only

    def decomp_reduc(self):
        """
        Decompose reducible representation.

        Decompose a reducible representation for a specific point group and
        returns the number of each irreducible representation in the reducible.

        Parameters
        ----------
        gamma: array_like
               Reducible representation
        group: str
            Point group of representation in Schoelflies notation (e.g., 'C2v')

        Returns
        -------
        Numpy array with  number of each irreducible representation, in
        order, in the provided reducible representation. Use get_mulliken()
        for listing of irreducible representations and order.

        Examples
        --------
        >>> rep = ReduceRep([9, -1, 3, 1], 'c2v', vibe_only=False)
        >>> rep.decomp_reduc([9, -1, 3, 1], 'C2v')
        array([3, 1, 3, 2])

        >>> rep = ReduceRep([15, 0, 0, 7, -2, -2], 'c3h, vibe_only=True)
        >>> rep.decomp_reduc([15, 0, 0, 7, -2, -2], 'C3h')
        array([3, 4, 1, 1])

        """
        table = tables[self.group]
        gamma = np.array(self.gamma)

        if self.group == 'c1':
            n_i = self.gamma
        else:
            mask = np.array(masks[self.group], dtype=bool)
            n_i = gamma.dot(np.linalg.inv(table)).real[mask]

        return np.rint(n_i).astype(int)

    def vibe_modes(self):
        """Return vibrational modes.

        Return the number of vibrational modes after rotation and translation
        are subtracted out.

        Parameters
        ----------
        None

        Returns
        -------
        Numpy array

        Examples
        --------
        >>> rep = ReduceRep([9, -1, 3, 1], 'c2v', vibe_only=False)
        >>> rep.vide_modes()
        array([2, 0, 1, 0])
        """
        if self.vibe_only is True:
            return self.decomp_reduc()
            print('Already only contains vibrational modes.')
        else:
            rot_trans = rot_trans_modes[self.group.lower()]
            irreducibles = self.decomp_reduc()

        return np.array(irreducibles) - np.array(rot_trans)

    def ir_active(self):
        """Retrun IR active modes.

        Return the number of each irreducible representation that are IR active
        in the given reducible representation. If vibe_only=False for the
        reducible representation, the rotational and translational modes are
        automatically subtracted out.

        Returns
        -------
        Numpy array

        Examples
        --------
        >>> rep = ReduceRep([9, -1, 3, 1], 'c2v', vibe_only=False)
        >>> rep.ir_active([3, 1, 3, 2, 'C2v')
        array([2, 0, 1, 0])

        """
        irreducibles = self.decomp_reduc()

        if self.vibe_only is True:
            return irreducibles * np.array(IR_active[self.group])

        if self.vibe_only is False:
            rot_trans = np.array(rot_trans_modes[self.group.lower()])

            vibrations = np.array(irreducibles) - rot_trans

            return vibrations * np.array(IR_active[self.group.lower()])

    def raman_active(self):
        """Return Raman active modes.

        Return the number of each irreducible representation that are Raman
        active in the given reducible representation. If vibe_only=False for
        the reducible representation, the rotational and translational modes
        are automatically subtracted out.

        Returns
        -------
        Numpy array

        Examples
        --------
        >>> rep = ReduceRep([9, -1, 3, 1], 'c2v', vibe_only=False)
        >>> rep.raman_active([3, 1, 3, 2, 'C2v')
        array([2, 0, 1, 0])

        """
        rot_trans = np.array(rot_trans_modes[self.group.lower()])
        irreducibles = self.decomp_reduc()

        vibrations = np.array(irreducibles) - rot_trans

        if np.any((vibrations % 1) != 0):
            print(""""Invalid reducible representation - all values must
                  be whole numbers.""")
        else:
            return vibrations * np.array(Raman_active[self.group.lower()])

    @classmethod
    def from_irred(cls, n_irred, group, vibe_only=True):
        """Create reducible for number of irreducibiel representations.

        Alternative constructor that returns a reducible representation
        given the number of each irreducible representations that comprise
        the reducible representation and the point group.

        Parameters
        ----------
        n_irred: array_like
            Number of each irreducible representations in the returned
            reducible representation
        group: str
            Point group in Schoelflies notation (e.g., 'C2v')
        vibe_only: bool
            Whether the resulting reducible representation only represents
            vibrational modes

        Returns
        -------
        Reducible object

        Examples
        --------
        >>> rep = Reducible.from_irred([1, 0, 1, 0], 'c2v')
        >>> rep.gamma
        array([2, 0, 2, 0])

        >>> rep = Reducible.from_irred([3, 1, 3, 2], 'C2v')
        >>> rep.gamma
        array([9, -1, 3, 1])

        >>> rep = Reducible.from_irred([3, 1, 1], 'C3v')
        >>> rep.gamma
        array([6, 3, 2])

        """
        irred_sum = np.sum((tables[group.lower()].T * n_irred).T, axis=0)

        return cls(np.rint(irred_sum).astype(int), group, vibe_only=True)

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
            Point group in Schoelflies notation (e.g., 'C2v')


        Returns
        -------
        Reducible object
        """
        n_atoms = np.array(n_atoms)

        if np.all(np.mod(n_atoms, 1) == 0):
            gamma = np.rint(n_atoms * np.array(atom_contribution[group]))
            return cls(gamma.astype(int), group, vibe_only=False)

        else:
            print("""Number of stationary atoms (n_atoms) must be an integer
                  value.""")
