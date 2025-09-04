# -*- coding: utf-8 -*-

from .tables import (
    tables,
    symmetry_func_dict,
    column_coeffs,
    mulliken
)

"""
Contains chemical group theory functions for calculating symmetry adapted
linear combinations (SALCs) or group orbitals using either the projection
operator method or using the symmetry functions in the character tables.
"""

from math import cos, sin, radians, isclose
import numpy as np
import sympy
sympy.init_printing(pretty_print=False)


def return_dict(func):
    """
    Return results as a dictionary.

    Return a list or array as a dictionary with Mulliken symbols as
    the keys.

    Parameters
    ----------
    func : function returning SALCs
        Function that returns SALCs as a list by default.

    Returns
    -------
    Dictionary.

    """
    def wrapper(*args, **kwargs):
        if kwargs.get('to_dict'):
            print(args)
            keys = mulliken[args[1].lower()]
            values = func(*args, **kwargs)
            return dict(zip(keys, values))
        else:
            return func(*args, **kwargs)
    return wrapper


# PROJECTION OPERATOR METHOD

def _expand_irreducible(irred, group):
    """
    Return expanded irreducible.

    Expands irreducible from the common condensed form to the full form. This
    entails repeating characters for every equivalent operation. For example,
    C3v returns the rotation twice and the reflection thrice.

    Parameters
    ----------
    irred : tuple
        Condensed irreducible representation from character table.
    group : str
        Point group Schoenflies notation (e.g., 'C2v').  This is
        case-insensitive.

    Returns
    -------
    Expanded irreducible representation as a list.

    Example
    ------
    >>> _expand_irreducible((2, -1, 0), 'c3v')
    >>> [2, -1, -1, 0, 0, 0]

    """
    expanded_irred = []
    for i in range(len(irred)):
        expanded_irred.extend([irred[i]] * column_coeffs[group.lower()][i])

    return expanded_irred


@return_dict
def calc_salcs_projection(projection, group, to_dict=False):
    """
    Return SALCs using projection operator method.

    Given the projections of orbitals as a result of a point group symmetry
    operations, returns the SALCs. This is a two-step process.
    1. Provide all ligands or outer atoms with SymPy variable names.
    2. Track an orbital to see how it transforms after each symmetry operation
    3. Provide a list of the results from step 2.

    Note: The projection operator method only turns one SALC for E and T
    point groups.

    Parameters
    ----------
    projection : List, tuple, or array of SymPy symbols
        Results of projection operations for symmetry operations of group.
    group : str
        Point group Schoenflies notation (e.g., 'C2v').  This is
        case-insensitive.

    Returns
    -------
    List or nested list of strings of the SALCs for each irreducible
    representation. Returns 0 for irreducibles with no SALC. If to_dict=True,
    returns a dictionary.

    Example
    -------
    >>> import sympy
    >>> a, b, c = sympy.symbols('a b c')
    >>> calc_salcs_projection([a, b, c, a, b, c], 'c3v')
    >>> [2*a + 2*b + 2*c, , 0, 2*a - b - c]
    >>> calc_salcs_projection([a, b, c, a, b, c], 'c3v', to_dict=True)
    >>> {'A1': 2*a + 2*b + 2*c, 'A2': 0, 'E': 2*a - b - c}

    """
    salcs = []

    for irred in tables[group.lower()]:
        product = np.array(_expand_irreducible(irred, group.lower()) *
                           np.array(projection))
        salcs.append(np.sum(product))

    return salcs


# USING SYMMETRY FUNCTIONS

def _angles_to_vectors(ligand_angles):
    """
    Calculate xyz vectors from angles around central atom.

    Given angles for outer ligands/atoms with respect to x-axis and z-axis,
    this function returns a list of [x, y, z] unit vectors.

    Parameters
    ----------
    ligand_angles : nested list or tuple of numbers in (theta, phi) order
        Angle of each outer atom/ligand in degrees with respect
        azimuth (theta) and polar (phi) coordinates. Theta is the angle from
        the positive x-axis on the xy-plane and phi is the angle from the
        positive z-axis.

    Returns
    -------
    Nested list containing xyz vectors.

    Example
    -------
    >>> _angles_to_vectors([[0, 90], [90, 90], [180, 90], [-90, 90]])
    >>> [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
         [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]])

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Spherical_coordinate_system
    """
    all_vectors = []
    for orbital in ligand_angles:
        azimuth, polar = radians(orbital[0]), radians(orbital[1])
        x = sin(polar) * cos(azimuth)
        y = sin(polar) * sin(azimuth)
        z = cos(polar)

        # convert to int if close
        vector = []
        for val in (x, y, z):
            if isclose(val, int(val), abs_tol=1e-3):
                vector.append(int(val))
            else:
                vector.append(round(val, 3))

        all_vectors.append(vector)

    return all_vectors


def _eval_sym_func(coords, funcs):
    """
    Evaluate symmetry functions for an irreducible representation.

    Evaluates symmetry functions using the supplied series of xyz coordinates.
    If all values evaluate as zero, the irreducible has no SALC, and 0 is
    returned.

    Parameters
    ----------
    coords : List, tuple, or array containing values in threes
        xyz coordinates of ligand unit vectors.
    funcs : str
        The symmetry function supplied as a string or tuple of strings
        (e.g., 'x**2-y**2' or ('z**2', 'x**2+y**2')).

    Returns
    -------
    List or 0.

    """
    salcs = []

    for func in funcs:
        ligand_contribs = []
        for unit_vector in coords:
            x, y, z = unit_vector[0], unit_vector[1], unit_vector[2]
            ligand_contrib = eval(func, {'x': x, 'y': y, 'z': z})
            ligand_contribs.append(round(ligand_contrib, 2))

        if np.any(ligand_contribs):
            salcs.append(ligand_contribs)

    if not salcs:
        return 0
    elif len(salcs) == 1:
        return salcs[0]
    else:
        return salcs


def _normalize_salcs(salcs):
    """
    Normalize SALC.

    Normalizes SALC by dividing each SALC through by the largest value in
    that SALC.

    Parameters
    ----------
    salcs : List or nested list
        Nested list of SALCs.

    Returns
    -------
    np.array

    """
    normalized_values = []
    for value in salcs:
        if isinstance(value, list):
            normalized_values.append(_normalize_salcs(value))
        elif isinstance(value, int):
            normalized_values.append(value)
        else:
            normalized_values.append(round(value / max(salcs, key=abs), 2))

    return normalized_values


def _weights_to_symbols(weights, symbols):
    """
    Convert ligand weights to symbolic representations.

    Convert ligand or outer atom weights (e.g., [2, -0.5, -0.5]) to a symbolic
    representations (e.g., [2*a, -0.5*b, -0.5*c]).

    Parameters
    ----------
    weights : list
        List or nested list containing weights from each ligand or outer atom.
    symbols : list of Sympy symbols
        List of SymPy symbols the user provided to represent the ligands or
        outer atoms.

    Returns
    -------
    List containing symbolic weights of each ligand or outer atom.

    Examples
    --------
    >>> import sympy
    >>> a, b, c = sympy.symbols('a b c')
    >>> _weights_to_symbols([[1, 2, 1], [0, 1, -1], [0, 0, 0]], [a, b, c])
    >>> [a + 2*b + c, b - c, 0]

    """
    symbolic_wt = []
    for weight in weights:
        if weight == 0:
            symbolic_wt.append(0)
        else:
            try:
                sym = np.array(weight).dot(np.array(symbols))
                if isinstance(sym, np.ndarray):
                    symbolic_wt.append(sym.tolist())
                else:
                    symbolic_wt.append(sym)
            except ValueError:
                symbolic_wt.append(_weights_to_symbols(weight, symbols))

    return symbolic_wt


@return_dict
def calc_salcs_func(ligands, group, symbols, mode='vector', to_dict=False):
    """
    Return SALCs using symmetry functions in character table.

    Returns SALCs (symmetry adapted linear combinations) using the character
    functions for a given point group. This requires the user to give the
    position of each ligand atomic orbital with respect to the x-axis, y-axis,
    and z-axis (mode='vector') or the spherical coordinate angles in degrees
    in pairs (azimuthal, polar) where the azimuthal is the angle from the
    positive x-axis on the xy-plane and polar is the angle from the positive
    z-axis (mode='angle').

    Parameters
    ----------
    ligand : list or nested list
        Nested list of ligand positions as xyz coordinates (mode='vector')
        or angles (mode='angle').
    symbols : SymPy symbols
        SymPy symbols representing outer ligands or atoms.
    group : str
        Point group Schoenflies notation (e.g., 'C2v').  This is
        case-insensitive.
    mode : 'vector' or 'angle'
        Whether the position of ligands or outer atoms are provided in xyz
        coordinates ('vector') or [theta, phi] angles ('angle').

    Returns
    -------
    List of sympy symbols indicating the weight and sign of each atomic
    orbital contribution to the SALC. There may be redundant SALCs returned
    due to multiple symmetry functions with an irreducible representation
    returning the same SALC.

    Examples
    --------
    >>> import sympy
    >>> a, b, c, d = sympy.symbols('a b c d')
    # for a square planar complex
    >>> calc_salcs_func([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
                            'd4h', [a, b, c, d], mode='vector')
    >>> [a + b + c + d, 0, a - b + c - d, 0, 0, 0, 0, 0, 0, [a - c, b - d]]
    # for a trigonal planar complex
    >>> calc_salcs_func([[0, -90], [120, -90], [240, -90]], 'd3h', [a, b, c],
                        mode='angle')
    >>> [1.0*a + 1.0*b + 1.0*c, 0,
         [1.0*a - 0.5*b - 0.5*c,
          1.0*b - 1.0*c, 1.0*a - 0.5*b - 0.5*c,
          1.0*b - 1.0*c], 0, 0, 0]

    References
    ----------
    [1] Kim, S. K. Group Theoretical Methods and Applications to Molecules
    and Crystals; Cambridge University Press: Cambridge, 1999, 155-157.
    """
    if mode == 'angle':
        ligand_vectors = _angles_to_vectors(ligands)
    elif mode == 'vector':
        ligand_vectors = ligands
    else:
        raise Exception("Invalid mode input: must be 'angle' or 'vector'")

    salcs = []
    for sym_func in symmetry_func_dict[group]:
        if sym_func == 0:
            salcs.append(0)
        elif sym_func == 1:
            salcs.append(1)
        else:
            salcs.append(_eval_sym_func(ligand_vectors, sym_func))

    return _weights_to_symbols(_normalize_salcs(salcs), symbols)
