# -*- coding: utf-8 -*-

"""
Contains chemical group theory functions for calculating symmetry adapted
linear combinations (SALCs) or group orbitals using either the projection
operator method or using the symmetry functions in the character tables.
"""

from math import cos, sin, radians, isclose
from functools import wraps
import numpy as np
import sympy

from .tables import tables, symmetry_func_dict, column_coeffs, mulliken, masks


def _return_dict(func):
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

    @wraps(func)
    def wrapper(*args, **kwargs):
        if kwargs.get("to_dict"):
            group = kwargs.get("group", args[1] if len(args) > 1
                               else kwargs['group'])
            keys = mulliken[group.lower()]
            values = func(*args, **kwargs)
            return dict(zip(keys, values, strict=True))

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
    [2, -1, -1, 0, 0, 0]

    """
    expanded_irred = []
    for i in range(len(irred)):
        expanded_irred.extend([irred[i]] * column_coeffs[group.lower()][i])

    return expanded_irred


def _set_largest_coeff_pos(salcs):
    """Flip signs of coefficients in SALC if largest coefficient is negative.

    Flips all coefficient signs if the largest magnitude coefficient is
    negative. If the coefficients are imaginary, the salc is returned
    unchanged.

    Parameters
    ----------
    salcs : list or nested list or tuple of SymPy expressions
        SALC as a SymPy expression

    Returns
    -------
    list or nested list of salcs

    Example
    -------
    >>> from sympy import I, pi, exp
    >>> a, b, c, d = sympy.symbols('a b c d')
    >>> _set_largest_coeff_pos([-2*a + b + c])
    [2*a - b - c]
    >>> _set_largest_coeff_pos([2*a- b - c])
    [2*a - b - c]
    >>> _set_largest_coeff_pos([-a + b])
    [a - b]
    >>> _set_largest_coeff_pos([a + b*exp(2*I*pi/3) + c*exp(-2*I*pi/3)])
    [a + b*exp(2*I*pi/3) + c*exp(-2*I*pi/3)]
    >>> _set_largest_coeff_pos([-a + -b + -c, 0, [-a + c, -b + d]])
    [a + b + c, 0, [a - c, b - d]]
    """
    results = []
    for salc in salcs:
        if isinstance(salc, list):
            results.append(_set_largest_coeff_pos(salc))
        elif salc == 0:
            results.append(salc)
        else:
            symbols = sorted(list(salc.free_symbols), key=str)
            coeffs = [salc.coeff(s) for s in symbols]
            if not all(coeff.is_real for coeff in coeffs):
                results.append(salc)
            elif max(coeffs, key=abs) < 0:
                results.append(-salc)
            else:
                results.append(salc)
    return results


def _normalize_salcs_expr(salcs, symbols, normalize_by='largest'):
    """
    Normalize SALC composed of sympy expressions.

    Normalizes SALC by dividing each SALC through by the largest value in
    that SALC.

    Parameters
    ----------
    salcs : List or nested list of sympy expressions.
        Nested list of SALCs.
    symbols : SymPy symbols
        SymPy symbols representing outer ligands or atoms.
    normalize_by : 'largest' or 'smallest'
        Coefficient used to divide all other coefficients during
        normalization. The default is 'largest'.

    Returns
    -------
    list of normalized SALCs

    """
    normalized_values = []
    for salc in salcs:
        if isinstance(salc, list):
            normalized_values.append(_normalize_salcs_expr(
                salc, symbols, normalize_by=normalize_by))
        elif salc == 0:
            normalized_values.append(salc)
        else:
            coeffs = [salc.coeff(symbol) for symbol in symbols]
            if normalize_by == 'largest':
                norm_coeff = max(coeffs, key=lambda x: (
                    sympy.Abs(x), sympy.re(x), sympy.im(x)))
            elif normalize_by == 'smallest':
                nz_coeffs = (c for c in coeffs if c != 0)
                norm_coeff = min(nz_coeffs, key=lambda x: (sympy.Abs(x),
                                                           sympy.re(x), sympy.im(x)))
            normalized_values.append(salc / norm_coeff)

    return normalized_values


def _combine_conjugates(salcs, mask):
    """Group together SALCs from doubly-degenerate complex conjugates.

    Parameters
    ----------
    salcs : List
        Non-nested list of SALCs.
    mask : tuple
        Tuple of 1 and 0 values.

    Returns
    -------
    List or nested list
    """
    r = []
    for i in range(len(salcs)):
        try:
            if mask[i] and mask[i + 1]:
                r.append(salcs[i])
            elif mask[i] and not mask[i + 1]:
                r.append([salcs[i], salcs[i + 1]])
        except IndexError:
            r.append(salcs[i])
    return r


@_return_dict
def calc_salcs_projection(projection, group, *, to_dict=False,
                          normalize_by='largest'):
    """
    Return SALCs using projection operator method.

    Given the projections of orbitals as a result of a point group symmetry
    operations, returns the SALCs. This is a three-step process.
    1. Provide all ligands or outer atoms with SymPy variable names.
    2. Track an orbital to see how it transforms after each symmetry operation
    3. Provide a list of the results from step 2.

    Note: The projection operator method only returns one SALC for E and T
    point groups.

    Parameters
    ----------
    projection : List, tuple, or array of SymPy symbols
        Results of projection operations for symmetry operations of group.
    group : str
        Point group Schoenflies notation (e.g., 'C2v').  This is
        case-insensitive.
    to_dict : bool
        True causes the function to return a dictionary with Mulliken
        symbols as the keys.
    normalize_by : 'largest' or 'smallest'
        Coefficient used to divide all other coefficients during
        normalization. The default is 'largest'.

    Returns
    -------
    List or nested list of sympy expressions of the SALCs for each irreducible
    representation. Returns 0 for irreducibles with no SALC. If to_dict=True,
    returns a dictionary.

    Example
    -------
    >>> import sympy
    >>> a, b, c = sympy.symbols('a b c')
    >>> calc_salcs_projection([a, b, c, a, b, c], 'c3v')
    [a + b + c, 0, a - b/2 - c/2]
    >>> calc_salcs_projection([a, b, c, a, b, c], 'c3v', to_dict=True)
    {'A1': a + b + c, 'A2': 0, 'E': a - b/2 - c/2}

    """
    if not isinstance(group, str):
        raise ValueError('group needs to be a string.')
    if group.lower() not in tables.keys():
        raise ValueError('Invalid point group.')
    if normalize_by not in ('smallest', 'largest'):
        raise ValueError('normalize_by must be "largest" or "smallest".')
    if not all([isinstance(proj, sympy.Expr) and not
                isinstance(proj, sympy.Number) for proj in projection]):
        raise ValueError('Projection must be all SymPy symbols.')
    n_ops = sum(column_coeffs[group.lower()])
    if n_ops != len(projection):
        raise ValueError(f'The projection length must be {n_ops} for {group}.')

    salcs = []

    for irred in tables[group.lower()]:
        product = np.array(
            _expand_irreducible(irred, group.lower()) * np.array(projection)
        )
        salcs.append(np.sum(product))
    salcs = _combine_conjugates(
        salcs, masks[group.lower()])  # next complx conj
    unique_sym = set().union(*(p.free_symbols for p in projection))
    symbols = sorted(unique_sym, key=str)

    return _set_largest_coeff_pos(
        _normalize_salcs_expr(salcs, symbols, normalize_by=normalize_by))


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
    [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]

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


def _eval_sym_func(coords, exprs):
    """
    Evaluate symmetry functions for an irreducible representation.

    Evaluates symmetry functions using the supplied series of xyz coordinates.
    If all values evaluate as zero, the irreducible has no SALC, and 0 is
    returned.

    Parameters
    ----------
    coords : List, tuple, or array containing values in threes
        xyz coordinates of ligand unit vectors.
    exprs : tuple of strings
        The symmetry function supplied as a string or tuple of strings
        (e.g., 'x**2-y**2' or ('z**2', 'x**2+y**2')).

    Returns
    -------
    List or 0.

    """
    salcs = []

    for expr in exprs:
        ligand_contribs = []
        for unit_vector in coords:
            x, y, z = unit_vector[0], unit_vector[1], unit_vector[2]
            ligand_contrib = eval(expr, {'x': x, 'y': y, 'z': z})
            ligand_contribs.append(round(ligand_contrib, 3))
        if np.any(ligand_contribs):
            salcs.append(ligand_contribs)

    if not salcs:
        return 0
    if len(salcs) == 1:
        return salcs[0]
    return salcs


def _normalize_salcs(salcs, normalize_by='largest'):
    """
    Normalize SALC.

    Normalizes SALC by dividing each SALC through by the largest value in
    that SALC.

    Parameters
    ----------
    salcs : List or nested list
        Nested list of SALCs.
    normalize_by : 'largest' or 'smallest'
        Coefficient used to divide all other coefficients during
        normalization. The default is 'largest'.

    Returns
    -------
    list of normalized SALCs

    """
    normalized_values = []
    for value in salcs:
        if isinstance(value, list):
            normalized_values.append(_normalize_salcs(value, normalize_by))
        elif value == 0:
            normalized_values.append(value)
        else:
            if normalize_by == 'largest':
                coeff = round(value / max(salcs, key=abs), 2)
            elif normalize_by == 'smallest':
                coeff = round(
                    value / min((s for s in salcs if s != 0), key=abs), 2)
            if coeff % 1 < 0.01:
                coeff = sympy.Integer(coeff)
            normalized_values.append(coeff)

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
    [a + 2*b + c, b - c, 0]

    """
    symbolic_wt = []
    for weight in weights:
        if weight == 0:
            symbolic_wt.append(0)
        elif isinstance(weight, (list, tuple)) and isinstance(weight[0], (list, tuple)):
            symbolic_wt.append(_weights_to_symbols(weight, symbols))
        else:
            sym = np.array(weight).dot(np.array(symbols))
            if isinstance(sym, np.ndarray):
                symbolic_wt.append(sym.tolist())
            else:
                symbolic_wt.append(sym)
    return symbolic_wt


@_return_dict
def calc_salcs_func(ligands, group, symbols, *, mode="vector", to_dict=False,
                    normalize_by='largest'):
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
    ligands : list or nested list
        Nested list of ligand positions as xyz coordinates (mode='vector')
        or angles (mode='angle').
    group : str
        Point group Schoenflies notation (e.g., 'C2v').  This is
        case-insensitive.
    symbols : SymPy symbols
        SymPy symbols representing outer ligands or atoms.
    mode : 'vector' or 'angle'
        Whether the position of ligands or outer atoms are provided in xyz
        coordinates ('vector') or [azimuthal, polar] angles ('angle'). The
        azimuthal angle is the angle from the positive x-axis on the xy-plane
        and the polar angle is the angle from the positive z-axis. Default is
        'vector'.
    to_dict : bool
        True causes the function to return a dictionary with Mulliken
        symbols as the keys.
    normalize_by : 'largest' or 'smallest'
        Coefficient used to divide all other coefficients during
        normalization. The default is 'largest'.

    Returns
    -------
    List of sympy symbols indicating the weight and sign of each atomic
    orbital contribution to the SALC. There may be redundant SALCs returned
    due to multiple symmetry functions with an irreducible representation
    returning the same SALC. Not all SALCs are guaranteed to be returned in
    this method.

    Examples
    --------
    >>> import sympy
    >>> a, b, c, d = sympy.symbols('a b c d')
    >>> coords = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    >>> calc_salcs_func(coords, 'd4h', [a, b, c, d], mode='vector')
    [a + b + c + d, 0, a - b + c - d, 0, 0, 0, 0, 0, 0, [a - c, b - d]]
    >>> coords = [[0, -90], [120, -90], [240, -90]]
    >>> calc_salcs_func(coords, 'd3h', [a, b, c], mode='angle')
    [a + b + c, 0, [a - 0.5*b - 0.5*c, b - c, a - 0.5*b - 0.5*c, b - c], \
0, 0, 0]

    Notes
    -----
    C1 has only E as a symmetry operation, so each orbital is its own SALC.

    References
    ----------
    [1] Kim, S. K. Group Theoretical Methods and Applications to Molecules
    and Crystals; Cambridge University Press: Cambridge, 1999, 155-157.
    """
    if not isinstance(group, str):
        raise ValueError('group needs to be a string.')
    if group.lower() not in tables.keys():
        raise ValueError('Invalid point group.')
    if group.lower() == 'c1':
        raise ValueError('Each orbital is its own SALC for C1.')
    if normalize_by not in ('smallest', 'largest'):
        raise ValueError('normalize_by must be "largest" or "smallest".')
    if len(ligands) != len(symbols):
        raise ValueError('The number of symbols must equal the ligands.')

    if mode == "angle":
        ligand_vectors = _angles_to_vectors(ligands)
    elif mode == "vector":
        ligand_vectors = ligands
    else:
        raise ValueError("Invalid mode input: must be 'angle' or 'vector'")

    salcs = []
    for sym_func in symmetry_func_dict[group.lower()]:
        if sym_func == 0:
            salcs.append(0)
        else:
            salcs.append(_eval_sym_func(ligand_vectors, sym_func))

    normalized = _normalize_salcs(salcs, normalize_by=normalize_by)

    return _set_largest_coeff_pos(_weights_to_symbols(normalized, symbols))
