# -*- coding: utf-8 -*-
"""
Utility functions related to stoichiometry.
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np


def get_coeff_mtx(substances, stoichs):
    """
    Create a net stoichiometry matrix from reactions
    described by pairs of dictionaries.

    Parameters
    ----------
    substances: sequence of keys in stoichs dict pairs
    stoichs: sequence of pairs of dicts
        pairs of reactant and product dicts mapping substance keys
        to stoichiometric coefficients (integers)

    Returns
    -------
    2 dimensional array of shape (len(substances), len(stoichs))

    """
    A = np.zeros((len(substances), len(stoichs)), dtype=int)
    for ri, sb in enumerate(substances):
        for ci, (reac, prod) in enumerate(stoichs):
            A[ri, ci] = prod.get(sb, 0) - reac.get(sb, 0)
    return A


def decompose_yields(yields, stoichs, atol=1e-10):
    """ Decomposes yields into mass-action reactions

    This function offers a way to express a reaction with non-integer
    stoichiometric coefficients as a linear combination of production reactions
    with integer coefficients.

    Ak = y

    A is (n_species x n_reactions) matrix, k is "rate coefficient", y is yields


    Parameters
    ----------
    yields: OrderedDict
        specie names as keys and yields as values
    stoichs: iterable dictionary pairs
        dict keys must match those of `yields`each pair
        of dictionaries gives stoichiometry
        (1st is reactant, 2nd is products)
    atol: float
        absolute tolerance for residuals

    Raises
    ------
    ValueError
        When atol is exceeded
    numpy.LinAlgError
        When numpy.linalg.lstsq fails to converge

    Returns
    -------
    1-dimensional array of effective rate coefficients.

    """
    # Sanity check:
    for ys in yields.keys():
        present = False
        for reac, prod in stoichs:
            if ys in reac or ys in prod:
                present = True
        assert present

    sbstncs = yields.keys()
    y = np.array(list(yields.values()))
    A = get_coeff_mtx(sbstncs, stoichs)
    k, residuals, rank, s = np.linalg.lstsq(A, y)
    if len(residuals) > 0:
        if np.any(residuals > atol):
            raise ValueError("atol not satisfied")
    return k
