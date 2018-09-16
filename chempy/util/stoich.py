# -*- coding: utf-8 -*-
"""
Utility functions related to stoichiometry.
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from chempy.units import unit_of, to_unitless


def get_coeff_mtx(substances, stoichs):
    """
    Create a net stoichiometry matrix from reactions
    described by pairs of dictionaries.

    Parameters
    ----------
    substances : sequence of keys in stoichs dict pairs
    stoichs : sequence of pairs of dicts
        Pairs of reactant and product dicts mapping substance keys
        to stoichiometric coefficients (integers).

    Returns
    -------
    2 dimensional array of shape (len(substances), len(stoichs))

    """
    A = np.zeros((len(substances), len(stoichs)), dtype=int)
    for ri, sb in enumerate(substances):
        for ci, (reac, prod) in enumerate(stoichs):
            A[ri, ci] = prod.get(sb, 0) - reac.get(sb, 0)
    return A


def decompose_yields(yields, rxns, atol=1e-10):
    """ Decomposes yields into mass-action reactions

    This function offers a way to express a reaction with non-integer
    stoichiometric coefficients as a linear combination of production reactions
    with integer coefficients.

    Ak = y

    A is (n_species x n_reactions) matrix, k is "rate coefficient", y is yields


    Parameters
    ----------
    yields : OrderedDict
        Specie names as keys and yields as values.
    rxns : iterable :class:`Reaction` instances
        Dict keys must match those of ``yields`` each pair
        of dictionaries gives stoichiometry
        (1st is reactant, 2nd is products).
    atol : float
        Absolute tolerance for residuals.


    Examples
    --------
    >>> from chempy import Reaction
    >>> h2a = Reaction({'H2O': 1}, {'H2': 1, 'O': 1})
    >>> h2b = Reaction({'H2O': 1}, {'H2': 1, 'H2O2': 1}, inact_reac={'H2O': 1})
    >>> decompose_yields({'H2': 3, 'O': 2, 'H2O2': 1}, [h2a, h2b])
    array([2., 1.])

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
    from chempy import ReactionSystem
    # Sanity check:
    rxn_keys = set.union(*(rxn.keys() for rxn in rxns))
    for key in yields.keys():
        if key not in rxn_keys:
            raise ValueError("Substance key: %s not in reactions" % key)
    rsys = ReactionSystem(rxns, rxn_keys)
    A = rsys.net_stoichs(yields.keys())
    b = list(yields.values())
    unit = unit_of(b[0])
    x, residuals, rank, s = np.linalg.lstsq(
        np.asarray(A.T, dtype=np.float64), to_unitless(b, unit), rcond=2e-16*max(A.shape))
    if len(residuals) > 0:
        if np.any(residuals > atol):
            raise ValueError("atol not satisfied")
    return x*unit
