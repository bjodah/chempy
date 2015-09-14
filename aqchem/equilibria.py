from __future__ import division

import numpy as np
from scipy.optimize import brentq


def equilibrium_quotient(c, stoich):
    tot = 1
    for idx, nr in enumerate(stoich):
        tot *= c[idx]**nr
    return tot


def equilibrium_residual(rc, c, stoich, K, activity_product=None):
    """
    Parameters
    ---------
    rc: float
        Reaction coordinate
    c: array_like of reals
        concentrations
    stoich: tuple
        per specie stoichiometry coefficient
    K: float
        equilibrium constant
    activity_product: callable
        callback for calculating the activity product taking
        concentration as single parameter.
    """
    c = c0 + rc*stoich
    Q = equilibrium_quotient(c, stoich)
    if activity_product is not None:
        Q *= activity_product(c)
    return K - Q


def get_rc_interval(stoich, c0):
    """ get reaction coordinate interval """
    limits = c0/stoich
    if np.any(limits < 0):
        upper = -np.max(limits[np.argwhere(limits < 0)])
    else:
        upper = 0

    if np.any(limits > 0):
        lower = -np.min(limits[np.argwhere(limits > 0)])
    else:
        lower = 0

    if lower is 0 and upper is 0:
        raise ValueError("0-interval")
    else:
        return lower, upper


def solve_equilibrium(c0, stoich, K, activity_product=None, delta_frac=1e-16):
    """
    Solve equilibrium concentrations by using scipy.optimize.brentq

    Parameters
    ----------
    c0: array_like
        Initial guess of equilibrium concentrations
    stoich: tuple
        per specie stoichiometry coefficient (law of mass action)
    K: float
        equilibrium constant
    activity_product: callable
        see ``equilibrium_residual``
    delta_frac: float
        to avoid division by zero the span of searched values for
        the reactions coordinate (rc) is shrunk by 2*delta_frac*span(rc)
    """
    lower, upper = get_rc_interval(np.array(stoich), np.array(c0))
    span = upper - lower
    rc = brentq(
        equilibrium_residual,
        lower + delta_frac*span,
        upper - delta_frac*span,
        (np.array(c0), np.array(stoich), K, activity_product)
    )
    return np.array(c0) + rc*np.array(stoich)
