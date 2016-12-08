"""
This module collects special cases of systems for which analytic solutions
are available.
"""
from __future__ import (absolute_import, division, print_function)

from .._util import get_backend

# TODO
# ----
#
# - Add documentation
# - Rename esoteric parameter names
# - Derive more general expressions (e.g. which allows finite initial dimer concentration)


def dimerization_irrev(t, kf, initial_C, P0=1, t0=0):
    return 1/(1/initial_C + 2*kf*(t-t0))


def pseudo_irrev(t, kf, P0, t0, excess_C, limiting_C, eps_l, backend=None):
    be = get_backend(backend)
    return P0*eps_l*limiting_C*(1 - be.exp(-excess_C*kf*(t-t0)))
pseudo_irrev.name = 'Pseudo first order irreversible'


def pseudo_rev(t, kf, P0, t0, excess_C, limiting_C, eps_l, beta, backend=None):
    be = get_backend(backend)
    kb = kf/beta
    return P0*eps_l*limiting_C*excess_C*kf/(excess_C*kf + kb)*(
        1 - be.exp(-(excess_C*kf+kb)*(t-t0)))
pseudo_rev.name = 'Pseudo first order reversible'


def binary_irrev(t, kf, P0, t0, excess_C, limiting_C, eps_l, backend=None):
    be = get_backend(backend)
    return P0*eps_l*excess_C*(1 - be.exp(-kf*(excess_C-limiting_C)*(t-t0)))/(
        excess_C/limiting_C - be.exp(-kf*(t-t0)*(excess_C-limiting_C)))
binary_irrev.name = 'Second order irreversible'


def binary_rev(t, kf, kb, prod, major, minor, backend=None):
    """ Calculates the analytic transient concentration of a complex.

    From second order reversible kinetics.

    Parameters
    ----------
    t : float, Symbol or array_like
        time
    kf : number or Symbol
        Forward rate constant.
    kb : number or Symbol
        Backward rate constant.
    prod : number or Symbol
        Initial concentration of the complex
    major : number or Symbol
        Initial concentration of the more abundant reactant
    minor : number or Symbol
        Initial concentration of the less abundant reactant
    backend : module or str
        Default is 'numpy', can also be e.g. ``sympy``.

    """
    # see _integrated.ipynb for derivation
    be = get_backend(backend)
    X, Y, Z = prod, major, minor
    x0 = Y*kf
    x1 = Z*kf
    x2 = kb + x0 + x1
    x3 = 2*X*kf
    x4 = -kb - x0 - x1
    x5 = -x3 + x4
    x6 = be.sqrt(-4*kf*(X**2*kf + X*x0 + X*x1 + Z*x0) + x5**2)
    x7 = -x6
    x8 = (x4 + x7)*be.exp(t*x6)
    return (-x5*x8 - x6*x8 - (x4 + x6)*(x2 + x3 + x6))/(2*kf*(x2 + x7 + x8))
binary_rev.name = 'Second order reversible'
