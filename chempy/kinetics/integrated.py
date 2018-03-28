"""
This module collects a few analytic solutions of initial value problems (IVPs) in
chemical kinetics (i.e. integrated rate expressions in closed form).
The expressions are useful in e.g. regression or for comparison with numerical
solution of the corresponding ODE system.
"""
from __future__ import (absolute_import, division, print_function)

from .._util import get_backend


def dimerization_irrev(t, kf, initial_C, P0=1, t0=0):
    return 1/(1/initial_C + 2*kf*(t-t0))


def pseudo_irrev(t, kf, prod, major, minor, backend=None):
    """ Analytic product transient of a irreversible pseudo first order reaction.

    Product concentration vs time from pseudo-first order irreversible kinetics.

    Parameters
    ----------
    t : float, Symbol or array_like
        Time.
    kf : number or Symbol
        Forward (bimolecular) rate constant.
    kb : number or Symbol
        Backward (unimolecular) rate constant.
    prod : number or Symbol
        Initial concentration of the complex.
    major : number or Symbol
        Initial concentration of the more abundant reactant.
    minor : number or Symbol
        Initial concentration of the less abundant reactant.
    backend : module or str
        Default is 'numpy', can also be e.g. ``sympy``.

    """
    be = get_backend(backend)
    return prod + minor*(1 - be.exp(-major*kf*t))
pseudo_irrev.name = 'Pseudo first order irreversible'


def pseudo_rev(t, kf, kb, prod, major, minor, backend=None):
    """ Analytic product transient of a reversible pseudo first order reaction.

    Product concentration vs time from pseudo-first order reversible kinetics.

    Parameters
    ----------
    t : float, Symbol or array_like
        Time.
    kf : number or Symbol
        Forward (bimolecular) rate constant.
    kb : number or Symbol
        Backward (unimolecular) rate constant.
    prod : number or Symbol
        Initial concentration of the complex.
    major : number or Symbol
        Initial concentration of the more abundant reactant.
    minor : number or Symbol
        Initial concentration of the less abundant reactant.
    backend : module or str
        Default is 'numpy', can also be e.g. ``sympy``.

    """
    be = get_backend(backend)
    return (
        -kb*prod + kf*major*minor + (kb*prod - kf*major*minor) *
        be.exp(-t*(kb + kf*major))
    )/(kb + kf*major)
pseudo_rev.name = 'Pseudo first order reversible'


def binary_irrev(t, kf, prod, major, minor, backend=None):
    """ Analytic product transient of a irreversible 2-to-1 reaction.

    Product concentration vs time from second order irreversible kinetics.

    Parameters
    ----------
    t : float, Symbol or array_like
    kf : number or Symbol
        Forward (bimolecular) rate constant.
    prod : number or Symbol
        Initial concentration of the complex.
    major : number or Symbol
        Initial concentration of the more abundant reactant.
    minor : number or Symbol
        Initial concentration of the less abundant reactant.
    backend : module or str
        Default is 'numpy', can also be e.g. ``sympy``.

    """
    be = get_backend(backend)
    return prod + major*(1 - be.exp(-kf*(major-minor)*t))/(
        major/minor - be.exp(-kf*t*(major-minor)))
binary_irrev.name = 'Second order irreversible'


def binary_rev(t, kf, kb, prod, major, minor, backend=None):
    """ Analytic product transient of a reversible 2-to-1 reaction.

    Product concentration vs time from second order reversible kinetics.

    Parameters
    ----------
    t : float, Symbol or array_like
        Time.
    kf : number or Symbol
        Forward (bimolecular) rate constant.
    kb : number or Symbol
        Backward (unimolecular) rate constant.
    prod : number or Symbol
        Initial concentration of the complex.
    major : number or Symbol
        Initial concentration of the more abundant reactant.
    minor : number or Symbol
        Initial concentration of the less abundant reactant.
    backend : module or str
        Default is 'numpy', can also be e.g. ``sympy``.

    """
    # see _integrated.ipynb for derivation
    be = get_backend(backend)
    X, Y, Z = prod, major, minor
    x0 = Y*kf
    x1 = Z*kf
    x2 = 2*X*kf
    x3 = -kb - x0 - x1
    x4 = -x2 + x3
    x5 = be.sqrt(-4*kf*(X**2*kf + X*x0 + X*x1 + Z*x0) + x4**2)
    x6 = kb + x0 + x1 + x5
    x7 = (x3 + x5)*be.exp(-t*x5)
    x8 = x3 - x5
    return (x4*x8 + x5*x8 + x7*(x2 + x6))/(2*kf*(x6 + x7))

binary_rev.name = 'Second order reversible'


def unary_irrev_cstr(t, k, r, p, fr, fp, fv, backend=None):
    """ Analytic solution for ``A -> B`` in a CSTR.

    Analytic solution for a first order process in a continuously
    stirred tank reactor (CSTR).

    Parameters
    ----------
    t : array_like
    k : float_like
        Rate constant
    r : float_like
        Initial concentration of reactant.
    p : float_like
        Initial concentration of product.
    fr : float_like
        Concentration of reactant in feed.
    fp : float_like
        Concentration of product in feed.
    fv : float_like
        Feed rate / tank volume ratio.
    backend : module or str
        Default is 'numpy', can also be e.g. ``sympy``.

    Returns
    -------
    length-2 tuple
        concentrations of reactant and product

    """
    # See _kinetics_cstr.ipynb
    be = get_backend(backend)
    x0 = fr*fv
    x1 = fv + k
    x2 = 1/x1
    x3 = fv*r + k*r - x0
    x4 = fr*k
    x5 = be.exp(-fv*t)
    return (
        x0*x2 + x2*x3*be.exp(-t*x1),
        -x2*x3*x5*(-1 + be.exp(-k*t)) + x2*x5*(-fp*fv - fp*k + fv*p + k*p - x4) + x2*(fp*x1 + x4)
    )


def binary_irrev_cstr(t, k, r, p, fr, fp, fv, n=1, backend=None):
    """ Analytic solution for ``2 A -> n B`` in a CSTR.

    Parameters
    ----------
    t : array_like
    k : float_like
        Rate constant
    r : float_like
        Initial concentration of reactant.
    p : float_like
        Initial concentration of product.
    fr : float_like
        Concentration of reactant in feed.
    fp : float_like
        Concentration of product in feed.
    fv : float_like
        Feed rate / tank volume ratio.
    n : int
    backend : module or str
        Default is 'numpy', can also be e.g. ``sympy``.

    Returns
    -------
    length-2 tuple
        concentrations of reactant and product

    """
    # Mathematica source:
    # FortranForm[
    # DSolve[{Derivative[1][A][t] == a0 f - f A[t] - 2 k A[t]^2,
    #   Derivative[1][B][t] == b0 f + n*k A[t]^2 - f B[t], A[0] == x,
    #   B[0] == y}, {A[t], B[t]}, {t}]]
    # Post processed using sympy's cse function.
    # (see _derive_analytic_cstr_bireac.ipynb)
    be = get_backend(backend)
    atanh = getattr(be, 'atanh', be.arctanh)
    three = 3*be.cos(0)

    x0 = 1/k
    x1 = be.sqrt(fv)
    x2 = 8*k
    x3 = fr*x2
    x4 = be.sqrt(fv + x3)
    x5 = x1*x4
    x6 = x1*x4/2
    x7 = atanh((-fv**(three/2)*x4 - 4*k*r*x5)/(fv**2 + fv*x3))
    x8 = fv*t
    x9 = fp*x2
    x10 = 4*k*n
    x11 = fr*x10
    x12 = be.exp(x8)
    x13 = n*x12
    return (
        x0*(-fv + x5*be.tanh(t*x6 - x7))/4,
        x0*(fv*x13 + 8*k*p + r*x10 - x1*x13*x4*be.tanh(
            x6*(t - 2*x7/(x1*x4))
        ) + x11*x12 - x11 + x12*x9 - x9)*be.exp(-x8)/8
    )
