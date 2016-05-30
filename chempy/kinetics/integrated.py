"""
This module collects special cases of systems for which analytic solutions
are available.
"""
from __future__ import (absolute_import, division, print_function)

from .._util import get_backend

# Add documentation
# Rename esoteric parameter names


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


def binary_rev(t, kf, P0, t0, excess_C, limiting_C, eps_l, beta, backend=None):
    be = get_backend(backend)
    one = backend.pi**0
    kb = kf/beta
    a = kf
    b = -excess_C*kf - limiting_C*kf - kb
    c = excess_C*limiting_C*kf
    P = (b**2 - 4*a*c)**(one/2)
    Q = P + b
    R = P - b
    return P0*eps_l*Q*(1 - be.exp(P*(t-t0)))/(2*a*(Q/R + be.exp(P*(t-t0))))
binary_rev.name = 'Second order reversible'
