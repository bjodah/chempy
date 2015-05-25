from __future__ import division


def dimerization_irrev(t, kf, P0, t0, initial_C):
    pass


def pseudo_irrev(t, kf, P0, t0, excess_C, limiting_C, eps_l, exp=None):
    if exp is None:
        try:
            from numpy import exp
        except ImportError:
            from math import exp
    return P0*eps_l*limiting_C*(1 - exp(-excess_C*kf*(t-t0)))
pseudo_irrev.name = 'Pseudo first order irreversible'


def pseudo_rev(t, kf, P0, t0, excess_C, limiting_C, eps_l, beta, exp=None):
    if exp is None:
        try:
            from numpy import exp
        except ImportError:
            from math import exp
    kb = kf/beta
    return P0*eps_l*limiting_C*excess_C*kf/(excess_C*kf + kb)*(
        1 - exp(-(excess_C*kf+kb)*(t-t0)))
pseudo_rev.name = 'Pseudo first order reversible'


def binary_irrev(t, kf, P0, t0, excess_C, limiting_C, eps_l, exp=None):
    if exp is None:
        try:
            from numpy import exp
        except ImportError:
            from math import exp
    return P0*eps_l*excess_C*(1 - exp(-kf*(excess_C-limiting_C)*(t-t0)))/(
        excess_C/limiting_C - exp(-kf*(t-t0)*(excess_C-limiting_C)))
binary_irrev.name = 'Second order irreversible'


def binary_rev(t, kf, P0, t0, excess_C, limiting_C, eps_l, beta,
               exp=None, one=1):
    if exp is None:
        try:
            from numpy import exp
        except ImportError:
            from math import exp
    kb = kf/beta
    a = kf
    b = -excess_C*kf - limiting_C*kf - kb
    c = excess_C*limiting_C*kf
    P = (b**2 - 4*a*c)**(one/2)
    Q = P + b
    R = P - b
    return P0*eps_l*Q*(1 - exp(P*(t-t0)))/(2*a*(Q/R + exp(P*(t-t0))))
binary_rev.name = 'Second order reversible'
