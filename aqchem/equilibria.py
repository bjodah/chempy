from __future__ import division, absolute_import

from functools import reduce
from operator import mul

import numpy as np
from scipy.optimize import brentq

from .chemistry import Reaction, ReactionSystem

def prodexp(factors, exponents):
    return reduce(mul, [factor ** exponent for factor, exponent in zip(factors, exponents)])


def equilibrium_quotient(c, stoich):
    tot = 1
    for idx, nr in enumerate(stoich):
        tot *= c[idx]**nr
    return tot


def equilibrium_residual(rc, c0, stoich, K, activity_product=None):
    """
    Parameters
    ---------
    rc: float
        Reaction coordinate
    c0: array_like of reals
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
    stoich = np.asarray(stoich)
    c0 = np.asarray(c0)
    mask = np.argwhere(stoich != 0)
    stoich_m = stoich[mask]
    c0_m = c0[mask]
    lower, upper = get_rc_interval(np.array(stoich_m), np.array(c0_m))
    span = upper - lower
    if True:
        rc = brentq(
            equilibrium_residual,
            lower, # + delta_frac*span,
            upper, # - delta_frac*span,
            (c0_m, stoich_m, K, activity_product)
        )
    else:
        rc = brentq(
            equilibrium_residual,
            lower + delta_frac*span,
            upper - delta_frac*span,
            (c0_m, stoich_m, K, activity_product)
        )
    return c0 + rc*stoich


class Equilibrium(Reaction):
    """
    Represents equilibrium reaction, subclass of Reatcion.
    Extra attributes:
    - solubility_product: bool (default: False)
    """


    latex_arrow = '\\rightleftharpoons'

    # def __init__(self, *args, **kwargs):
    #     self.solubility_product = kwargs.pop('solubility_product', False)
    #     super(self.__class__, self).__init__(*args, **kwargs)

    def K(self, T=None):
        if callable(self.params):
            if T is None:
                raise ValueError("No T provided")
            return self.params(T)
        else:
            if T is not None:
                raise ValueError("T provided but params not callble")
            return self.params

    def solid_factor(self, substances, concs):
        factor = 1
        for r, n in self.reac.items():
            if r.solid:
                factor *= concs[substances.index(r)]**-n
        for p, n in self.prod.items():
            if p.solid:
                factor *= concs[substances.index(p)]**n
        return factor

    def dimensionality(self):
        result = 0
        for r, n in self.reac.items():
            if r.solid:
                continue
            result -= n
        for p, n in self.prod.items():
            if p.solid:
                continue
            result += n
        return result


def atom_balance(substances, concs, atom_number):
    res = 0
    for s, c in zip(substances, concs):
        res += s.elemental_composition.get(atom_number, 0)*c
    return res

def charge_balance(substances, concs):
    res = 0
    for s, c in zip(substances, concs):
        res += s.charge*c
    return res


def dot(iter_a, iter_b):
    return sum([a*b for a, b in zip(iter_a, iter_b)])


class EqSystem(ReactionSystem):

    def equilibrium_quotients(self, c):
        return [equilibrium_quotient(
            c, eq.net_stoich(
                self.substances)) for eq in self.rxns]

    def charge_balance(self, concs):
        return charge_balance(self.substances, concs)

    def charge_balance_vector(self):
        return [s.charge for s in self.substances]

    def atom_balance(self, concs, atom_number):
        # Convenience method
        return atom_balance(self.substances, concs, atom_number)

    def atom_balance_vectors(self):
        atom_numbers = set()
        for s in self.substances:
            for key in s.elemental_composition:
                atom_numbers.add(key)
        vectors = []
        sorted_atom_numbers = sorted(atom_numbers)
        for atm_nr in sorted_atom_numbers:
            vectors.append([s.elemental_composition.get(atm_nr, 0) for
                            s in self.substances])
        return vectors, sorted_atom_numbers

    def rref(self):
        """
        Calculates the Reduced Row Echelon Form of the linear equations
        system for charge and atom balance.

        Returns
        -------
        Length-2 tuple:
           - list of lists with coefficients for the substances
           - pivot array
        """
        from sympy import Matrix
        atom_vecs, atm_nrs = self.atom_balance_vectors()
        M = Matrix([self.charge_balance_vector()] + atom_vecs)
        rref, pivot = M.rref()
        #indep_atm_nrs =
        return rref.tolist()[:len(pivot)], pivot #indep_atm_nrs

    def independent_atoms_from_pivot(self, pivot):
        """ atom number 0 represents charge """
        atom_vecs, atm_nrs = self.atom_balance_vectors()
        return [([0] + atm_nrs)[idx] for idx in pivot]

    def qk(self, concs, scaling=1.0):
        return [q - scaling**rxn.dimensionality()*rxn.params*rxn.solid_factor(
            self.substances, concs) for q, rxn in zip(
                self.equilibrium_quotients(concs), self.rxns)]

    def preserved(self, concs, init_concs):
        rref, pivot = self.rref()
        p0s = [dot(row, init_concs) for row in rref]
        return [dot(row, concs) - p0 for row, p0 in zip(rref, p0s)], pivot

    def f(self, concs, init_concs, scaling=1, reduced=None):
        pres, pivot = self.preserved(concs, init_concs)
        if reduced is None:
            return self.qk(concs, scaling) + pres
        else:
            if reduced:
                subs = []
                for idx, p in zip(pivot, pres)[::-1]:
                    c = concs[idx]
                    subs.append((c, (c-p.subs(subs)).simplify()))
                qk = [expr.subs(subs) for expr in self.qk(concs, scaling)]
                return qk, pivot
            else:
                return self.qk(concs, scaling) + pres, None

    def num_cb_factory(self, init_concs, jac=False, scaling=1.0, logC=False,
                       square=False, reduced=None):
        import sympy as sp
        y = sp.symarray('y', len(self.substances))
        f, elim = self.f(y, init_concs, scaling, reduced=bool(reduced))
        if reduced is not None:
            if reduced:
                y = [y[idx] for idx in range(len(y)) if idx not in elim]
                print(y)
        if square:
            f = [_.subs([(yi, yi*yi) for yi in y]) for _ in f]
            print(f) ## DEBUG
        if logC:
            f = [_.subs([(yi, sp.exp(yi)) for yi in y]) for _ in f]
        f_lmbd = sp.lambdify(y, f, modules='numpy')
        def f_cb(arr):
            print(arr) ## DEBUG
            return f_lmbd(*arr)
        j_cb = None
        if jac:
            print(len(f), len(y))
            j = sp.Matrix(1, len(y), lambda _, q: f[q]).jacobian(y)
            print(j) ## DEBUG
            j_lmbd = sp.lambdify(y, j, modules='numpy')
            def j_cb(arr):
                return j_lmbd(*arr)
        if reduced is None:
            return f_cb, j_cb
        else:
            return f_cb, j_cb, elim

    def root(self, init_concs, scaling=1.0, logC=False, square=False,
             delta=None, reduced=False, **kwargs):
        import numpy as np
        from scipy.optimize import root
        init_concs = np.asarray([init_concs[k] for k in self.substances] if
                             isinstance(init_concs, dict) else init_concs)
        f, j, elim = self.num_cb_factory(init_concs*scaling, jac=True,
                                         scaling=scaling, logC=logC,
                                         square=square, reduced=reduced)
        if delta is None:
            delta = kwargs.get('tol', 1e-12)
        x0 = self.initial_guess(init_concs*scaling + delta)
        if reduced:
            x0 = np.array([x for idx, x in enumerate(x0) if idx not in elim])
        if square:
            x0 = x0**0.5
        if logC:
            x0 = np.log(x0)
        result = root(f, x0, jac=j)
        if result.success:
            # TODO: back-transform result
            x = result.x
            if logC:
                x = np.exp(x)
            if square:
                x *= x
            x /= scaling
            print('Success: ', x)
        if reduced:
            return result # TODO: reintroduce eliminated...
        else:
            return result

    def initial_guess(self, init_concs, weights=range(1, 5)):
        for w in weights:
            for eq in self.rxns:
                new_init_concs = solve_equilibrium(
                    init_concs, eq.net_stoich(self.substances), eq.params)
                init_concs = (w*init_concs + new_init_concs)/(w+1)
        return init_concs
