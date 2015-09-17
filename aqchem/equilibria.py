from __future__ import division, absolute_import

from functools import reduce
from operator import mul

import numpy as np
from scipy.optimize import brentq

from .chemistry import Reaction, ReactionSystem


def prodexp(factors, exponents):
    return reduce(mul, [factor ** exponent for factor, exponent
                        in zip(factors, exponents)])


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


def _solve_equilibrium_coord(c0, stoich, K, activity_product=None):
    mask, = np.nonzero(stoich)
    stoich_m = stoich[mask]
    c0_m = c0[mask]
    lower, upper = get_rc_interval(stoich_m, c0_m)
    # span = upper - lower
    return brentq(
        equilibrium_residual,
        lower,  # + delta_frac*span,
        upper,  # - delta_frac*span,
        (c0_m, stoich_m, K, activity_product)
    )


def solve_equilibrium(c0, stoich, K, activity_product=None): #, delta_frac=1e-16):
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
    stoich = np.array(stoich)
    c0 = np.array(c0)
    rc = _solve_equilibrium_coord(c0, stoich, K, activity_product)
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

class EqSystemBase(ReactionSystem):

    def __init__(self, *args, **kwargs):
        super(EqSystemBase, self).__init__(*args, **kwargs)
        self.stoichs = np.array([eq.net_stoich(self.substances)
                                 for eq in self.rxns]).transpose()

    def equilibrium_quotients(self, concs):
        return [equilibrium_quotient(concs, self.stoichs[:, ri])
                for ri in range(self.nr)]

    def charge_balance(self, concs):
        return charge_balance(self.substances, concs)

    def charge_balance_vector(self):
        return [s.charge for s in self.substances]

    def atom_balance(self, concs, atom_number):
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

    def qk(self, concs, scaling=1.0):
        return [q - scaling**rxn.dimensionality()*rxn.params*rxn.solid_factor(
            self.substances, concs) for q, rxn in zip(
                self.equilibrium_quotients(concs), self.rxns)]

    def multiple_root(self, init_concs, varied, values, **kwargs):
        nval = len(values)
        if isinstance(init_concs, dict):
            init_concs = np.array([init_concs[k] for k in self.substances])
        else:
            init_concs = np.asarray(init_concs)

        if init_concs.ndim == 1:
            init_concs = np.tile(init_concs, (nval, 1))
        elif init_concs.ndim == 2:
            if init_concs.shape[0] != nval:
                raise ValueError("Illegal dimension of init_conc")
            if init_concs.shape[1] != self.ns:
                raise ValueError("Illegal dimension of init_conc")
        else:
            raise ValueError("init_conc has too many dimensions")

        if not isinstance(varied, int):
            varied = self.substances.index(varied)

        x = np.empty((nval, self.ns))
        success = np.empty(nval, dtype=bool)
        for idx in range(nval):
            x0 = init_concs[idx, :]
            x0[varied] = values[idx]
            resx, res = self.root(x0, **kwargs)
            success[idx] = resx is not None
            x[idx, :] = resx
        return x, success

    def plot(self, init_concs, varied, values, ax=None, fail_vline=True,
             plot_kwargs=None, **kwargs):
        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.subplot(1, 1, 1, xscale='log', yscale='log')
        if plot_kwargs is None:
            plot_kwargs = {}
        x, success = self.multiple_root(init_concs, varied, values, **kwargs)
        if not isinstance(varied, int):
            varied = self.substances.index(varied)
        for idx_s in range(self.ns):
            if idx_s == varied:
                continue
            ax.plot(values[success], x[success, idx_s],
                    label=self.substances[idx_s].name, **plot_kwargs)
        ax.legend(loc='best')
        if fail_vline:
            for i, s in enumerate(success):
                if not s:
                    ax.axvline(values[i], c='k', ls='--')
        return x, success


class EqSystem(EqSystemBase):

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
        return rref.tolist()[:len(pivot)], pivot

    def independent_atoms_from_pivot(self, pivot):
        """ atom number 0 represents charge """
        atom_vecs, atm_nrs = self.atom_balance_vectors()
        return [([0] + atm_nrs)[idx] for idx in pivot]

    def preserved(self, concs, init_concs):
        rref, pivot = self.rref()
        p0s = [dot(row, init_concs) for row in rref]
        return [dot(row, concs) - p0 for row, p0 in zip(rref, p0s)], pivot

    def f(self, concs, init_concs, scaling=1, reduced=None):
        pres, pivot = self.preserved(concs, init_concs*scaling)
        if reduced is None:
            return self.qk(concs, scaling) + pres
        else:
            if reduced:
                import sympy as sp
                subs = []
                for idx, p in zip(pivot, pres)[::-1]:
                    c = concs[idx]
                    subs.append((c, (c-p.subs(subs)).simplify()))
                qk = [expr.subs(subs) for expr in self.qk(concs, scaling)]
                reduced_cbs = [sp.lambdify(
                    [y for idx, y in enumerate(concs) if idx not in pivot],
                    expr) for c, expr in subs]
                return qk, pivot, reduced_cbs
            else:
                return self.qk(concs, scaling) + pres, None, None

    def num_cb_factory(self, init_concs, jac=False, scaling=1.0, logC=False,
                       square=False, reduced=None):
        import sympy as sp
        y = sp.symarray('y', self.ns)
        f, elim, red_cbs = self.f(y, init_concs, scaling,
                                  reduced=bool(reduced))
        if reduced is not None:
            if reduced:
                y = [y[idx] for idx in range(len(y)) if idx not in elim]
        if square:
            f = [_.subs([(yi, yi*yi) for yi in y]) for _ in f]
        if logC:
            f = [_.subs([(yi, sp.exp(yi)) for yi in y]) for _ in f]
        f_lmbd = sp.lambdify(y, f, modules='numpy')

        def f_cb(arr):
            return f_lmbd(*arr)
        j_cb = None
        if jac:
            j = sp.Matrix(1, len(y), lambda _, q: f[q]).jacobian(y)
            j_lmbd = sp.lambdify(y, j, modules='numpy')

            def j_cb(arr):
                return j_lmbd(*arr)
        if reduced is None:
            return f_cb, j_cb
        else:
            return f_cb, j_cb, elim, red_cbs

    def root(self, init_concs, scaling=1.0, logC=False, square=False,
             delta=None, reduced=False, init_iter=20, **kwargs):
        import numpy as np
        from scipy.optimize import root
        init_concs = np.array([init_concs[k] for k in self.substances] if
                              isinstance(init_concs, dict) else init_concs)
        f, j, elim, red_cbs = self.num_cb_factory(
            init_concs, jac=True, scaling=scaling, logC=logC,
            square=square, reduced=reduced)
        if delta is None:
            delta = kwargs.get('tol', 1e-12)
        x0 = self.initial_guess(init_concs + delta, scaling=scaling,
                                repetitions=init_iter)
        if reduced:
            x0 = np.array([x for idx, x in enumerate(x0) if idx not in elim])
        if square:
            x0 = x0**0.5
        if logC:
            x0 = np.log(x0)
        result = root(f, x0, jac=j, **kwargs)
        if result.success:
            x = result.x
            if logC:
                x = np.exp(x)
            if square:
                x *= x
            if reduced:
                idx_red = 0
                idx_elm = 0
                new_x = []
                print('elim: ', elim)
                for idx in range(len(init_concs)):
                    if idx in elim:
                        new_x.append(red_cbs[idx_elm](*x))
                        idx_elm += 1
                    else:
                        new_x.append(x[idx_red])
                        idx_red += 1
                x = np.array(new_x)
            if np.any(x < 0):
                return None, result
            x /= scaling
            return x, result
        else:
            return None, result


    def initial_guess(self, init_concs, scaling=1, repetitions=50, steffensens=0):
        # Fixed point iteration
        def f(x):
            xnew = np.zeros_like(x)
            for ri, eq in enumerate(self.rxns):
                xnew += solve_equilibrium(
                    x, self.stoichs[:, ri],
                    eq.params * scaling**eq.dimensionality())
            return xnew/self.nr
        init_concs *= scaling
        for _ in range(repetitions):
            init_concs = f(init_concs)

        for _ in range(steffensens):
            # Steffensen's method
            f_x = f(init_concs)
            f_xf = f(init_concs + f_x)
            g = (f_xf - f_x)/f_x
            f_over_g = f_x / g
            factor = np.max(np.abs(f_over_g/init_concs))
            if factor > 0.99:
                init_concs = init_concs - f_over_g/(1.1*factor)
            else:
                init_concs = init_concs - f_over_g

        return init_concs


class REqSystem(EqSystem):

    def _c(self, coords, init_concs, scaling):
        return [scaling*ic + sum([coords[ri]*self.stoichs[cidx, ri] for ri in range(self.nr)])
                for cidx, ic in enumerate(init_concs)]

    def f(self, coords, init_concs, scaling=1):
        return self.qk(self._c(coords, init_concs, scaling), scaling)

    def num_cb_factory(self, init_concs, jac=False, scaling=1.0):
        import sympy as sp
        r = sp.symarray('r', self.nr)
        f = self.f(r, init_concs, scaling)
        f_lmbd = sp.lambdify(r, f, modules='numpy')
        def f_cb(arr):
            return f_lmbd(*arr)
        j_cb = None
        if jac:
            j = sp.Matrix(1, len(r), lambda _, q: f[q]).jacobian(r)
            j_lmbd = sp.lambdify(r, j, modules='numpy')
            def j_cb(arr):
                return j_lmbd(*arr)
        return f_cb, j_cb

    def initial_guess(self, init_concs, scaling, repetitions=10):
        def f(x):
            c0 = init_concs[:]
            xnew = np.zeros_like(x)
            for ri, eq in enumerate(self.rxns):
                xnew = _solve_equilibrium_coord(c0, self.stoichs[:, ri], eq.params)
                xnew = (x + xnew)/2
                c0 = c0 + np.dot(self.stoichs, xnew)
            return xnew
        x0 = np.zeros(self.nr)
        for idx in range(repetitions):
            print(x0) #DEBUG
            x0 = (idx*x0 + f(x0))/(idx + 1)
        print(x0) #DEBUG
        return x0*0

    def root(self, init_concs, scaling=1.0, **kwargs):
        import numpy as np
        from scipy.optimize import root
        init_concs = np.array([init_concs[k] for k in self.substances] if
                              isinstance(init_concs, dict) else init_concs)
        f, j = self.num_cb_factory(init_concs, jac=True, scaling=scaling)
        x0 = self.initial_guess(init_concs, scaling)
        result = root(f, x0, jac=j, **kwargs)
        if result.success:
            x = result.x
            C = np.array(self._c(x, init_concs, scaling))
            if np.any(C < 0):
                return None, result
            C /= scaling
            return C, result
        else:
            return None, result
