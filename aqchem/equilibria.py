from __future__ import division, absolute_import

import copy
from collections import defaultdict
from functools import reduce
from operator import mul

import numpy as np

from .chemistry import Reaction, ReactionSystem


def prodexp(factors, exponents):
    return reduce(mul, [factor ** exponent for factor, exponent
                        in zip(factors, exponents)])


def equilibrium_quotient(concs, stoich):
    if not hasattr(concs, 'ndim') or concs.ndim == 1:
        tot = 1
    else:
        tot = np.ones(concs.shape[0])
        concs = concs.T
    for nr, conc in zip(stoich, concs):
        tot *= conc**nr
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
    from scipy.optimize import brentq
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


def solve_equilibrium(c0, stoich, K, activity_product=None):
    # , delta_frac=1e-16):
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
                raise ValueError("T provided but params not callable")
            return self.params

    def solid_factor(self, substances, sc_concs):
        factor = 1
        for r, n in self.reac.items():
            if r.solid:
                factor *= sc_concs[substances.index(r)]**-n
        for p, n in self.prod.items():
            if p.solid:
                factor *= sc_concs[substances.index(p)]**n
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
    if not hasattr(concs, 'ndim') or concs.ndim == 1:
        res = 0
    elif concs.ndim == 2:
        res = np.zeros(concs.shape[0])
        concs = concs.T
    else:
        raise NotImplementedError
    for s, c in zip(substances, concs):
        res += s.charge*c
    return res


def dot(iter_a, iter_b):
    return sum([a*b for a, b in zip(iter_a, iter_b)])


class EqSystemBase(ReactionSystem):

    def as_per_substance_array(self, cont):
        if isinstance(cont, np.ndarray):
            pass
        elif isinstance(cont, dict):
            cont = [cont[k] for k in self.substances]
        cont = np.asarray(cont)
        if cont.shape[-1] != self.ns:
            raise ValueError("Incorrect size")
        return cont

    def as_substance_index(self, sbstnc):
        if isinstance(sbstnc, int):
            return sbstnc
        else:
            return self.substances.index(sbstnc)

    def __init__(self, *args, **kwargs):
        super(EqSystemBase, self).__init__(*args, **kwargs)

    @property
    def stoichs(self):
        return np.array([eq.net_stoich(self.substances)
                         for eq in self.rxns]).transpose()


    def upper_conc_bounds(self, init_concs):
        init_concs_arr = self.as_per_substance_array(init_concs)
        atom_conc = defaultdict(float)
        for conc, subst in zip(init_concs_arr, self.substances):
            for atm_nr, coeff in subst.elemental_composition.items():
                atom_conc[atm_nr] += coeff*conc
        bounds = []
        for subst in self.substances:
            upper = float('inf')
            for atm_nr, coeff in subst.elemental_composition.items():
                upper = min(upper, atom_conc[atm_nr]/coeff)
            bounds.append(upper)
        return bounds

    def equilibrium_quotients(self, concs):
        return [equilibrium_quotient(concs, self.stoichs[:, ri])
                for ri in range(self.nr)]

    def charge_balance(self, concs):
        return charge_balance(self.substances,
                              self.as_per_substance_array(concs))

    def charge_balance_vector(self):
        return [s.charge for s in self.substances]

    def atom_balance(self, concs, atom_number):
        return atom_balance(self.substances,
                            self.as_per_substance_array(concs),
                            atom_number)

    def atom_balance_vectors(self, skip_atom_nrs=()):
        atom_numbers = set()
        for s in self.substances:
            for key in s.elemental_composition:
                atom_numbers.add(key)
        vectors = []
        sorted_atom_numbers = sorted(atom_numbers)
        for atm_nr in sorted_atom_numbers:
            if atm_nr in skip_atom_nrs:
                continue
            vectors.append([s.elemental_composition.get(atm_nr, 0) for
                            s in self.substances])
        return vectors, sorted_atom_numbers

    def atom_conservation(self, concs, init_concs):
        atom_vecs, atm_nrs = self.atom_balance_vectors()
        A = np.array(atom_vecs)
        return (atm_nrs,
                np.dot(A, self.as_per_substance_array(concs).T),
                np.dot(A, self.as_per_substance_array(init_concs).T))

    def qk(self, sc_concs, scaling=1.0, norm=False):
        vec = []
        qs = self.equilibrium_quotients(sc_concs)
        for idx, rxn in enumerate(self.rxns):
            k = (scaling**rxn.dimensionality()*rxn.params *
                 rxn.solid_factor(self.substances, sc_concs))
            if norm:
                vec.append(qs[idx]/k - 1)
            else:
                vec.append(qs[idx] - k)
        return vec

    def multiple_root(self, init_concs, varied=None, values=None, carry=False,
                      init_guess=None, **kwargs):
        if carry:
            if init_guess is not None:
                raise NotImplementedError
        nval = len(values)
        if isinstance(init_concs, dict):
            init_concs = np.array([init_concs[k] for k in self.substances])
        else:
            init_concs = np.asarray(init_concs)

        if init_concs.ndim == 1:
            init_concs = np.tile(init_concs, (nval, 1))
            init_concs[:, self.as_substance_index(varied)] = values
        elif init_concs.ndim == 2:
            if init_concs.shape[0] != nval:
                raise ValueError("Illegal dimension of init_concs")
            if init_concs.shape[1] != self.ns:
                raise ValueError("Illegal dimension of init_concs")
        else:
            raise ValueError("init_concs has too many dimensions")

        x0 = None
        x = np.empty((nval, self.ns))
        success = np.empty(nval, dtype=bool)
        res = None  # silence pyflakes
        for idx in range(nval):
            root_kw = kwargs.copy()
            g0 = init_guess if init_guess is None else init_guess[idx, :]
            if carry and idx > 0 and res.success:
                x0 = res.x
            resx, res = self.root(init_concs[idx, :], init_guess=g0,
                                  x0=x0, **root_kw)
            success[idx] = resx is not None
            x[idx, :] = resx
        return x, init_concs, success

    def plot(self, init_concs, varied, values, ax=None, fail_vline=True,
             plot_kwargs=None, subplot_kwargs=None, **kwargs):
        if ax is None:
            import matplotlib.pyplot as plt
            if subplot_kwargs is None:
                subplot_kwargs = dict(xscale='log', yscale='log')
            ax = plt.subplot(1, 1, 1, **subplot_kwargs)
        if plot_kwargs is None:
            plot_kwargs = {}
        x, new_init_concs, success = self.multiple_root(
            init_concs, varied, values, **kwargs)
        ls, c = '- -- : -.'.split(), 'krgbcmy'
        extra_kw = {}
        for idx_s in range(self.ns):
            if idx_s == self.as_substance_index(varied):
                continue
            if 'ls' not in plot_kwargs and 'linestyle' not in plot_kwargs:
                extra_kw['ls'] = ls[idx_s % len(ls)]
            if 'c' not in plot_kwargs and 'color' not in plot_kwargs:
                extra_kw['c'] = c[idx_s % len(c)]
            try:
                lbl = '$' + self.substances[idx_s].latex_name + '$'
            except TypeError:
                lbl = self.substances[idx_s].name
            extra_kw.update(plot_kwargs)
            ax.plot(values[success], x[success, idx_s],
                    label=lbl, **extra_kw)

        # ax.legend(loc='best')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        # Put a legend to the right of the current axis
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        ax.set_xlabel(str(varied))
        ax.set_ylabel('Concentration')
        if fail_vline:
            for i, s in enumerate(success):
                if not s:
                    ax.axvline(values[i], c='k', ls='--')
        return x, new_init_concs, success

    def plot_errors(self, concs, init_concs, varied, axes=None,
                    charge=True, atoms=True, Q=True,
                    subplot_kwargs=None):
        if axes is None:
            import matplotlib.pyplot as plt
            if subplot_kwargs is None:
                subplot_kwargs = dict(xscale='log')
            fig, axes = plt.subplots(1, 2, figsize=(10, 4),
                                     subplot_kw=subplot_kwargs)

        if charge:
            q1 = self.charge_balance(concs)
            q2 = self.charge_balance(init_concs)
            axes[0].plot(concs[:, self.as_substance_index(varied)],
                         q1-q2, label='Abs charge error')
            axes[1].plot(concs[:, self.as_substance_index(varied)],
                         (q1-q2)/np.abs(q2), label='Rel charge error')

        if atoms:
            atm_nrs, m1, m2 = self.atom_conservation(concs, init_concs)
            for atm_nr, a1, a2 in zip(atm_nrs, m1, m2):
                axes[0].plot(concs[:, self.as_substance_index(varied)],
                             a1-a2, label='Abs ' + str(atm_nr))
                axes[1].plot(concs[:, self.as_substance_index(varied)],
                             (a1-a2)/np.abs(a2), label='Rel ' + str(atm_nr))

        if Q:
            qs = self.equilibrium_quotients(concs)
            ks = [rxn.params*rxn.solid_factor(self.substances, concs) for
                  rxn in self.rxns]
            for idx, (q, k) in enumerate(zip(qs, ks)):
                axes[0].plot(concs[:, self.as_substance_index(varied)],
                             q-k, label='Abs R:' + str(idx))
                axes[1].plot(concs[:, self.as_substance_index(varied)],
                             (q-k)/k, label='Rel R:' + str(idx))

        axes[0].legend(loc='best')
        axes[1].legend(loc='best')


class EqSystem(EqSystemBase):

    def rref(self, charge=True, skip_atom_nrs=()):
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
        chg = [self.charge_balance_vector()] if charge else []
        M = Matrix(chg + [v for v, n in zip(atom_vecs, atm_nrs)
                          if n not in skip_atom_nrs])
        rref, pivot = M.rref()
        return rref.tolist()[:len(pivot)], pivot

    def independent_atoms_from_pivot(self, pivot):
        """ atom number 0 represents charge """
        atom_vecs, atm_nrs = self.atom_balance_vectors()
        return [([0] + atm_nrs)[idx] for idx in pivot]

    def preserved(self, sc_concs, sc_init_concs, charge=True,
                  skip_atom_nrs=(), presw=1, norm=False, rref=True):
        res = []
        if rref:
            vecs, pivot = self.rref(charge, skip_atom_nrs)
        else:
            atom_vecs, atm_nrs = self.atom_balance_vectors()
            chg = [self.charge_balance_vector()] if charge else []
            vecs = chg + [v for v, n in zip(atom_vecs, atm_nrs)
                          if n not in skip_atom_nrs]
            pivot = None
        s0s = [dot(row, sc_init_concs) for row in vecs]
        for row, s0 in zip(vecs, s0s):
            if norm:
                res.append(presw*(dot(row, sc_concs)/s0 - 1))
            else:
                res.append(presw*(dot(row, sc_concs) - s0))
        return res, pivot

    def f(self, sc_concs, init_concs, scaling=1, reduced=False, norm=False,
          pres_norm=False, pres1st=False, presw=1, const_indices=(),
          rref=True, charge=None, extra_pres_sq=False):
        sc_concs = copy.copy(sc_concs)
        sc_init_concs = scaling*init_concs
        skip_atom_nrs = set()
        if charge is None:
            charge = True
        for cidx in const_indices:
            for k in self.substances[cidx].elemental_composition:
                skip_atom_nrs.add(k)
            if self.substances[cidx].charge != 0:
                charge = False
            sc_concs[cidx] = sc_init_concs[cidx]
        qk = self.qk(sc_concs, scaling, norm)
        if pres_norm and reduced:
            raise NotImplementedError
        pres, pivot = self.preserved(
            sc_concs, sc_init_concs, charge=charge,
            skip_atom_nrs=skip_atom_nrs, presw=presw, norm=pres_norm,
            rref=rref)
        if extra_pres_sq:
            pv, _ = self.preserved(sc_concs*sc_concs,
                                   sc_init_concs*sc_init_concs,
                                   charge=charge, skip_atom_nrs=skip_atom_nrs,
                                   presw=presw, norm=pres_norm, rref=False)
            qk += pv
        if reduced:
            import sympy as sp
            subs = []
            for idx, p in zip(pivot, pres)[::-1]:
                c = sc_concs[idx]
                subs.append((c, (c-p.subs(subs)).simplify()))
            res = [expr.subs(subs) for expr in qk]
            reduced_cbs = [sp.lambdify(
                [y for idx, y in enumerate(sc_concs) if idx not in pivot],
                expr) for _, expr in subs[::-1]]
        else:
            res = (pres + qk) if pres1st else (qk + pres)
            pivot, reduced_cbs = None, None
        return res, pivot, reduced_cbs

    def num_cb_factory(self, init_concs, jac=False, scaling=1.0, logC=False,
                       square=False, reduced=False, norm=False, pres1st=False,
                       pres_norm=False, presw=1, const_indices=(), tanh_b=None,
                       rref=True, charge=None, extra_pres_sq=False):
        import sympy as sp
        y = sp.symarray('y', self.ns)
        f, elim, red_cbs = self.f(
            y, init_concs, scaling, reduced=reduced, norm=norm,
            pres1st=pres1st, pres_norm=pres_norm, presw=presw,
            const_indices=const_indices, rref=rref, charge=charge,
            extra_pres_sq=extra_pres_sq)

        if elim is not None:
            for eidx in elim:
                if eidx in const_indices:
                    raise NotImplementedError
        else:
            elim = []
        if reduced or len(const_indices) > 0:
            y = [y[idx] for idx in range(len(y)) if
                 idx not in elim and idx not in const_indices]

        if len(y) > len(f):
            raise ValueError("Under-determined system of equations")

        if tanh_b:
            tanh_subs = [(yi, m + s*sp.tanh((yi - m)/s)) for
                         yi, m, s in zip(y, tanh_b[0], tanh_b[1])]
            f = [_.subs(tanh_subs) for _ in f]
        if square:
            sq_subs = [(yi, yi**2) for yi in y]
            f = [_.subs(sq_subs) for _ in f]
        if logC:
            logC_subs = [(yi, sp.exp(yi)) for yi in y]
            f = [_.subs(logC_subs).powsimp() for _ in f]
        f_lmbd = sp.lambdify(y, f, modules='numpy')

        def f_cb(arr):
            return f_lmbd(*arr)
        j_cb = None
        if jac:
            j = sp.Matrix(1, len(f), lambda _, q: f[q]).jacobian(y)
            j_lmbd = sp.lambdify(y, j, modules='numpy')

            def j_cb(arr):
                return j_lmbd(*arr)
            return f_cb, j_cb, elim, red_cbs

    def root(self, init_concs, scaling=1.0, logC=False, square=False,
             tanh=False, delta=None, reduced=False, norm=False, init_iter=20,
             pres_norm=False, init_guess=None, x0=None, pres1st=False, presw=1,
             const=(), rref=True, charge=None, extra_pres_sq=False, **kwargs):
        from scipy.optimize import root
        init_concs = self.as_per_substance_array(init_concs)
        const_indices = list(map(self.as_substance_index, const))
        sc_upper_bounds = np.array(self.upper_conc_bounds(
            init_concs*scaling))
        if tanh:
            factor = 0.5
            sc_lower_bounds = -factor*sc_upper_bounds
            sc_bounds_span = sc_upper_bounds - sc_lower_bounds
            sc_middle = sc_lower_bounds + 0.5*sc_bounds_span
            tanh_b = (sc_middle, sc_bounds_span)
        else:
            tanh_b = None
        f, j, elim, red_cbs = self.num_cb_factory(
            init_concs, jac=True, scaling=scaling, logC=logC, square=square,
            tanh_b=tanh_b, reduced=reduced, norm=norm, pres_norm=pres_norm,
            pres1st=pres1st, presw=presw, const_indices=const_indices,
            rref=rref, charge=charge, extra_pres_sq=extra_pres_sq)
        if delta is None:
            delta = kwargs.get('tol', 1e-12)
        if x0 is None:
            if init_guess is None:
                init_guess = init_concs
            x0 = self.initial_guess(init_guess + delta, scaling=scaling,
                                    repetitions=init_iter)
            if reduced or len(const) > 0:
                x0 = np.array([x for idx, x in enumerate(x0) if
                               idx not in elim and idx not in const_indices])
            if tanh:
                # c = m + s*tanh((x - m)/s)
                # x = m + s*atanh((c - m)/s)
                x0 = sc_middle + sc_bounds_span*np.arctanh(
                    (x0 - sc_middle)/sc_bounds_span)
            if square:
                x0 = x0**0.5
            if logC:
                x0 = np.log(x0)
        else:
            if init_guess is not None:
                raise ValueError("x0 and init_guess passed")
        result = root(f, x0, jac=j, **kwargs)
        if result.success:
            x = result.x.copy()
            if logC:
                x = np.exp(x)
            if square:
                x *= x
            if tanh:
                # c = m + s*tanh((x - m)/s)
                # x = m + s*atanh((c - m)/s)
                x = sc_middle + sc_bounds_span*np.tanh(
                    (x0 - sc_middle)/sc_bounds_span)
            if reduced or len(const) > 0:
                idx_red = 0
                idx_elm = 0
                new_x = []
                for idx in range(len(init_concs)):
                    if idx in elim:
                        new_x.append(red_cbs[idx_elm](*x))
                        idx_elm += 1
                    elif idx in const_indices:
                        new_x.append(scaling*init_concs[idx])
                    else:
                        new_x.append(x[idx_red])
                        idx_red += 1
                x = np.array(new_x)
            # Sanity checks:
            neg_conc, too_much = np.any(x < 0), np.any(x > sc_upper_bounds*(1 + 1e-12))
            if neg_conc or too_much:
                print("neg_conc, too_much", neg_conc, too_much) #DEBUG
                return None, result
            x /= scaling
            return x, result
        else:
            return None, result

    def initial_guess(self, init_concs, scaling=1, repetitions=50,
                      steffensens=0):
        # Fixed point iteration
        def f(x):
            xnew = np.zeros_like(x)
            for ri, eq in enumerate(self.rxns):
                xnew += solve_equilibrium(
                    x, self.stoichs[:, ri],
                    eq.params * scaling**eq.dimensionality())
            return xnew/self.nr
        init_concs = init_concs*scaling
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
        return [scaling*ic + sum([coords[ri]*self.stoichs[cidx, ri]
                                  for ri in range(self.nr)])
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
            j = sp.Matrix(1, len(f), lambda _, q: f[q]).jacobian(r)
            j_lmbd = sp.lambdify(r, j, modules='numpy')

            def j_cb(arr):
                return j_lmbd(*arr)
        return f_cb, j_cb

    def initial_guess(self, init_concs, scaling, repetitions=10):
        def f(x):
            c0 = init_concs*scaling + np.dot(self.stoichs, x)
            xnew = np.zeros_like(x)
            for ri, eq in enumerate(self.rxns):
                xnew += _solve_equilibrium_coord(c0, self.stoichs[:, ri],
                                                 eq.params)
                # xnew = (xnew + x)/2
                # c0 = c0 + np.dot(self.stoichs, xnew)
                # xnew = _solve_equilibrium_coord(c0, self.stoichs[:, ri],
                #                                 eq.params)
                # xnew = (xnew + x)/2
                # c0 = c0 + np.dot(self.stoichs, xnew)
            return xnew/(2*self.nr+1)
        x0 = np.zeros(self.nr)
        for idx in range(repetitions):
            print(x0)
            x0 = (x0 + f(x0))/2  # (idx*x0 + f(x0))/(idx + 1)
        return x0

    def root(self, init_concs, scaling=1.0, init_iter=20,
             init_guess=None, x0=None, **kwargs):
        from scipy.optimize import root
        init_concs = self.as_per_substance_array(init_concs)
        f, j = self.num_cb_factory(init_concs, jac=True, scaling=scaling)
        if x0 is None:
            if init_guess is None:
                x0 = self.initial_guess(init_concs, scaling,
                                        repetitions=init_iter)
            else:
                x0 = init_guess
        print(x0)
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
