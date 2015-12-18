from __future__ import division, absolute_import

import math
import os
import warnings
from collections import defaultdict
from functools import reduce, partial
from operator import mul, add


import numpy as np

from .chemistry import ReactionSystem, equilibrium_quotient


def reducemap(args, reduce_op, map_op):
    return reduce(reduce_op, map(map_op, *args))


def vec_dot_vec(vec1, vec2):
    # return np.dot(vec1, vec2)
    # return np.add.reduce(np.multiply(vec1, vec2))
    return reducemap((vec1, vec2), add, mul)


def prodpow(bases, exponents):
    exponents = np.asarray(exponents)
    return np.multiply.reduce(bases**exponents, axis=-1)


def mat_dot_vec(iter_mat, iter_vec, iter_term=None):  # pure python (slow)
    if iter_term is None:
        return [vec_dot_vec(row, iter_vec) for row in iter_mat]
    else:
        # daxpy
        return [vec_dot_vec(row, iter_vec) + term for row, term
                in zip(iter_mat, iter_term)]


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
    if not hasattr(stoich, 'ndim') or stoich.ndim == 1:
        c = c0 + stoich*rc
    else:
        c = c0 + np.dot(stoich, rc)
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


def composition_balance(substances, concs, composition_number):
    if not hasattr(concs, 'ndim') or concs.ndim == 1:
        res = 0
    elif concs.ndim == 2:
        res = np.zeros(concs.shape[0])
        concs = concs.T
    else:
        raise NotImplementedError
    for s, c in zip(substances, concs):
        res += s.composition.get(composition_number, 0)*c
    return res


class _NumSys(object):

    small = 0  # precipitation limit
    pre_processor = None
    post_processor = None

    def __init__(self, eqsys, rref_equil=False, rref_preserv=False, ln=None,
                 exp=None, precipitates=()):
        if ln is None:
            ln = math.log
        if exp is None:
            exp = math.exp
        self.eqsys = eqsys
        self.rref_equil = rref_equil
        self.rref_preserv = rref_preserv
        self.ln = ln
        self.exp = exp
        self.precipitates = precipitates


class NumSysLin(_NumSys):

    def internal_x0(self, init_concs):
        return init_concs

    def f(self, yvec, params):
        from pyneqsys.symbolic import linear_exprs
        init_concs, eq_params = params[:self.eqsys.ns], params[self.eqsys.ns:]
        non_precip_rids = self.eqsys.non_precip_rids(self.precipitates)
        A, ks = self.eqsys.stoichs_constants(
            self.eqsys.eq_constants(non_precip_rids, eq_params, self.small),
            self.rref_equil, ln=self.ln, exp=self.exp,
            non_precip_rids=non_precip_rids)
        # yvec == C
        f1 = [q/k-1 if k != 0 else q for q, k in zip(prodpow(yvec, A), ks)]
        B, comp_nrs = self.eqsys.composition_balance_vectors()
        f2 = linear_exprs(B, yvec, mat_dot_vec(B, init_concs))
        # import sympy as sp
        # f3 = [sp.Piecewise((yi**2, yi < 0), (0, True)) for yi in yvec]
        # f3 = [sp.ITE(yi < 0, yi**2, 0) for yi in yvec]
        return f1 + f2  # + f3


class NumSysSquare(NumSysLin):

    small = 1e-30

    @staticmethod
    def pre_processor(x, params):
        return (np.sqrt(np.abs(x)), params)

    @staticmethod
    def post_processor(x, params):
        return x**2, params

    def internal_x0(self, init_concs):
        return np.sqrt(np.abs(init_concs))

    def f(self, yvec, params):
        ysq = [yi*yi for yi in yvec]
        return NumSysLin.f(self, ysq, params)


class NumSysLog(_NumSys):

    small = math.exp(-60)  # anything less than `small` is insignificant

    @staticmethod
    def pre_processor(x, params):
        return (np.log(np.asarray(x) + NumSysLog.small)/10,  # 10: damping
                params)  # zero conc. ~= small

    @staticmethod
    def post_processor(x, params):
        return np.exp(x), params

    def internal_x0(self, init_concs):
        #return [1]*len(init_concs)
        return [0.1]*len(init_concs)
        # np.log(np.abs(init_concs))/10 # [0]*len(init_concs)

    def f(self, yvec, params):
        from pyneqsys.symbolic import linear_exprs
        init_concs, eq_params = params[:self.eqsys.ns], params[self.eqsys.ns:]
        non_precip_rids = self.eqsys.non_precip_rids(self.precipitates)
        A, ks = self.eqsys.stoichs_constants(
            self.eqsys.eq_constants(non_precip_rids, eq_params, self.small),
            self.rref_equil, ln=self.ln, exp=self.exp,
            non_precip_rids=non_precip_rids)
        f1 = mat_dot_vec(A, yvec, [-self.ln(k) for k in ks])  # yvec == ln(C)
        B, comp_nrs = self.eqsys.composition_balance_vectors()
        f2 = linear_exprs(B, list(map(self.exp, yvec)),
                          mat_dot_vec(B, init_concs),
                          rref=self.rref_preserv)
        return f1 + f2


class EqSystem(ReactionSystem):

    def eq_constants(self, non_precip_rids=(), eq_params=None, small=0):
        if eq_params is None:
            eq_params = [eq.param for eq in self.rxns]
        return np.array([small if idx in non_precip_rids else
                         eq_params[idx] for idx, eq in enumerate(eq_params)])

    def upper_conc_bounds(self, init_concs):
        init_concs_arr = self.as_per_substance_array(init_concs)
        composition_conc = defaultdict(float)
        for conc, s_obj in zip(init_concs_arr, self.substances.values()):
            for comp_nr, coeff in s_obj.composition.items():
                if comp_nr == 0:
                    continue
                composition_conc[comp_nr] += coeff*conc
        bounds = []
        for s_obj in self.substances.values():
            upper = float('inf')
            for comp_nr, coeff in s_obj.composition.items():
                if comp_nr == 0:
                    continue
                upper = min(upper, composition_conc[comp_nr]/coeff)
            bounds.append(upper)
        return bounds

    def equilibrium_quotients(self, concs):
        stoichs = self.stoichs()
        return [equilibrium_quotient(concs, stoichs[ri, :])
                for ri in range(self.nr)]

    def stoichs_constants(self, eq_params, rref=False, Matrix=None,
                          ln=None, exp=None, non_precip_rids=()):
        if rref:
            from pyneqsys.symbolic import linear_rref
            ln = ln or math.log
            rA, rb = linear_rref(self.stoichs(non_precip_rids),
                                 map(ln, eq_params),
                                 Matrix)
            exp = exp or math.exp
            return rA.tolist(), list(map(exp, rb))
        else:
            return (self.stoichs(non_precip_rids),
                    eq_params)

    def composition_balance_vectors(self):
        composition_keys = set()
        for s in self.substances.values():
            for key in s.composition:
                composition_keys.add(key)
        vs = []
        sorted_composition_keys = sorted(composition_keys)
        for key in sorted_composition_keys:
            vs.append([s.composition.get(key, 0) for s in self.substances.values()])
        return vs, sorted_composition_keys

    def composition_conservation(self, concs, init_concs):
        composition_vecs, comp_keys = self.composition_balance_vectors()
        A = np.array(composition_vecs)
        return (comp_keys,
                np.dot(A, self.as_per_substance_array(concs).T),
                np.dot(A, self.as_per_substance_array(init_concs).T))

    @property
    def solid_substance_idxs(self):
        return [idx for idx, s in enumerate(self.substances.values()) if s.solid]

    @property
    def solid_rxn_idxs(self):
        return [idx for idx, rxn in enumerate(self.rxns)
                if rxn.has_solids(self.substances)]

    def fw_cond_factory(self, ri):
        """ """
        rxn = self.rxns[ri]

        def fw_cond(x, p):
            solid_stoich_coeff = rxn.solid_stoich(self.substances)[1]
            q = rxn.Q(self.substances, x)
            k = rxn.K()
            QoverK = q/k
            if solid_stoich_coeff > 0:
                return QoverK < 1
            elif solid_stoich_coeff < 0:
                return QoverK > 1
            else:
                raise NotImplementedError
        return fw_cond

    def bw_cond_factory(self, ri, small):
        rxn = self.rxns[ri]

        def bw_cond(x, p):
            solid_idx = rxn.solid_stoich(self.substances)[2]
            if x[solid_idx] < small:
                return False
            else:
                return True
        return bw_cond

    def get_neqsys_x0(self, init_concs, rref_equil=False, rref_preserv=False,
                      NumSys=(NumSysLin,), **kwargs):

        if len(NumSys) > 1:
            from pyneqsys import ChainedNeqSys
            neqsys_x0_pairs = [
                self.get_neqsys_x0(init_concs, rref_equil,
                                   rref_preserv, (_NS,), **kwargs)
                for _NS in NumSys
            ]
            return (
                ChainedNeqSys(zip(*neqsys_x0_pairs)[0]),
                neqsys_x0_pairs[0][1]
            )

        import sympy as sp

        def _get_numsys_kwargs(precipitates):
            numsys = NumSys[0](
                self, ln=sp.log, exp=sp.exp, rref_equil=rref_equil,
                rref_preserv=rref_preserv, precipitates=precipitates)
            new_kwargs = kwargs.copy()
            if numsys.pre_processor is not None:
                new_kwargs['pre_processors'] = [numsys.pre_processor]
            if numsys.post_processor is not None:
                new_kwargs['post_processors'] = [numsys.post_processor]
            return numsys, new_kwargs

        from pyneqsys import SymbolicSys

        if len(self.solid_rxn_idxs) == 0:
            numsys, new_kw = _get_numsys_kwargs(())
            return (
                SymbolicSys.from_callback(
                    numsys.f, self.ns, nparams=self.ns+self.nr, **new_kw),
                numsys.internal_x0(init_concs)
            )
        else:
            # we have multiple equation systems corresponding
            # to the permutations of presence of each solid phase
            from pyneqsys import ConditionalNeqSys

            def factory(conds):
                numsys, new_kw = _get_numsys_kwargs(conds)
                return SymbolicSys.from_callback(
                    numsys.f, self.ns, nparams=self.ns+self.nr, **new_kw)
            cond_cbs = [(self.fw_cond_factory(ri),
                         self.bw_cond_factory(ri, NumSys[0].small)) for
                        ri in self.solid_rxn_idxs]
            return (ConditionalNeqSys(cond_cbs, factory),
                    _get_numsys_kwargs(())[0].internal_x0(init_concs))

    def get_neqsys(self, init_concs, rref_equil=False, rref_preserv=False,
                   NumSys=(NumSysLin,), **kwargs):
        from pyneqsys import ConditionalNeqSys, ChainedNeqSys

        def factory(conds):
            numsys, new_kw = _get_numsys_kwargs(conds)
            return ChainedNeqSys.from_callback(NumSys,
                numsys.f, self.ns, nparams=self.ns+self.nr, **new_kw)
        cond_cbs = [(self.fw_cond_factory(ri),
                     self.bw_cond_factory(ri, NumSys[0].small)) for
                    ri in self.solid_rxn_idxs]
        return ConditionalNeqSys(cond_cbs, factory)


    def non_precip_rids(self, precipitates):
        return [idx for idx, precip in zip(
            self.solid_rxn_idxs, precipitates) if not precip]

    def _result_is_sane(self, init_concs, x):
        sc_upper_bounds = np.array(self.upper_conc_bounds(init_concs))
        neg_conc, too_much = np.any(x < 0), np.any(
            x > sc_upper_bounds*(1 + 1e-12))
        if neg_conc or too_much:
            if neg_conc:
                warnings.warn("Negative concentration")
            if too_much:
                warnings.warn("Too much of at least one component")
            return False
        return True

    def root(self, init_concs, x0=None, neqsys=None,
             solver=None, NumSys=(NumSysLin,), **kwargs):
        init_concs = self.as_per_substance_array(init_concs)
        params = np.concatenate((init_concs, self.eq_constants()))
        internal_x0 = None
        if neqsys is None:
            neqsys, internal_x0 = self.get_neqsys_x0(
                init_concs,
                rref_equil=kwargs.pop('rref_equil', False),
                rref_preserv=kwargs.pop('rref_preserv', False),
                NumSys=NumSys
            )
            if x0 is None:
                x0, _ = neqsys.post_process(internal_x0, params)
        x, sol = neqsys.solve(solver, x0, params, internal_x0, **kwargs)
        if not sol['success']:
            warnings.warn("Root finding indicated as failed by solver.")
        sane = self._result_is_sane(init_concs, x)
        return x, sol, sane

    @staticmethod
    def _get_default_plot_ax(subplot_kwargs=None):
        import matplotlib.pyplot as plt
        if subplot_kwargs is None:
            subplot_kwargs = dict(xscale='log', yscale='log')
        return plt.subplot(1, 1, 1, **subplot_kwargs)

    def substance_labels(self, latex=False):
        if latex:
            result = ['$' + s.latex_name + '$' for s in self.substances.values()]
            return result
        else:
            return [s.name for s in self.substances.values()]

    def plot(self, xres, varied_data, conc_unit_str='M', subplot_kwargs=None,
             latex_names=False, ax=None, **kwargs):
        """ plots results from roots() """
        if ax is None:
            ax = self._get_default_plot_ax(subplot_kwargs)
        from pyneqsys.plotting import plot_series, mpl_outside_legend
        plot_series(xres, varied_data, labels=self.substance_labels(
            latex_names), ax=ax, **kwargs)
        mpl_outside_legend(ax)
        xlbl = '$[' + varied.latex_name + ']$' if latex_names else str(varied)
        ax.set_xlabel(xlbl + ' / %s' % conc_unit_str)
        ax.set_ylabel('Concentration / %s' % conc_unit_str)

    def roots(self, init_concs, varied, varied_data, x0=None, solver=None,
              NumSys=(NumSysLin,), plot_kwargs=None, **kwargs):
        plot = plot_kwargs is not None
        if plot:
            if plot_kwargs is True:
                plot_kwargs = {}
            latex_names = plot_kwargs.pop('latex_names', False)
            conc_unit_str = plot_kwargs.pop('conc_unit_str', 'M')
            if 'ax' not in plot_kwargs:
                plot_kwargs['ax'] = self._get_default_plot_ax()

        new_kwargs = kwargs.copy()
        init_concs = self.as_per_substance_array(init_concs)
        neqsys, x0_ = self.get_neqsys_x0(
            init_concs,
            rref_equil=new_kwargs.pop('rref_equil', False),
            rref_preserv=new_kwargs.pop('rref_preserv', False),
            NumSys=NumSys
        )
        if x0 is None:
            x0 = x0_

        if plot:
            cb = neqsys.solve_and_plot_series
            if 'plot_series_ax' not in new_kwargs:
                new_kwargs['plot_series_ax'] = plot_kwargs.pop('ax', None)
            if 'plot_series_kwargs' not in new_kwargs:
                new_kwargs['plot_series_kwargs'] = {}
            if 'labels' not in new_kwargs['plot_series_kwargs']:
                new_kwargs['plot_series_kwargs']['labels'] = (
                    self.substance_labels(latex_names))
            if len(plot_kwargs) > 0:
                raise KeyError("Unhandled kwarg keys: %s" % str(
                    plot_kwargs.keys()))
        else:
            cb = neqsys.solve_series

        params = np.concatenate((init_concs, self.eq_constants()))
        xvecs, sols = cb(solver,
                         x0, params, varied_data,
                         self.as_substance_index(varied), **new_kwargs)
        sanity = [self._result_is_sane(init_concs, x) for x in xvecs]

        if plot:
            import matplotlib.pyplot as plt
            from pyneqsys.plotting import mpl_outside_legend
            mpl_outside_legend(plt.gca())
            varied_subst = self.substances[varied]
            xlbl = ('$[' + varied_subst.latex_name + ']$' if latex_names
                    else str(varied_subst))
            plt.gca().set_xlabel(xlbl + ' / %s' % conc_unit_str)
            plt.gca().set_ylabel('Concentration / %s' % conc_unit_str)

        return xvecs, sols, sanity

    def plot_errors(self, concs, init_concs, varied_data, varied, axes=None,
                    compositions=True, Q=True, subplot_kwargs=None):
        if axes is None:
            import matplotlib.pyplot as plt
            if subplot_kwargs is None:
                subplot_kwargs = dict(xscale='log')
            fig, axes = plt.subplots(1, 2, figsize=(10, 4),
                                     subplot_kw=subplot_kwargs)
        varied_idx = self.as_substance_index(varied)
        ls, c = '- -- : -.'.split(), 'krgbcmy'
        all_inits = np.tile(self.as_per_substance_array(init_concs),
                            (len(varied_data), 1))
        all_inits[:, varied_idx] = varied_data
        if compositions:
            cmp_nrs, m1, m2 = self.composition_conservation(concs, all_inits)
            for cidx, (cmp_nr, a1, a2) in enumerate(zip(cmp_nrs, m1, m2)):
                axes[0].plot(concs[:, varied_idx],
                             a1-a2, label='Comp ' + str(cmp_nr),
                             ls=ls[cidx % len(ls)], c=c[cidx % len(c)])
                axes[1].plot(concs[:, varied_idx],
                             (a1-a2)/np.abs(a2), label='Comp ' + str(cmp_nr),
                             ls=ls[cidx % len(ls)], c=c[cidx % len(c)])

        if Q:
            qs = self.equilibrium_quotients(concs)  # TODO: handle solid phases
            ks = [rxn.param for rxn in self.rxns]
            for idx, (q, k) in enumerate(zip(qs, ks)):
                axes[0].plot(concs[:, varied_idx],
                             q-k, label='K R:' + str(idx),
                             ls=ls[(idx+cidx) % len(ls)],
                             c=c[(idx+cidx) % len(c)])
                axes[1].plot(concs[:, varied_idx],
                             (q-k)/k, label='K R:' + str(idx),
                             ls=ls[(idx+cidx) % len(ls)],
                             c=c[(idx+cidx) % len(c)])

        from pyneqsys.plotting import mpl_outside_legend
        mpl_outside_legend(axes[0])
        mpl_outside_legend(axes[1])
        axes[0].set_title("Absolute errors")
        axes[1].set_title("Relative errors")
