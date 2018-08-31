# -*- coding: utf-8 -*-
"""
Module collecting classes and functions for dealing with (multiphase) chemical
equilibria.

.. Note::

  This module is provisional at the moment, i.e. the API is not stable and may
  break without a deprecation cycle.

"""
from __future__ import division, absolute_import

import warnings

import numpy as np

from .chemistry import equilibrium_quotient, Equilibrium, Species
from .reactionsystem import ReactionSystem
from ._util import get_backend
from .util.pyutil import deprecated
from ._eqsys import EqCalcResult, NumSysLin, NumSysLog, NumSysSquare as _NumSysSquare


NumSysSquare = deprecated()(_NumSysSquare)


class EqSystem(ReactionSystem):

    _BaseReaction = Equilibrium
    _BaseSubstance = Species

    def html(self, *args, **kwargs):
        k = 'color_categories'
        kwargs[k] = kwargs.get(k, False)
        return super(EqSystem, self).html(*args, **kwargs)

    def eq_constants(self, non_precip_rids=(), eq_params=None, small=0):
        if eq_params is None:
            eq_params = [eq.param for eq in self.rxns]
        return [small if idx in non_precip_rids else
                eq for idx, eq in enumerate(eq_params)]

    def equilibrium_quotients(self, concs):
        stoichs = self.stoichs()
        return [equilibrium_quotient(concs, stoichs[ri, :])
                for ri in range(self.nr)]

    def stoichs_constants(self, eq_params=None, rref=False, Matrix=None,
                          backend=None, non_precip_rids=()):
        if eq_params is None:
            eq_params = self.eq_constants()
        if rref:
            from pyneqsys.symbolic import linear_rref
            be = get_backend(backend)
            rA, rb = linear_rref(self.stoichs(non_precip_rids),
                                 list(map(be.log, eq_params)),
                                 Matrix)
            return rA.tolist(), list(map(be.exp, rb))
        else:
            return (self.stoichs(non_precip_rids),
                    eq_params)

    def composition_conservation(self, concs, init_concs):
        composition_vecs, comp_keys = self.composition_balance_vectors()
        A = np.array(composition_vecs)
        return (comp_keys,
                np.dot(A, self.as_per_substance_array(concs).T),
                np.dot(A, self.as_per_substance_array(init_concs).T))

    def other_phase_species_idxs(self, phase_idx=0):
        return [idx for idx, s in enumerate(
            self.substances.values()) if s.phase_idx != phase_idx]

    @property
    @deprecated(last_supported_version='0.3.1', will_be_missing_in='0.8.0',
                use_instead=other_phase_species_idxs)
    def precipitate_substance_idxs(self):
        return [idx for idx, s in enumerate(
            self.substances.values()) if s.precipitate]

    def phase_transfer_reaction_idxs(self, phase_idx=0):
        return [idx for idx, rxn in enumerate(self.rxns)
                if rxn.has_precipitates(self.substances)]

    @property
    @deprecated(last_supported_version='0.3.1', will_be_missing_in='0.8.0',
                use_instead=phase_transfer_reaction_idxs)
    def precipitate_rxn_idxs(self):
        return [idx for idx, rxn in enumerate(self.rxns)
                if rxn.has_precipitates(self.substances)]

    def dissolved(self, concs):
        """ Return dissolved concentrations """
        new_concs = concs.copy()
        for r in self.rxns:
            if r.has_precipitates(self.substances):
                net_stoich = np.asarray(r.net_stoich(self.substances))
                s_net, s_stoich, s_idx = r.precipitate_stoich(self.substances)
                new_concs -= new_concs[s_idx]/s_stoich * net_stoich
        return new_concs

    def _fw_cond_factory(self, ri, rtol=1e-14):
        rxn = self.rxns[ri]

        def fw_cond(x, p):
            precip_stoich_coeff, precip_idx = rxn.precipitate_stoich(self.substances)[1:3]
            q = rxn.Q(self.substances, self.dissolved(x))
            k = rxn.equilibrium_constant()
            if precip_stoich_coeff > 0:
                return q*(1+rtol) < k
            elif precip_stoich_coeff < 0:
                return q > k*(1+rtol)
            else:
                raise NotImplementedError
        return fw_cond

    def _bw_cond_factory(self, ri, small):
        rxn = self.rxns[ri]

        def bw_cond(x, p):
            precipitate_idx = rxn.precipitate_stoich(self.substances)[2]
            if x[precipitate_idx] < small:
                return False
            else:
                return True
        return bw_cond

    def _SymbolicSys_from_NumSys(self, NS, conds, rref_equil, rref_preserv,
                                 new_eq_params=True):
        from pyneqsys.symbolic import SymbolicSys
        import sympy as sp
        ns = NS(self, backend=sp, rref_equil=rref_equil,
                rref_preserv=rref_preserv, precipitates=conds,
                new_eq_params=new_eq_params)
        symb_kw = {}
        if ns.pre_processor is not None:
            symb_kw['pre_processors'] = [ns.pre_processor]
        if ns.post_processor is not None:
            symb_kw['post_processors'] = [ns.post_processor]
        if ns.internal_x0_cb is not None:
            symb_kw['internal_x0_cb'] = ns.internal_x0_cb
        return SymbolicSys.from_callback(
            ns.f, self.ns, nparams=self.ns + (self.nr if new_eq_params else 0),
            **symb_kw)

    def get_neqsys_conditional_chained(self, rref_equil=False,
                                       rref_preserv=False, NumSys=NumSysLin, **kwargs):
        from pyneqsys import ConditionalNeqSys, ChainedNeqSys

        def factory(conds):
            return ChainedNeqSys([self._SymbolicSys_from_NumSys(
                NS, conds, rref_equil, rref_preserv, **kwargs
            ) for NS in NumSys])

        cond_cbs = [(self._fw_cond_factory(ri),
                     self._bw_cond_factory(ri, NumSys[0].small)) for
                    ri in self.phase_transfer_reaction_idxs()]
        return ConditionalNeqSys(cond_cbs, factory)

    def get_neqsys_chained_conditional(self, rref_equil=False,
                                       rref_preserv=False,
                                       NumSys=NumSysLin, **kwargs):
        from pyneqsys import ConditionalNeqSys, ChainedNeqSys

        def mk_factory(NS):
            def factory(conds):
                return self._SymbolicSys_from_NumSys(NS, conds, rref_equil,
                                                     rref_preserv, **kwargs)
            return factory

        return ChainedNeqSys(
            [ConditionalNeqSys(
                [(self._fw_cond_factory(ri),
                  self._bw_cond_factory(ri, NS.small)) for
                 ri in self.phase_transfer_reaction_idxs()],
                mk_factory(NS)
            ) for NS in NumSys])

    def get_neqsys_static_conditions(self, rref_equil=False,
                                     rref_preserv=False,
                                     NumSys=(NumSysLin,), precipitates=None, **kwargs):
        if precipitates is None:
            precipitates = (False,)*len(self.phase_transfer_reaction_idxs())
        from pyneqsys import ChainedNeqSys
        return ChainedNeqSys([self._SymbolicSys_from_NumSys(
            NS, precipitates, rref_equil, rref_preserv, **kwargs) for NS in NumSys])

    def get_neqsys(self, neqsys_type, NumSys=NumSysLin, **kwargs):
        new_kw = {'rref_equil': False, 'rref_preserv': False}
        if neqsys_type == 'static_conditions':
            new_kw['precipitates'] = None
        for k in new_kw:
            if k in kwargs:
                new_kw[k] = kwargs.pop(k)

        try:
            NumSys[0]
        except TypeError:
            new_kw['NumSys'] = (NumSys,)
        else:
            new_kw['NumSys'] = NumSys

        return getattr(self, 'get_neqsys_' + neqsys_type)(**new_kw)

    def non_precip_rids(self, precipitates):
        return [idx for idx, precip in zip(
            self.phase_transfer_reaction_idxs(), precipitates) if not precip]

    def _result_is_sane(self, init_concs, x, rtol=1e-9):
        sc_upper_bounds = np.array(self.upper_conc_bounds(init_concs))
        neg_conc, too_much = np.any(x < 0), np.any(
            x > sc_upper_bounds*(1 + rtol))
        if neg_conc or too_much:
            if neg_conc:
                warnings.warn("Negative concentration")
            if too_much:
                warnings.warn("Too much of at least one component")
            return False
        return True

    def _solve(self, init_concs, x0=None, NumSys=(NumSysLog, NumSysLin),
               neqsys='chained_conditional', **kwargs):
        if isinstance(neqsys, str):
            neqsys = self.get_neqsys(
                neqsys, NumSys=NumSys,
                rref_equil=kwargs.pop('rref_equil', False),
                rref_preserv=kwargs.pop('rref_preserv', False),
                precipitates=kwargs.pop('precipitates', None))
        if x0 is None:
            x0 = init_concs
        params = np.concatenate((init_concs, [float(elem) for elem
                                              in self.eq_constants()]))
        x, sol = neqsys.solve(x0, params, **kwargs)
        if not sol['success']:
            warnings.warn("Root-finding indicated as failed by solver.")
        sane = self._result_is_sane(init_concs, x)
        return x, sol, sane

    def solve(self, init_concs, varied=None, **kwargs):
        results = EqCalcResult(self, init_concs, varied)
        results.solve()
        return results

    def root(self, init_concs, x0=None, neqsys=None, NumSys=NumSysLog,
             neqsys_type='chained_conditional', **kwargs):
        init_concs = self.as_per_substance_array(init_concs)
        params = np.concatenate((init_concs, [float(elem) for elem
                                              in self.eq_constants()]))
        if neqsys is None:
            neqsys = self.get_neqsys(
                neqsys_type, NumSys=NumSys,
                rref_equil=kwargs.pop('rref_equil', False),
                rref_preserv=kwargs.pop('rref_preserv', False),
                precipitates=kwargs.pop('precipitates', None))
        if x0 is None:
            x0 = init_concs
        x, sol = neqsys.solve(x0, params, **kwargs)
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
            result = ['$' + s.latex_name + '$'
                      for s in self.substances.values()]
            return result
        else:
            return [s.name for s in self.substances.values()]

    def roots(self, init_concs, varied_data, varied, x0=None,
              NumSys=NumSysLog, plot_kwargs=None,
              neqsys_type='chained_conditional', **kwargs):
        """
        Parameters
        ----------
        init_concs : array or dict
        varied_data : array
        varied_idx : int or str
        x0 : array
        NumSys : _NumSys subclass
            See :class:`NumSysLin`, :class:`NumSysLog`, etc.
        plot_kwargs : dict
            See py:meth:`pyneqsys.NeqSys.solve`. Two additional keys
            are intercepted here:
                latex_names: bool (default: False)
                conc_unit_str: str (default: 'M')
        neqsys_type : str
            what method to use for NeqSys construction (get_neqsys_*)
        \\*\\*kwargs :
            Keyword argumetns passed on to py:meth:`pyneqsys.NeqSys.solve_series`.

        """
        _plot = plot_kwargs is not None
        if _plot:
            latex_names = plot_kwargs.pop('latex_names', False)
            conc_unit_str = plot_kwargs.pop('conc_unit_str', 'M')
            if 'ax' not in plot_kwargs:
                plot_kwargs['ax'] = self._get_default_plot_ax()

        init_concs = self.as_per_substance_array(init_concs)
        neqsys = self.get_neqsys(
            neqsys_type, NumSys=NumSys,
            rref_equil=kwargs.pop('rref_equil', False),
            rref_preserv=kwargs.pop('rref_preserv', False),
            precipitates=kwargs.pop('precipitates', None))
        if x0 is None:
            x0 = init_concs

        if _plot:
            cb = neqsys.solve_and_plot_series
            if 'plot_kwargs' not in kwargs:
                kwargs['plot_kwargs'] = plot_kwargs
            if 'labels' not in kwargs['plot_kwargs']:
                kwargs['plot_kwargs']['labels'] = (
                    self.substance_labels(latex_names))
            if 'substances' in plot_kwargs:
                if 'indices' in plot_kwargs:
                    raise ValueError("Now I am confused..")
                kwargs['plot_kwargs']['indices'] = map(
                    self.as_substance_index, plot_kwargs.pop('substances'))
                print(kwargs['plot_kwargs']['indices'])
        else:
            cb = neqsys.solve_series

        params = np.concatenate((init_concs, self.eq_constants()))
        xvecs, info_dicts = cb(
            x0, params, varied_data, self.as_substance_index(varied),
            propagate=False, **kwargs)
        sanity = [self._result_is_sane(init_concs, x) for x in xvecs]

        if _plot:
            import matplotlib.pyplot as plt
            from pyneqsys.plotting import mpl_outside_legend
            mpl_outside_legend(plt.gca())
            varied_subst = self.substances[varied]
            xlbl = ('$[' + varied_subst.latex_name + ']_0$' if latex_names
                    else str(varied_subst))
            plt.gca().set_xlabel(xlbl + ' / ' + conc_unit_str)
            plt.gca().set_ylabel('Concentration / ' + conc_unit_str)

        return xvecs, info_dicts, sanity

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
            # TODO: handle precipitate phases in plotting Q error
            qs = self.equilibrium_quotients(concs)
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
