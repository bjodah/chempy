# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from itertools import product
import math

try:
    import numpy as np
except ImportError:
    np = None

from .printing import number_to_scientific_html
from ._util import get_backend, mat_dot_vec, prodpow


class EqCalcResult(object):

    attrs = {
        'sane': bool, 'success': bool,
        'nfev': int, 'njev': int,
        'time_cpu': float, 'time_wall': float
    }

    def __init__(self, eqsys, init_concs, varied):
        self.eqsys = eqsys
        self.all_inits, self.varied_keys = self.eqsys.per_substance_varied(init_concs, varied)
        self.conc = np.empty_like(self.all_inits)
        for k, v in self.attrs.items():
            setattr(self, k, np.zeros(self.all_inits.shape[:-1], dtype=v))

    def solve(self, **kwargs):
        for index in product(*map(range, self.all_inits.shape[:-1])):
            slc = tuple(index) + (slice(None),)
            self.conc[slc], nfo, sane = self.eqsys._solve(self.all_inits[slc], **kwargs)
            self.sane[index] = sane

            def _get(k):
                try:
                    return nfo[k]
                except TypeError:
                    return nfo[-1][k]

            for k in self.attrs:
                if k == 'sane':
                    continue
                try:
                    getattr(self, k)[index] = _get(k)
                except KeyError:
                    pass

    def _repr_html_(self):
        def fmt(num):
            return number_to_scientific_html(num, fmt=5)
        if len(self.varied_keys) == 0:
            raise NotImplementedError()
        elif len(self.varied_keys) == 1:
            var_html = self.eqsys.substances[self.varied_keys[0]].html_name
            header = ["[%s]<sub>0</sub>" % var_html] + ["[%s]" % s.html_name for s in self.eqsys.substances.values()]

            def row(i):
                j = self.eqsys.as_substance_index(self.varied_keys[0])
                return map(fmt, [self.all_inits[i, j]] + self.conc[i, :].tolist())
            pre = "  <td style='font-weight: bold;'>\n      "
            linker = "\n    </td>\n    <td>\n      "
            post = "\n    </td>"
            rows = [pre + linker.join(row(i)) + post for i in range(self.all_inits.shape[0])]
            template = """<table>\n  <tr>\n    <th>\n    %s\n    </th>\n  </tr>\n  <tr>\n  %s\n  </tr>\n</table>"""
            head_linker = "\n    </th>\n    <th>\n      "
            row_linker = "\n  </tr>\n  <tr>\n  "
            return template % (head_linker.join(header), row_linker.join(rows))
        else:
            raise NotImplementedError()

    def plot(self, ls=('-', '--', ':', '-.'), c=('k', 'r', 'g', 'b', 'c', 'm', 'y'), latex=None):
        import matplotlib.pyplot as plt
        if latex is None:
            latex = next(iter(self.eqsys.substances.values())).latex_name is not None
        if len(self.varied_keys) == 0:
            raise NotImplementedError()
        elif len(self.varied_keys) == 1:
            x = self.all_inits[:, self.eqsys.as_substance_index(self.varied_keys[0])]
            for idx, (k, v) in enumerate(self.eqsys.substances.items()):
                lbl = (r'$\mathrm{' + v.latex_name + '}$') if latex else v.name
                plt.plot(x, self.conc[:, idx], label=lbl, ls=ls[idx % len(ls)], c=c[idx % len(c)])

            ax = plt.gca()

            # Log-log
            ax.set_xscale('log')
            ax.set_yscale('log')

            # Axis labels
            var_latex = self.eqsys.substances[self.varied_keys[0]].latex_name
            ax.set_xlabel((r"$[\mathrm{%s}]_0$" if latex else "[%s]0") % var_latex)
            ax.set_ylabel("Concentration")

            # Outside legend
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
            # Put a legend to the right of the current axis
            ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        else:
            raise NotImplementedError()


class _NumSys(object):

    small = 0  # precipitation limit
    pre_processor = None
    post_processor = None
    internal_x0_cb = None

    def __init__(self, eqsys, rref_equil=False, rref_preserv=False,
                 backend=None, precipitates=(), new_eq_params=True):
        self.eqsys = eqsys
        self.rref_equil = rref_equil
        self.rref_preserv = rref_preserv
        self.backend = get_backend(backend)
        self.precipitates = precipitates
        self.new_eq_params = new_eq_params

    def _get_A_ks(self, eq_params):
        non_precip_rids = self.eqsys.non_precip_rids(self.precipitates)
        return self.eqsys.stoichs_constants(
            self.eqsys.eq_constants(non_precip_rids, eq_params, self.small),
            self.rref_equil, backend=self.backend, non_precip_rids=non_precip_rids)

    def _inits_and_eq_params(self, params):
        eq_params = params[self.eqsys.ns:]
        if not self.new_eq_params:
            assert not eq_params, "Adjust number of parameters accordingly"
            eq_params = None  # use those of eqsys
        return params[:self.eqsys.ns], eq_params


class NumSysLin(_NumSys):

    def internal_x0_cb(self, init_concs, params):
        # reduce risk of stationary starting point
        return (99*init_concs + self.eqsys.dissolved(init_concs))/100

    def f(self, yvec, params):
        from pyneqsys.symbolic import linear_exprs
        init_concs, eq_params = self._inits_and_eq_params(params)
        A, ks = self._get_A_ks(eq_params)
        # yvec == C
        f_equil = [q/k - 1 if k != 0 else q for q, k
                   in zip(prodpow(yvec, A), ks)]
        B, comp_nrs = self.eqsys.composition_balance_vectors()
        f_preserv = linear_exprs(B, yvec, mat_dot_vec(B, init_concs),
                                 rref=self.rref_preserv)
        return f_equil + f_preserv


class _NumSysLinNegPenalty(NumSysLin):

    def f(self, yvec, params):
        import sympy as sp
        f_penalty = [sp.Piecewise((yi**2, yi < 0), (0, True)) for yi in yvec]
        return super(_NumSysLinNegPenalty, self).f(yvec, params) + f_penalty


class NumSysLinRel(NumSysLin):

    def max_concs(self, params, min_=min, dtype=np.float64):
        init_concs = params[:self.eqsys.ns]
        return self.eqsys.upper_conc_bounds(init_concs, min_=min_, dtype=dtype)

    def pre_processor(self, x, params):
        return x/self.max_concs(params), params

    def post_processor(self, x, params):
        return x*self.max_concs(params), params

    def f(self, yvec, params):
        import sympy as sp
        return NumSysLin.f(self, [m*yi for m, yi in zip(
            self.max_concs(params, min_=lambda x: sp.Min(*x), dtype=object),
            yvec)], params)


class NumSysSquare(NumSysLin):

    small = 1e-35

    def pre_processor(self, x, params):
        return (np.sqrt(np.abs(x)), params)

    def post_processor(self, x, params):
        return x**2, params

    def internal_x0_cb(self, init_concs, params):
        return np.sqrt(np.abs(init_concs))

    def f(self, yvec, params):
        ysq = [yi*yi for yi in yvec]
        return NumSysLin.f(self, ysq, params)


class NumSysLinTanh(NumSysLin):

    def pre_processor(self, x, params):
        ymax = self.eqsys.upper_conc_bounds(params[:self.eqsys.ns])
        return np.arctanh((8*x/ymax - 4) / 5), params

    def post_processor(self, x, params):
        ymax = self.eqsys.upper_conc_bounds(params[:self.eqsys.ns])
        return ymax*(4 + 5*np.tanh(x))/8, params

    def internal_x0_cb(self, init_concs, params):
        return self.pre_processor(init_concs, init_concs)[0]

    def f(self, yvec, params):
        import sympy
        ymax = self.eqsys.upper_conc_bounds(
            params[:self.eqsys.ns],
            min_=lambda a, b: sympy.Piecewise((a, a < b), (b, True)))
        ytanh = [yimax*(4 + 5*sympy.tanh(yi))/8
                 for yimax, yi in zip(ymax, yvec)]
        return NumSysLin.f(self, ytanh, params)


class NumSysLog(_NumSys):

    small = math.exp(-36)  # anything less than `small` is insignificant

    def pre_processor(self, x, params):
        return (np.log(np.asarray(x) + NumSysLog.small),  # 10: damping
                params)  # zero conc. ~= small

    def post_processor(self, x, params):
        return np.exp(x), params

    def internal_x0_cb(self, init_concs, params):
        # return [1]*len(init_concs)
        return [0.1]*len(init_concs)

    def f(self, yvec, params):
        from pyneqsys.symbolic import linear_exprs
        init_concs, eq_params = self._inits_and_eq_params(params)
        A, ks = self._get_A_ks(eq_params)
        # yvec == ln(C)
        f_equil = mat_dot_vec(A, yvec, [-self.backend.log(k) for k in ks])
        B, comp_nrs = self.eqsys.composition_balance_vectors()
        f_preserv = linear_exprs(B, list(map(self.backend.exp, yvec)),
                                 mat_dot_vec(B, init_concs),
                                 rref=self.rref_preserv)
        return f_equil + f_preserv
