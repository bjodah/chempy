# -*- coding: utf-8 -*-
"""
Non-public API (expect changes without notice).

Helper functions for using native code generation together with pyodesys.
"""
from __future__ import print_function, absolute_import, division

from collections import OrderedDict
import sys

try:
    from pyodesys.native import native_sys
except ImportError:
    native_sys = None

from .. import Substance

_anon = """
    template <typename T>
    constexpr T vecmin(T&& a){
        return std::forward<T>(a);
    }
    template <typename T1, typename T2>
    constexpr typename std::common_type<T1, T2>::type vecmin(T1&& a, T2&& b){
        return (a < b) ? std::forward<T1>(a) : std::forward<T2>(b);
    }
    template <typename T, typename... Ts>
    constexpr typename std::common_type<T, Ts...>::type vecmin(T&& a, Ts&&... args){
        return vecmin(std::forward<T>(a), vecmin(std::forward<Ts>(args)...));
    }
    std::vector<double> upper_conc_bounds(const double * const y){
        auto bounds = std::vector<double>(${odesys.ny});
        double cc[${ncomp}];
      % for ci in range(ncomp):
        cc[${ci}] = ${' + '.join([('%d*y[%d]' % (v, k)) if v != 1 else 'y[%d]' % k for k, v in comp_conc[ci].items()])};
      % endfor
      % for si, subst_key in enumerate(odesys.names):
        bounds[${si}] = vecmin(${', '.join([('cc[%d]/%d' % (ci, n)) if n != 1 else 'cc[%d]' % ci for ci, n in subst_comp[si].items()])});
      % endfor
        return bounds;
    }
"""  # noqa


_first_step = """
    double max_euler_step;
    auto bounds = upper_conc_bounds(y);
    auto fvec = std::vector<double>(%(ny)d);
    auto hvec = std::vector<double>(%(ny)d);
    rhs(t, y, &fvec[0]);
    for (int idx=0; idx<%(ny)d; ++idx){
        if (fvec[idx] == 0) {
            hvec[idx] = std::numeric_limits<double>::infinity();
        } else if (fvec[idx] > 0) {
            hvec[idx] = (bounds[idx] - y[idx])/fvec[idx];
        } else { // fvec[idx] < 0
            hvec[idx] = -y[idx]/fvec[idx];
        }
    }
    max_euler_step = *std::min_element(std::begin(hvec), std::end(hvec));
    return m_rtol*std::min(max_euler_step, 1.0);
"""


def _get_comp_conc(rsys, odesys, comp_keys):
    comp_conc = []
    for comp_key in comp_keys:
        if comp_key == 0:
            continue  # see Substance.__doc__
        _d = OrderedDict()
        for si, subst_key in enumerate(odesys.names):
            coeff = rsys.substances[subst_key].composition.get(comp_key, 0)
            if coeff != 0:
                _d[si] = coeff
        comp_conc.append(_d)
    return comp_conc


def _get_subst_comp(rsys, odesys, comp_keys):
    subst_comp = []
    for subst_key in odesys.names:
        _d = OrderedDict()
        for k, v in rsys.substances[subst_key].composition.items():
            _d[comp_keys.index(k)] = v
        subst_comp.append(_d)
    return subst_comp


def get_native(rsys, odesys, integrator):
    from mako.template import Template
    from mako.exceptions import text_error_template
    comp_keys = Substance.composition_keys(rsys.substances.values())

    try:
        anon = Template(_anon).render(odesys=odesys, ncomp=len(comp_keys),
                                      comp_conc=_get_comp_conc(rsys, odesys, comp_keys),
                                      subst_comp=_get_subst_comp(rsys, odesys, comp_keys))
    except:
        sys.stderr.write(text_error_template().render())
        raise

    return native_sys[integrator].from_other(odesys, namespace_override={
        'p_anon': anon,
        'p_first_step': _first_step % {'ny': odesys.ny}
    }, namespace_extend={
        'p_includes': ["<algorithm>", "<limits>", "<type_traits>",  "<vector>"]
    })
