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
    PartiallySolvedSystem = None
else:
    from pyodesys.symbolic import PartiallySolvedSystem


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
      % for si, subst_key in enumerate(getattr(odesys, 'free_names', odesys.names)):
       % if len(subst_comp[si]) > 0:
        bounds[${si}] = vecmin(${', '.join([('cc[%d]/%d' % (ci, n)) if n != 1 else 'cc[%d]' % ci for ci, n in subst_comp[si].items()])});
       % else:
        bounds[${si}] = INFINITY;
       % endif
      % endfor
        return bounds;
    }
"""  # noqa


_first_step = """
    m_upper_bounds = upper_conc_bounds(${init_conc});
    m_lower_bounds.resize(${odesys.ny});
    return m_rtol*std::min(get_dx_max(x, y), 1.0);
"""  # if (m_upper_bounds.size() == 0)


def _get_comp_conc(rsys, odesys, comp_keys, skip_keys):
    comp_conc = []
    for comp_key in comp_keys:
        if comp_key in skip_keys:
            continue  # see Substance.__doc__
        _d = OrderedDict()
        for si, subst_key in enumerate(odesys.names):
            coeff = rsys.substances[subst_key].composition.get(comp_key, 0)
            if coeff != 0:
                _d[si] = coeff
        comp_conc.append(_d)
    return comp_conc


def _get_subst_comp(rsys, odesys, comp_keys, skip_keys):
    subst_comp = []
    for subst_key in odesys.names:
        _d = OrderedDict()
        for k, v in rsys.substances[subst_key].composition.items():
            if k in skip_keys:
                continue
            _d[comp_keys.index(k)] = v
        subst_comp.append(_d)
    return subst_comp


def _render(tmpl, **kwargs):
    from mako.template import Template
    from mako.exceptions import text_error_template
    try:
        return str(Template(tmpl).render(**kwargs))
    except:
        sys.stderr.write(text_error_template().render())
        raise


def get_native(rsys, odesys, integrator, skip_keys=(0,)):
    comp_keys = Substance.composition_keys(rsys.substances.values(), skip_keys=skip_keys)
    if isinstance(odesys, PartiallySolvedSystem):
        init_conc = '&m_p[%d]' % (len(odesys.params) - len(odesys.original_dep))
    else:
        init_conc = 'y'

    return native_sys[integrator].from_other(odesys, namespace_override={
        'p_anon': _render(_anon, odesys=odesys, ncomp=len(comp_keys),
                          comp_conc=_get_comp_conc(rsys, odesys, comp_keys, skip_keys),
                          subst_comp=_get_subst_comp(rsys, odesys, comp_keys, skip_keys)),
        'p_first_step': _render(_first_step, init_conc=init_conc, odesys=odesys),
        'p_get_dx_max': True
    }, namespace_extend={
        'p_includes': ["<type_traits>",  "<vector>"]
    })
