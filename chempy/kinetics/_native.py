# -*- coding: utf-8 -*-
"""
Non-public API (expect changes without notice).

Helper functions for using native code generation together with pyodesys.
"""
from __future__ import print_function, absolute_import, division

from collections import OrderedDict

try:
    from pyodesys.native import native_sys
except ImportError:
    native_sys = None
    PartiallySolvedSystem = None
    render_mako = None
else:
    from pyodesys.symbolic import PartiallySolvedSystem
    from pyodesys.native.util import render_mako


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
        bounds[${si}] = vecmin(${', '.join(['INFINITY' if n == 0 else ('cc[%d]/%d' % (ci, n)) if n != 1 else 'cc[%d]' % ci for ci, n in subst_comp[si].items()])});
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
"""

_roots_ss = """
    const int ny = get_ny();
    std::vector<double> f(ny);
    double tot=0.0;
    rhs(x, y, &f[0]);
    for (int i=0; i<ny; ++i){
        tot += std::min(std::abs(f[i]/m_atol[i]), std::abs(f[i]/y[i]/m_rtol));  // m_atol needs to be of size ny!
    }
    out[0] = tot/ny - m_special_settings[0];
    this->nrev++;
    return AnyODE::Status::success;
"""

_constr_special_settings = r"""
    if (m_special_settings.size() == 0){
         std::cerr << __FILE__ << ":" << __LINE__ << ": no special_settings passed, using default [%(factor)s]\n";
         m_special_settings = {%(factor)s};
    } else {
         // std::cerr << __FILE__ << ":" << __LINE__ << ": using special_settings:" << m_special_settings[0] << "\n";
    }
""" % {'factor': '1e2'}


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


def get_native(rsys, odesys, integrator, skip_keys=(0,), steady_state_root=False, conc_roots=None):
    comp_keys = Substance.composition_keys(rsys.substances.values(), skip_keys=skip_keys)
    if PartiallySolvedSystem is None:
        raise ValueError("Failed to import 'native_sys' from 'pyodesys.native'")
    elif isinstance(odesys, PartiallySolvedSystem):
        init_conc = '&m_p[%d]' % (len(odesys.params) - len(odesys.original_dep))
    else:
        init_conc = 'y'

    kw = dict(namespace_override={
        'p_get_dx_max': True,
    })
    if all(subst.composition is None for subst in rsys.substances.values()):
        pass
    else:
        kw['namespace_override']['p_anon'] = render_mako(
            _anon, odesys=odesys, ncomp=len(comp_keys),
            comp_conc=_get_comp_conc(rsys, odesys, comp_keys, skip_keys),
            subst_comp=_get_subst_comp(rsys, odesys, comp_keys, skip_keys))
        kw['namespace_override']['p_first_step'] = render_mako(
            _first_step, init_conc=init_conc, odesys=odesys)
    ns_extend = kw.get('namespace_extend', {})

    if steady_state_root or conc_roots:
        if not native_sys[integrator]._NativeCode._support_roots:
            raise ValueError("integrator '%s' does not support roots." % integrator)
        if odesys.roots is not None:
            raise ValueError("roots already set")
    if steady_state_root:
        assert conc_roots is None
        kw['namespace_override']['p_nroots'] = ' return 1; '
        kw['namespace_override']['p_roots'] = _roots_ss
        if 'p_constructor' not in ns_extend:
            ns_extend['p_constructor'] = []
        ns_extend['p_constructor'] += [_constr_special_settings]
    elif conc_roots:
        # This could (with some effort) be rewritten to take limits as parameters and have a
        # preprocessor in odesys.pre_processors do the dedimensionalization.
        assert not steady_state_root
        assert all(k in odesys.names and k in rsys.substances for k in conc_roots)
        kw['namespace_override']['p_nroots'] = ' return %d; ' % len(conc_roots)
        kw['namespace_override']['p_roots'] = (
            ''.join(['    out[%(i)d] = y[%(j)d] - m_special_settings[%(i)d];\n' %
                     dict(i=i, j=odesys.names.index(k)) for i, k in enumerate(conc_roots)]) +
            '    return AnyODE::Status::success;\n'
        )
        if 'p_constructor' not in ns_extend:
            ns_extend['p_constructor'] = []
        ns_extend['p_constructor'] += [
            'if (m_special_settings.size() != %d) throw std::runtime_error("special_settings missing");' %
            len(conc_roots)
        ]

    if 'p_includes' not in ns_extend:
        ns_extend['p_includes'] = set()
    ns_extend['p_includes'] |= {"<type_traits>",  "<vector>"}
    return native_sys[integrator].from_other(odesys, namespace_extend=ns_extend, **kw)
