# -*- coding: utf-8 -*-
"""
This module contains functions for formulating systems of Ordinary Differential
Equations (ODE-systems) which may be integrated numerically to model temporal
evolution of concentrations in reaction systems.
"""
from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict
from functools import reduce
from itertools import chain
from operator import mul
import math

try:
    import numpy as np
except ImportError:
    np = None

from ..units import to_unitless, get_derived_unit, rescale
from ..util._expr import Expr
from .rates import RateExpr, MassAction


def law_of_mass_action_rates(conc, rsys, variables=None):
    """ Returns a generator of reaction rate expressions

    Rates from the law of mass action (:attr:`Reaction.inact_reac` ignored)
    from a :class:`ReactionSystem`.

    Parameters
    ----------
    conc : array_like
        concentrations (floats or symbolic objects)
    rsys : ReactionSystem instance
        See :class:`ReactionSystem`
    variables : dict (optional)
        to override parameters in the rate expressions of the reactions

    Examples
    --------
    >>> from chempy import ReactionSystem, Reaction
    >>> line, keys = 'H2O -> H+ + OH- ; 1e-4', 'H2O H+ OH-'
    >>> rxn = Reaction.from_string(line, keys)
    >>> rsys = ReactionSystem([rxn], keys)
    >>> next(law_of_mass_action_rates([55.4, 1e-7, 1e-7], rsys))
    0.00554
    >>> from chempy.kinetics.rates import Arrhenius, MassAction
    >>> rxn.param = MassAction(Arrhenius({'A': 1e10, 'Ea_over_R': 9314}))
    >>> print('%.5g' % next(law_of_mass_action_rates([55.4, 1e-7, 1e-7], rsys, {'temperature': 293})))
    0.0086693

    """
    for idx_r, rxn in enumerate(rsys.rxns):
        if isinstance(rxn.param, RateExpr):
            if isinstance(rxn.param, MassAction):
                yield rxn.param(dict(chain(variables.items(), zip(rsys.substances.keys(), conc))), reaction=rxn)
            else:
                raise ValueError("Not mass-action rate in reaction %d" % idx_r)
        else:
            rate = 1
            for substance_key, coeff in rxn.reac.items():
                s_idx = rsys.as_substance_index(substance_key)
                rate *= conc[s_idx]**coeff
            yield rate*rxn.param


def dCdt_list(rsys, rates):
    """ Returns a list of the time derivatives of the concentrations

    Parameters
    ----------
    rsys: ReactionSystem instance
    rates: array_like
        rates (to be weighted by stoichiometries) of the reactions
        in ``rsys``

    Examples
    --------
    >>> from chempy import ReactionSystem, Reaction
    >>> line, keys = 'H2O -> H+ + OH- ; 1e-4', 'H2O H+ OH-'
    >>> rsys = ReactionSystem([Reaction.from_string(line, keys)], keys)
    >>> dCdt_list(rsys, [0.0054])
    [-0.0054, 0.0054, 0.0054]

    """
    f = [0]*rsys.ns
    net_stoichs = rsys.net_stoichs()
    for idx_s in range(rsys.ns):
        for idx_r in range(rsys.nr):
            f[idx_s] += net_stoichs[idx_r, idx_s]*rates[idx_r]
    return f


def get_odesys(rsys, include_params=True, substitutions=None, SymbolicSys=None, unit_registry=None,
               output_conc_unit=None, output_time_unit=None, cstr=False, **kwargs):
    """ Creates a :class:`pyneqsys.SymbolicSys` from a :class:`ReactionSystem`

    The parameters passed to RateExpr will contain the key ``'time'`` corresponding to the
    independent variable of the IVP.

    Parameters
    ----------
    rsys : ReactionSystem
        Each reaction of ``rsys`` will have their :meth:`Reaction.rate_expr()` invoked.
        Note that if :attr:`Reaction.param` is not a :class:`RateExpr` (or convertible to
        one through :meth:`as_RateExpr`) it will be used to construct a :class:`MassAction`
        instance.
    include_params : bool (default: True)
        Whether rate constants should be included into the rate expressions or
        left as free parameters in the :class:`pyneqsys.SymbolicSys` instance.
    substitutions : dict, optional
        Variable substitutions used by rate expressions (in respective Reaction.param).
        values are allowed to be values of instances of :class:`Expr`.
    SymbolicSys : class (optional)
        Default : :class:`pyneqsys.SymbolicSys`.
    unit_registry: dict (optional)
        See :func:`chempy.units.get_derived_units`.
    output_conc_unit : unit (Optional)
    output_time_unit : unit (Optional)
    cstr : bool
        Generate expressions for continuously stirred tank reactor.
    \*\*kwargs :
        Keyword arguemnts passed on to `SymbolicSys`.

    Returns
    -------
    pyodesys.symbolic.SymbolicSys
    extra : dict, with keys:
        - param_keys : list of str instances
        - unique : OrderedDict mapping str to value (possibly None)
        - p_units : list of units
        - max_euler_step_cb : callable or None
        - linear_dependencies : None or factory of solver callback
        - rate_exprs_cb : callable
        - cstr_fr_fc : None or (feed-ratio-key, subtance-key-to-feed-conc-key-map)

    Examples
    --------
    >>> from chempy import Equilibrium, ReactionSystem
    >>> eq = Equilibrium({'Fe+3', 'SCN-'}, {'FeSCN+2'}, 10**2)
    >>> substances = 'Fe+3 SCN- FeSCN+2'.split()
    >>> rsys = ReactionSystem(eq.as_reactions(kf=3.0), substances)
    >>> odesys, extra = get_odesys(rsys)
    >>> init_conc = {'Fe+3': 1.0, 'SCN-': .3, 'FeSCN+2': 0}
    >>> tout, Cout, info = odesys.integrate(5, init_conc)
    >>> Cout[-1, :].round(4)
    array([ 0.7042,  0.0042,  0.2958])

    """
    if SymbolicSys is None:
        from pyodesys.symbolic import SymbolicSys

    r_exprs = [rxn.rate_expr() for rxn in rsys.rxns]

    _ori_pk = set.union(*(ratex.all_parameter_keys() for ratex in r_exprs))
    _ori_uk = set.union(*(ratex.all_unique_keys() for ratex in r_exprs))
    _subst_pk = set()
    _active_subst = {}
    _passive_subst = {}
    substitutions = substitutions or {}

    unique = OrderedDict()
    unique_units = {}  # OrderedDict()

    cstr_fr_fc = (
        'feedratio',
        OrderedDict([(sk, 'fc_'+sk) for sk in rsys.substances])
    ) if cstr is True else cstr

    if cstr_fr_fc:
        _ori_pk.add(cstr_fr_fc[0])
        for k in cstr_fr_fc[1].values():
            _ori_pk.add(k)

    def _reg_unique_unit(k, arg_dim, idx):
        if unit_registry is None:
            return
        unique_units[k] = reduce(mul, [unit_registry[dim]**v for dim, v in arg_dim[idx].items()])

    def _reg_unique(expr, rxn=None):
        if not isinstance(expr, Expr):
            raise NotImplementedError("Currently only Expr sub classes are supported.")
        if unit_registry is None:
            arg_dim = None
        else:
            arg_dim = expr.args_dimensionality(reaction=rxn)
        if expr.args is None:
            for idx, k in enumerate(expr.unique_keys):
                if k not in substitutions:
                    unique[k] = None
                    _reg_unique_unit(k, arg_dim, idx)
        else:
            for idx, arg in enumerate(expr.args):
                if isinstance(arg, Expr):
                    _reg_unique(arg)
                elif expr.unique_keys is not None and idx < len(expr.unique_keys):
                    uk = expr.unique_keys[idx]
                    if uk not in substitutions:
                        unique[uk] = arg
                        _reg_unique_unit(uk, arg_dim, idx)

    for sk, sv in substitutions.items():
        if sk not in _ori_pk and sk not in _ori_uk:
            raise ValueError("Substitution: '%s' does not appear in any rate expressions." % sk)
        if isinstance(sv, Expr):
            _subst_pk.update(sv.parameter_keys)
            _active_subst[sk] = sv
            if not include_params:
                _reg_unique(sv)
        else:
            _passive_subst[sk] = sv
    all_pk = list(filter(lambda x: x not in substitutions and x != 'time', _ori_pk.union(_subst_pk)))

    if not include_params:
        for rxn, ratex in zip(rsys.rxns, r_exprs):
            _reg_unique(ratex, rxn)

    all_pk_with_unique = list(chain(all_pk, filter(lambda k: k not in all_pk, unique.keys())))
    if include_params:
        param_names_for_odesys = all_pk
    else:
        param_names_for_odesys = all_pk_with_unique

    if unit_registry is None:
        p_units = None
    else:
        # We need to make rsys_params unitless and create
        # a pre- & post-processor for SymbolicSys
        pk_units = [get_derived_unit(unit_registry, k) for k in all_pk]
        p_units = pk_units if include_params else (pk_units + [unique_units[k] for k in unique])
        new_r_exprs = []
        for ratex in r_exprs:
            if ratex.args is None:
                _new_ratex = ratex
            else:
                _pu, _new_ratex = ratex.dedimensionalisation(unit_registry)
            new_r_exprs.append(_new_ratex)
        r_exprs = new_r_exprs

        time_unit = get_derived_unit(unit_registry, 'time')
        conc_unit = get_derived_unit(unit_registry, 'concentration')

        def post_processor(x, y, p):
            time = x*time_unit
            if output_time_unit is not None:
                time = rescale(time, output_time_unit)
            conc = y*conc_unit
            if output_conc_unit is not None:
                conc = rescale(conc, output_conc_unit)
            return time, conc, np.array([elem*p_unit for elem, p_unit in zip(p.T, p_units)], dtype=object).T

        kwargs['to_arrays_callbacks'] = (
            lambda x: to_unitless(x, time_unit),
            lambda y: to_unitless(y, conc_unit),
            lambda p: np.array([to_unitless(px, p_unit) for px, p_unit in zip(
                p.T if hasattr(p, 'T') else p, p_units)]).T
        )
        kwargs['post_processors'] = kwargs.get('post_processors', []) + [post_processor]

    def dydt(t, y, p, backend=math):
        variables = dict(chain(y.items(), p.items()))
        if 'time' in variables:
            raise ValueError("Key 'time' is reserved.")
        variables['time'] = t
        for k, act in _active_subst.items():
            if unit_registry is not None:
                _, act = act.dedimensionalisation(unit_registry)
            variables[k] = act(variables, backend=backend)
        variables.update(_passive_subst)
        return rsys.rates(variables, backend=backend, ratexs=r_exprs, cstr_fr_fc=cstr_fr_fc)

    def reaction_rates(t, y, p, backend=math):
        variables = dict(chain(y.items(), p.items()))
        if 'time' in variables:
            raise ValueError("Key 'time' is reserved.")
        variables['time'] = t
        for k, act in _active_subst.items():
            if unit_registry is not None:
                _, act = act.dedimensionalisation(unit_registry)
            variables[k] = act(variables, backend=backend)
        variables.update(_passive_subst)
        return [ratex(variables, backend=backend, reaction=rxn) for
                rxn, ratex in zip(rsys.rxns, r_exprs)]

    names = [s.name for s in rsys.substances.values()]
    latex_names = [None if s.latex_name is None else ('$\\mathrm{' + s.latex_name + '}$')
                   for s in rsys.substances.values()]

    compo_vecs, compo_names = rsys.composition_balance_vectors()

    odesys = SymbolicSys.from_callback(
        dydt, dep_by_name=True, par_by_name=True, names=names,
        latex_names=latex_names, param_names=param_names_for_odesys,
        linear_invariants=None if len(compo_vecs) == 0 else compo_vecs,
        linear_invariant_names=None if len(compo_names) == 0 else compo_names,
        **kwargs)

    symbolic_ratexs = reaction_rates(
        odesys.indep, dict(zip(odesys.names, odesys.dep)),
        dict(zip(odesys.param_names, odesys.params)), backend=odesys.be)
    rate_exprs_cb = odesys._callback_factory(symbolic_ratexs)

    if rsys.check_balance(strict=True):
        # Composition available, we can provide callback for calculating
        # maximum allowed Euler forward step at start of integration.
        def max_euler_step_cb(x, y, p=()):
            _x, _y, _p = odesys.pre_process(*odesys.to_arrays(x, y, p))
            upper_bounds = rsys.upper_conc_bounds(_y)
            fvec = odesys.f_cb(_x[0], _y, _p)
            h = []
            for idx, fcomp in enumerate(fvec):
                if fcomp == 0:
                    h.append(float('inf'))
                elif fcomp > 0:
                    h.append((upper_bounds[idx] - _y[idx])/fcomp)
                else:  # fcomp < 0
                    h.append(-_y[idx]/fcomp)
            min_h = min(h)
            return min(min_h, 1)

        def linear_dependencies(preferred=None):
            if preferred is not None:
                if len(preferred) == 0:
                    raise ValueError("No preferred substance keys provided")
                if len(preferred) >= len(rsys.substances):
                    raise ValueError("Cannot remove all concentrations from linear dependencies")
                for k in preferred:
                    if k not in rsys.substances:
                        raise ValueError("Unknown substance key: %s" % k)

            def analytic_solver(x0, y0, p0, be):
                if preferred is None:
                    _preferred = None
                else:
                    _preferred = list(preferred)
                A = be.Matrix(compo_vecs)
                rA, pivots = A.rref()

                analytic_exprs = OrderedDict()
                for ri, ci1st in enumerate(pivots):
                    for idx in range(ci1st, odesys.ny):
                        key = odesys.names[idx]
                        if rA[ri, idx] == 0:
                            continue
                        if _preferred is None or key in _preferred:
                            terms = [rA[ri, di]*(odesys.dep[di] - y0[odesys.dep[di]])
                                     for di in range(ci1st, odesys.ny) if di != idx]
                            analytic_exprs[odesys[key]] = y0[odesys.dep[idx]] - sum(terms)/rA[ri, idx]
                            if _preferred is not None:
                                _preferred.remove(key)
                            break
                for k in reversed(list(analytic_exprs.keys())):
                    analytic_exprs[k] = analytic_exprs[k].subs(analytic_exprs)
                if _preferred is not None and len(_preferred) > 0:
                    raise ValueError("Failed to obtain analytic expression for: %s" % ', '.join(_preferred))
                return analytic_exprs

            return analytic_solver

    else:
        max_euler_step_cb = None
        linear_dependencies = None

    return odesys, {
        'param_keys': all_pk,
        'unique': unique,
        'p_units': p_units,
        'max_euler_step_cb': max_euler_step_cb,
        'linear_dependencies': linear_dependencies,
        'rate_exprs_cb': rate_exprs_cb,
        'cstr_fr_fc': cstr_fr_fc
    }


def chained_parameter_variation(odesys, durations, init_conc, varied_params, default_params, integrate_kwargs=None):
    """ Integrate an ODE-system for a serie of durations with some parameters changed in-between

    Parameters
    ----------
    odesys : :class:`pyodesys.ODESys` instance
    durations : iterable of floats
    init_conc : dict or array_like
    varied_params : dict mapping parameter name to array_like
        Each array_like need to be of same length as durations.
    default_params : dict or array_like
        Default values for the parameters of the ODE system.
    integrate_kwargs : dict
        Keyword arguments passed on to :meth:`pyodesys.ODESys.integrate`.

    """
    for k, v in varied_params.items():
        if len(v) != len(durations):
            raise ValueError("Mismathced lengths of durations and varied_params")
    integrate_kwargs = integrate_kwargs or {}
    touts = []
    couts = []
    infos = {}
    c0 = init_conc.copy()
    for idx, duration in enumerate(durations):
        params = default_params.copy()
        for k, v in varied_params.items():
            params[k] = v[idx]
        tout, cout, info = odesys.integrate(duration, c0, params, **integrate_kwargs)
        c0 = cout[-1, :]
        idx0 = 0 if idx == 0 else 1
        t_global = 0 if idx == 0 else touts[-1][-1]
        touts.append(tout[idx0:] + t_global)
        couts.append(cout[idx0:, ...])
        for k, v in info.items():
            if k.startswith('internal'):
                continue
            if k in infos:
                infos[k] += (v,)
            else:
                infos[k] = (v,)
    return np.concatenate(touts), np.concatenate(couts), infos
