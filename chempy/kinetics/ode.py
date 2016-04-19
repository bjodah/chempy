# -*- coding: utf-8 -*-
"""
This module contains functions for formulating systems of Ordinary Differential
Equations (ODE-systems) which may be integrated numerically to model temporal
evolution of concentrations in reaction systems.
"""
from __future__ import (absolute_import, division, print_function)

from itertools import chain
import math

from ..units import default_unit_in_registry, to_unitless, get_derived_unit
from ..util.pyutil import deprecated
from .rates import RateExpr, MassAction, law_of_mass_action_rates as _lomar


law_of_mass_action_rates = deprecated(
    use_instead='.rates.law_of_mass_action_rates')(_lomar)


def dCdt(rsys, rates):
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
    >>> dCdt(rsys, [0.0054])
    [-0.0054, 0.0054, 0.0054]

    """
    f = [0]*rsys.ns
    net_stoichs = rsys.net_stoichs()
    for idx_s in range(rsys.ns):
        for idx_r in range(rsys.nr):
            f[idx_s] += net_stoichs[idx_r, idx_s]*rates[idx_r]
    return f


def get_odesys(rsys, include_params=False, global_params=None,
               SymbolicSys=None,
               unit_registry=None, output_conc_unit=None,
               output_time_unit=None, **kwargs):
    """ Creates a :class:`pyneqsys.SymbolicSys` from a :class:`ReactionSystem`

    Parameters
    ----------
    rsys : ReactionSystem
    include_params : bool (default: False)
        whether rate constants should be included into the rate expressions or
        left as free parameters in the :class:`pyneqsys.SymbolicSys` instance.
    global_params : dict, (optional)
        shared state used by rate expressions (in respective Reaction.param).
    SymbolicSys : class (optional)
        default : :class:`pyneqsys.SymbolicSys`
    unit_registry: dict (optional)
        see :func:`chempy.units.get_derived_units`
    output_conc_unit : unit (Optional)
    output_time_unit : unit (Optional)
    \*\*kwargs :
        Keyword arguemnts pass on to `SymbolicSys`

    """
    if SymbolicSys is None:
        from pyodesys.symbolic import SymbolicSys

    substance_keys = list(rsys.substances.keys())

    if 'names' not in kwargs:
        kwargs['names'] = substance_keys

    if unit_registry is None:
        def _param(rxn):
            if isinstance(rxn.param, RateExpr):
                return rxn.param
            else:
                return MassAction([rxn.param], rxn)
        r_exprs = [_param(rxn) for rxn in rsys.rxns]
    else:
        # We need to make rsys_params unitless and create
        # a post- & pre-processor for SymbolicSys
        r_exprs = []
        p_units = []
        for ri, rxn in enumerate(rsys.rxns):
            rate_expr = rxn.param
            if not isinstance(rate_expr, RateExpr):
                rate_expr = MassAction([rate_expr], rxn)  # default
            _p = rate_expr.args

            _pu = [default_unit_in_registry(_, unit_registry) for _ in _p]
            p_units.extend(_pu)
            _rsys_params = [to_unitless(p, unit) for p, unit in zip(_p, _pu)]
            r_exprs.append(rate_expr.__class__(
                _rsys_params, rxn, rate_expr.arg_keys, rate_expr.ref))

    state_keys = list(set.union(*(set(ratex.state_keys) for ratex in r_exprs)))

    param_keys = []
    p_defaults = []
    if not include_params:
        for ratex in r_exprs:
            if ratex.arg_keys is not None:
                param_keys.extend(ratex.arg_keys)
                p_defaults.extend(ratex.args)

    # param_keys = chain(ratex.arg_keys for ratex in r_exprs)
    # p_defaults = chain(ratex.args for ratex in r_exprs)

    if unit_registry is None:
        def pre_processor(x, y, p):
            return (
                x,
                rsys.as_per_substance_array(y),
                [p[k] for k in state_keys] + [p[k] for k in param_keys]
            )

        def post_processor(x, y, p):
            return (
                x,
                y,  # dict(zip(substance_keys, y)),
                dict(zip(state_keys+param_keys, p))
            )
    else:
        time_unit = get_derived_unit(unit_registry, 'time')
        conc_unit = get_derived_unit(unit_registry, 'concentration')

        def pre_processor(x, y, p):
            return (
                to_unitless(x, time_unit),
                rsys.as_per_substance_array(to_unitless(y, conc_unit)),
                [to_unitless(elem, p_unit) for elem, p_unit in zip(p, p_units)]
            )

        def post_processor(x, y, p):
            time = x*time_unit
            if output_time_unit is not None:
                time = time.rescale(output_time_unit)
            conc = y*conc_unit
            if output_conc_unit is not None:
                conc = conc.rescale(output_conc_unit)
            return time, conc, [elem*p_unit for elem, p_unit
                                in zip(p, p_units)]

    kwargs['pre_processors'] = [pre_processor] + kwargs.get('pre_processors', [])
    kwargs['post_processors'] = kwargs.get('post_processors', []) + [post_processor]

    def dydt(t, y, p, backend=math):
        variables = dict(chain(
            zip(substance_keys, y),
            zip(state_keys, p[:len(state_keys)]),
            zip(param_keys, p[len(state_keys):])
        ))
        return dCdt(rsys, [rat(variables, backend=backend) for rat in r_exprs])

    return SymbolicSys.from_callback(
        dydt, len(substance_keys),
        len(state_keys) + (0 if include_params else len(param_keys)),
        **kwargs)
