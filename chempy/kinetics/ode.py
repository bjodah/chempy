# -*- coding: utf-8 -*-
"""
This module contains functions for formulating systems of Ordinary Differential
Equations (ODE-systems) which may be integrated numerically to model temporal
evolution of concentrations in reaction systems.
"""
from __future__ import (absolute_import, division, print_function)

from itertools import chain
import math

from ..units import to_unitless, get_derived_unit
from ..util.pyutil import deprecated
from ..util.expr import Expr
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


def get_odesys(rsys, include_params=False, substitutions=None,
               SymbolicSys=None,
               unit_registry=None, output_conc_unit=None,
               output_time_unit=None, **kwargs):
    """ Creates a :class:`pyneqsys.SymbolicSys` from a :class:`ReactionSystem`

    Parameters
    ----------
    rsys : ReactionSystem
        note that if :attr:`param` if not RateExpr it will be inspected for
        :meth:`_as_RateExpr`, lacking such it will be used to construct a
        :class:`MassAction` instance.
    include_params : bool (default: False)
        whether rate constants should be included into the rate expressions or
        left as free parameters in the :class:`pyneqsys.SymbolicSys` instance.
    substitutions : dict, optional
        variable substitutions used by rate expressions (in respective Reaction.param).
        values are allowed to be tuple like: (new_vars, callback)
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
        kwargs['names'] = list(rsys.substances.values())  # pyodesys>=0.5.3

    def _param(rxn):
        if isinstance(rxn.param, RateExpr):
            return rxn.param
        else:
            try:
                return rxn.param._as_RateExpr(rxn)
            except AttributeError:
                return MassAction([rxn.param], rxn=rxn)
    r_exprs = [_param(rxn) for rxn in rsys.rxns]

    _original_param_keys = set.union(*(set(ratex.parameter_keys) for ratex in r_exprs))
    _from_subst = set()
    _active_subst = {}
    _passive_subst = {}
    substitutions = substitutions or {}
    for key, v in substitutions.items():
        if key not in _original_param_keys:
            raise ValueError("Substitution: '%s' does not appear in any rate expressions.")
        if isinstance(v, Expr):
            _from_subst.update(v.parameter_keys)
            _active_subst[key] = v
        else:
            _passive_subst[key] = v
    param_keys = list(filter(lambda x: x not in substitutions, _original_param_keys.union(_from_subst)))

    arg_keys = []
    p_defaults = []
    if not include_params:
        for ratex in r_exprs:
            if ratex.arg_keys is not None:
                arg_keys.extend(ratex.arg_keys)
                p_defaults.extend(ratex.args)

    # arg_keys = chain(ratex.arg_keys for ratex in r_exprs)
    # p_defaults = chain(ratex.args for ratex in r_exprs)

    if unit_registry is None:
        def pre_processor(x, y, p):
            return (
                x,
                rsys.as_per_substance_array(y),
                [p[k] for k in param_keys] + [p[k] for k in arg_keys]
            )

        def post_processor(x, y, p):
            return (
                x,
                y,  # dict(zip(substance_keys, y)),
                dict(zip(param_keys+arg_keys, p))
            )
    else:
        # We need to make rsys_params unitless and create
        # a pre- & post-processor for SymbolicSys
        p_units = [get_derived_unit(unit_registry, k) for k in param_keys]
        new_r_exprs = []
        for ratex in r_exprs:
            _pu, _new_rates = ratex._dedimensionalisation(unit_registry)
            p_units.extend(_pu)
            new_r_exprs.append(_new_rates)
        r_exprs = new_r_exprs

        time_unit = get_derived_unit(unit_registry, 'time')
        conc_unit = get_derived_unit(unit_registry, 'concentration')

        def pre_processor(x, y, p):
            return (
                to_unitless(x, time_unit),
                rsys.as_per_substance_array(to_unitless(y, conc_unit)),
                # [to_unitless(elem, p_unit) for elem, p_unit in zip(p, p_units)]
                [to_unitless(p[k], p_unit) for k, p_unit in zip(chain(param_keys, arg_keys), p_units)]
            )

        def post_processor(x, y, p):
            time = x*time_unit
            if output_time_unit is not None:
                time = time.rescale(output_time_unit)
            conc = y*conc_unit
            if output_conc_unit is not None:
                conc = conc.rescale(output_conc_unit)
            return time, conc, [elem*p_unit for elem, p_unit in zip(p, p_units)]

    kwargs['pre_processors'] = [pre_processor] + kwargs.get('pre_processors', [])
    kwargs['post_processors'] = kwargs.get('post_processors', []) + [post_processor]

    def dydt(t, y, p, backend=math):
        variables = dict(chain(
            zip(substance_keys, y),
            zip(param_keys, p[:len(param_keys)]),
            zip(arg_keys, p[len(param_keys):])
        ))
        for k, act in _active_subst.items():
            if unit_registry is not None:
                _, act = act._dedimensionalisation(unit_registry)
            variables[k] = act(variables, backend=backend)
        variables.update(_passive_subst)
        return dCdt(rsys, [rat(variables, backend=backend) for rat in r_exprs])

    return SymbolicSys.from_callback(
        dydt, len(substance_keys),
        len(param_keys) + (0 if include_params else len(arg_keys)),
        **kwargs)
