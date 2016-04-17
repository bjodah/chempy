# -*- coding: utf-8 -*-
"""
This module contains functions for formulating systems of Ordinary Differential
Equations (ODE-systems) which may be integrated numerically to model temporal
evolution of concentrations in reaction systems.
"""
from __future__ import (absolute_import, division, print_function)

from .ode import dCdt

def get_odesys(rsys, include_params=False, global_params=None,
               SymbolicSys=None,
               unit_registry=None, output_conc_unit=None,
               output_time_unit=None, state=None, **kwargs):
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
    state : object (optional)
        argument for reaction parameters
    \*\*kwargs :
        Keyword arguemnts pass on to `SymbolicSys`

    """
    if SymbolicSys is None:
        from pyodesys.symbolic import SymbolicSys

    substance_keys = list(rsys.substances.keys())

    if 'names' not in kwargs:
        kwargs['names'] = substance_keys

    if unit_registry is None:
        rck = None if include_params else rsys.
        rate_exprs = [rxn._get_param_cb(rck) for rxn in rsys.rxns]
    else:
        # We need to make rsys_params unitless and create
        # a post- & pre-processor for SymbolicSys
        rate_exprs = []
        p_units = []
        for ri, rxn in enumerate(rsys.rxns):
            rate_expr = rxn.param
            if not isinstance(rate_expr, _RateExpr):
                rate_expr = MassAction([rate_expr])  # default
            _params = rate_expr.get_params()
            _p_units = [default_unit_in_registry(_, unit_registry) for _ in _params]
            p_units.extend(_p_units)
            _rsys_params = [to_unitless(p, unit) for p, unit in zip(_params, _p_units)]
            rate_exprs.append(rate_expr.rebuild(_rsys_params))

        time_unit = get_derived_unit(unit_registry, 'time')
        conc_unit = get_derived_unit(unit_registry, 'concentration')

        def pre_processor(x, y, p):
            return to_unitless(x, time_unit), to_unitless(y, conc_unit), [
                to_unitless(elem, p_unit) for elem, p_unit in zip(p, p_units)]

        def post_processor(x, y, p):
            time = x*time_unit
            if output_time_unit is not None:
                time = time.rescale(output_time_unit)
            conc = y*conc_unit
            if output_conc_unit is not None:
                conc = conc.rescale(output_conc_unit)
            return time, conc, [elem*p_unit for elem, p_unit
                                in zip(p, p_units)]
        kwargs['pre_processors'] = [pre_processor]
        kwargs['post_processors'] = [post_processor]


    param_keys = chain(ratex.arg_keys for ratex in rate_exprs)
    param_defaults = chain(ratex.args for ratex in rate_exprs)

    state_keys = set.union(*(ratex.state_keys for ratex in rate_exprs))
    nstate = len(state_keys)

    def dydt(t, y, p, backend=math):
        rates = []
        variables = dict(chain(
            zip(substance_keys, y),
            zip(state_keys, p[:len(state_keys)]),
            zip(param_keys, param_defaults) if include_params else zip(param_keys, p[len(state_keys):])
        ))
        return dCdt(rsys, [rat(variables, backend=backend) for rat in rate_exprs])

    return SymbolicSys.from_callback(
        dydt, len(substance_keys), len(state_keys) + (0 if include_params else len(param_keys)), **kwargs)
