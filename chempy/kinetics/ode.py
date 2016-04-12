# -*- coding: utf-8 -*-
"""
This module contains functions for formulating systems of Ordinary Differential
Equations (ODE-systems) which may be integrated numerically to model temporal
evolution of concentrations in reaction systems.
"""
from __future__ import (absolute_import, division, print_function)

from chempy.units import (
    get_derived_unit, to_unitless, default_unit_in_registry
)

from ..util.pyutil import deprecated

from .rates import _RateExpr, MassAction, law_of_mass_action_rates as _lomar

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


class _Always(object):
    __slots__ = ('value')

    def __init__(self, value):
        self.value = value

    def __getitem__(self, key):
        return self.value


def get_odesys(rsys, include_params=False, SymbolicSys=None,
               unit_registry=None, output_conc_unit=None,
               output_time_unit=None, state=None, **kwargs):
    """ Creates a :class:`pyneqsys.SymbolicSys` from a :class:`ReactionSystem`

    Parameters
    ----------
    rsys: ReactionSystem
    include_params: bool (default: False)
        whether rate constants should be included into the rate expressions or
        left as free parameters in the :class:`pyneqsys.SymbolicSys` instance.
    SymbolicSys: class (optional)
        default: :class:`pyneqsys.SymbolicSys`
    unit_registry: dict (optional)
        see :func:`chempy.units.get_derived_units`
    output_conc_unit: unit (Optional)
    output_time_unit: unit (Optional)
    state: object (optional)
        argument for reaction parameters
    \*\*kwargs:
        Keyword arguemnts pass on to `SymbolicSys`

    """
    if SymbolicSys is None:
        from pyodesys.symbolic import SymbolicSys

    if 'names' not in kwargs:
        kwargs['names'] = list(rsys.substances.keys())

    rate_exprs = []
    if unit_registry is not None:
        # We need to make rsys_params unitless and create
        # a post- & pre-processor for SymbolicSys

        rsys_params = []
        p_units = []
        for ri, rxn in enumerate(rsys.rxns):
            rate_expr = rxn.param
            if not isinstance(rate_expr, _RateExpr):
                rate_expr = MassAction([rate_expr])  # default
            _params = rate_expr.get_params()
            _p_units = [default_unit_in_registry(_, unit_registry) for _ in _params]
            p_units.extend(_p_units)
            _rsys_params = [to_unitless(p, unit) for p, unit in zip(_params, _p_units)]
            rsys_params.extend(_rsys_params)
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
    else:
        for ri, rxn in enumerate(rsys.rxns):
            param = rxn.param
            if isinstance(param, _RateExpr):
                rate_exprs.append(param)
            else:
                rate_exprs.append(MassAction([param]))
        rsys_params = rsys.params()

    def dydt(t, y, p):
        rates = []
        for ri, rate_expr in enumerate(rate_exprs):
            if unit_registry is None:
                rates.append(rate_expr.eval(rsys, ri, y, None))
            else:
                rates.append(rate_expr.eval(rsys, ri, y, None))
        return dCdt(rsys, rates)

    return SymbolicSys.from_callback(
        dydt, rsys.ns, 0 if include_params else rsys.nr, **kwargs)
