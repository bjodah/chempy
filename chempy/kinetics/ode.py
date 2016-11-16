# -*- coding: utf-8 -*-
"""
This module contains functions for formulating systems of Ordinary Differential
Equations (ODE-systems) which may be integrated numerically to model temporal
evolution of concentrations in reaction systems.
"""
from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict
from itertools import chain
import math

from ..units import to_unitless, get_derived_unit, default_unit_in_registry
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
    >>> from chempy.kinetics.rates import ArrheniusMassAction
    >>> rxn.param = ArrheniusMassAction({'A': 1e10, 'Ea_over_R': 9314}, rxn=rxn)
    >>> print('%.5g' % next(law_of_mass_action_rates([55.4, 1e-7, 1e-7], rsys, {'temperature': 293})))
    0.0086693

    """
    for idx_r, rxn in enumerate(rsys.rxns):
        if isinstance(rxn.param, RateExpr):
            if isinstance(rxn.param, MassAction):
                yield rxn.param(dict(chain(variables.items(), zip(rsys.substances.keys(), conc))))
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
               output_conc_unit=None, output_time_unit=None, **kwargs):
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
    include_params : bool (default: False)
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
    \*\*kwargs :
        Keyword arguemnts pass on to `SymbolicSys`.

    Returns
    -------
    pyodesys.symbolic.SymbolicSys
    param_keys : list of str instances
    unique : OrderedDict mapping str to value (possibly None)
    p_units : list of units

    Examples
    --------
    >>> from chempy import Equilibrium, ReactionSystem
    >>> eq = Equilibrium({'Fe+3', 'SCN-'}, {'FeSCN+2'}, 10**2)
    >>> substances = 'Fe+3 SCN- FeSCN+2'.split()
    >>> rsys = ReactionSystem(eq.as_reactions(kf=3.0), substances)
    >>> odesys = get_odesys(rsys)[0]
    >>> init_conc = {'Fe+3': 1.0, 'SCN-': .3, 'FeSCN+2': 0}
    >>> tout, Cout, info = odesys.integrate(5, init_conc)
    >>> Cout[-1, :].round(4)
    array([ 0.7042,  0.0042,  0.2958])

    """
    if SymbolicSys is None:
        from pyodesys.symbolic import SymbolicSys

    r_exprs = [rxn.rate_expr() for rxn in rsys.rxns]

    _ori_pk = set.union(*(set(ratex.parameter_keys) for ratex in r_exprs))  # parameter_keys
    _subst_pk = set()
    _active_subst = {}
    _passive_subst = {}
    substitutions = substitutions or {}

    unique = OrderedDict()

    def _reg_unique(expr):
        if not isinstance(expr, Expr):
            raise NotImplementedError("Currently only Expr sub classes are supported.")
        if expr.args is None:
            for k in expr.unique_keys:
                unique[k] = None
        else:
            for idx, arg in enumerate(expr.args):
                if isinstance(arg, Expr):
                    _reg_unique(arg)
                elif expr.unique_keys is not None and idx < len(expr.unique_keys):
                    unique[expr.unique_keys[idx]] = arg

    for sk, sv in substitutions.items():
        if sk not in _ori_pk:
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
        for ratex in r_exprs:
            _reg_unique(ratex)

    all_pk_with_unique = list(chain(all_pk, unique.keys()))
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
        unique_units = [default_unit_in_registry(unit_registry, uv) for uv in unique.values()]
        p_units = pk_units if include_params else (pk_units + unique_units)
        new_r_exprs = []
        for ratex in r_exprs:
            _pu, _new_rate = ratex._recursive_as_RateExpr().dedimensionalisation(unit_registry)
            new_r_exprs.append(_new_rate)
        r_exprs = new_r_exprs

        time_unit = get_derived_unit(unit_registry, 'time')
        conc_unit = get_derived_unit(unit_registry, 'concentration')

        def pre_processor(x, y, p):
            return (
                to_unitless(x, time_unit),
                to_unitless(y, conc_unit),
                [to_unitless(px, p_unit) for px, p_unit in zip(p, p_units)]
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
        variables = dict(chain(y.items(), p.items()))
        variables['time'] = t
        for k, act in _active_subst.items():
            if unit_registry is not None:
                _, act = act.dedimensionalisation(unit_registry)
            variables[k] = act(variables, backend=backend)
        variables.update(_passive_subst)
        return rsys.rates(variables, backend=backend, ratexs=r_exprs)

    names = [s.name for s in rsys.substances.values()]
    latex_names = [s.latex_name for s in rsys.substances.values()]
    return SymbolicSys.from_callback(
        dydt, dep_by_name=True, par_by_name=True, names=names,
        latex_names=latex_names, param_names=param_names_for_odesys, **kwargs), all_pk, unique, p_units
