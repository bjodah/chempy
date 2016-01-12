# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from collections import namedtuple

from chempy.units import get_derived_unit, to_unitless

def law_of_mass_action_rates(y, rsys, k=None):
    for idx, rxn in enumerate(rsys.rxns):
        rate = 1
        for substance, coeff in rxn.reac.items():
            idx = rsys.as_substance_index(substance)
            rate *= y[idx]**coeff
        if k is None:
            yield rate * rxn.param
        else:
            yield rate * k[idx]


def dCdt(rsys, rates):
    f = [0]*rsys.ns
    net_stoichs = rsys.net_stoichs()
    for idx_s, sbstnc in enumerate(rsys.substances):
        for idx_r, rxn in enumerate(rsys.rxns):
            f[idx_s] += net_stoichs[idx_r, idx_s]*rates[idx_r]
    return f


class Always(object):
    __slots__ = ('value')
    def __init__(self, value):
        self.value = value

    def __getitem__(self, key):
        return self.value


def get_odesys(rsys, include_params=False, SymbolicSys=None,
               unit_registry=None):
    if SymbolicSys is None:
        from pyodesys.symbolic import SymbolicSys

    rsys_params = rsys.params()

    if unit_registry is not None:
        # We need to make rsys_params unit less and create
        # a post- & pre-processor for SymbolicSys
        time_unit = get_derived_unit(unit_registry, 'time')
        conc_unit = get_derived_unit(unit_registry, 'concentration')
        p_units = list(law_of_mass_action_rates(Always(1/conc_unit), rsys,
                                              Always(conc_unit/time_unit)))
        rsys_params = [to_unitless(elem, p_unit) for elem, p_unit
                       in zip(rsys_params, p_units)]
        def pre_processor(x, y, p):
            return to_unitless(x, time_unit), to_unitless(y, conc_unit), [
                to_unitless(elem, p_unit) for elem, p_unit in zip(p, p_units)]
        def post_processor(x, y, p):
            return x*time_unit, y*conc_unit, [elem*p_unit for elem, p_unit
                                              in zip(p, p_units)]
        kwargs = {'pre_processors': [pre_processor],
                  'post_processors': [post_processor]}
    else:
        kwargs = {}

    def dydt(t, y, p):
        rates = list(law_of_mass_action_rates(
            y, rsys, rsys_params if include_params else p))
        return dCdt(rsys, rates)

    return SymbolicSys.from_callback(
        dydt, rsys.ns, 0 if include_params else rsys.nr, **kwargs)
