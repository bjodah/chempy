# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

try:
    import numpy as np
except ImportError:
    np = None

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

from .. import ReactionSystem, Equilibrium
from ..units import get_derived_unit, to_unitless, default_units as u


def _dominant_reaction_effects(substance_key, rsys, rates, linthreshy, eqk1, eqk2, eqs):
    tot = np.zeros(rates.shape[0])
    reaction_effects = rsys.per_reaction_effect_on_substance(substance_key)
    data = []
    for ri, n in reaction_effects.items():
        tot += n*rates[..., ri]
        if ri in eqk1:
            otheri = eqk2[eqk1.index(ri)]
            y = n*rates[..., ri] + reaction_effects[otheri]*rates[..., otheri]
            rxn = eqs[eqk1.index(ri)]
        elif ri in eqk2:
            continue
        else:
            y = n*rates[..., ri]
            rxn = rsys.rxns[ri]
        if np.all(np.abs(y) < linthreshy):
            continue
        data.append((y, rxn))
    return data, tot


def _combine_rxns_to_eq(rsys):
    eqk1, eqk2 = zip(*rsys.identify_equilibria())
    eqs = [Equilibrium(
        rsys.rxns[i1].reac,
        rsys.rxns[i1].prod,
        (rsys.rxns[i1].param, rsys.rxns[i2].param),
        inact_reac=rsys.rxns[i1].inact_reac,
        inact_prod=rsys.rxns[i1].inact_prod
    ) for i1, i2 in zip(eqk1, eqk2)]
    return eqk1, eqk2, eqs


def plot_reaction_contributions(
        xyp, rsys, rate_exprs_cb, substance_keys=None, varied=None, axes=None,
        total=False, linthreshy=1e-9, relative=False, xscale='log', yscale='symlog',
        xlabel='Time', ylabel=None, combine_equilibria=False, selection=slice(None),
        unit_registry=None):
    """ Plots per reaction contributions to concentration evolution of a substance.

    Parameters
    ----------
    xyp : ``pyodesys.results.Result`` instance or length 3 tuple or xout,yout,params
    result : pyodesys.results.Result
    substance_key : str

    """
    from pyodesys.results import Result
    if isinstance(xyp, Result):
        xyp = xyp.odesys.to_arrays(xyp.xout, xyp.yout, xyp.params, reshape=False)
    if varied is None:
        varied = xyp[0]
    if xyp[1].shape[-2] != varied.size:
        raise ValueError("Size mismatch between varied and yout")
    if substance_keys is None:
        substance_keys = rsys.substances.keys()
    if axes is None:
        _fig, axes = plt.subplots(len(substance_keys))
    rates = rate_exprs_cb(*xyp)
    if unit_registry is not None:
        time_unit = get_derived_unit(unit_registry, 'time')
        conc_unit = get_derived_unit(unit_registry, 'concentration')
        rates = to_unitless(rates*conc_unit/time_unit, u.molar/u.second)

    eqk1, eqk2, eqs = _combine_rxns_to_eq(rsys) if combine_equilibria else ([], [], [])

    for sk, ax in zip(substance_keys, axes):
        data, tot = _dominant_reaction_effects(sk, rsys, rates, linthreshy, eqk1, eqk2, eqs)
        factor = 1/xyp[1][:, rsys.as_substance_index(sk)] if relative else 1
        if total:
            ax.plot(varied, factor*tot, c='k', label='Total', linewidth=2, ls=':')
        for y, rxn in sorted(data, key=lambda args: args[0][-1], reverse=True):
            ax.plot(varied, factor*y,
                    label=r'$\mathrm{%s}$' % rxn.latex(rsys.substances))

        if rsys.substances[sk].latex_name is None:
            ttl = rsys.substances[sk].name
            ttl_template = '%s'
        else:
            ttl = rsys.substances[sk].latex_name
            ttl_template = r'\mathrm{$%s$}'

        if yscale == 'symlog':
            ax.axhline(linthreshy, linestyle='--', color='k', linewidth=.5)
            ax.axhline(-linthreshy, linestyle='--', color='k', linewidth=.5)
            ax.set_yscale(yscale, linthreshy=linthreshy)
        else:
            ax.set_yscale(yscale)

        if ylabel is None:
            ax.set_ylabel(r'$\frac{d}{dt}\left[%s\right]\ /\ M\cdot s^{-1}$' % ttl)
        else:
            ax.set_ylabel(ylabel)
            ax.set_title(ttl_template % ttl)

        ax.set_xlabel(xlabel)
        ax.set_xscale(xscale)
        ax.legend(loc='best')


def dominant_reactions_graph(concs, rate_exprs_cb, rsys, substance_key, linthreshy=1e-9,
                             fname='dominant_reactions_graph.png', relative=False,
                             combine_equilibria=False, **kwargs):
    from ..util.graph import rsys2graph
    rates = rate_exprs_cb(0, concs)
    eqk1, eqk2, eqs = _combine_rxns_to_eq(rsys) if combine_equilibria else ([], [], [])
    rrate, rxns = zip(*_dominant_reaction_effects(
        substance_key, rsys, rates, linthreshy, eqk1, eqk2, eqs)[0])
    rsys = ReactionSystem(rxns, rsys.substances, rsys.name)
    lg_rrate = np.log10(np.abs(rrate))
    rsys2graph(rsys, fname=fname, penwidths=1 + lg_rrate - np.min(lg_rrate), **kwargs)
