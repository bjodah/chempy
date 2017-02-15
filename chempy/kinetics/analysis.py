# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from itertools import chain

try:
    import numpy as np
except ImportError:
    np = None

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

from .. import Equilibrium


def plot_reaction_contributions(varied, concs, rate_exprs_cb, rsys, substance_keys=None, axes=None,
                                linthreshy=1e-9, relative=False, xscale='log', yscale='symlog',
                                xlabel='Time', ylabel=None, combine_equilibria=False):
    """ Plots per reaction contributions to concentration evolution of a substance.

    Parameters
    ----------
    result : pyodesys.results.Result
    substance_key : str

    """
    if concs.shape[0] != varied.size:
        raise ValueError("Size mismatch between varied and concs")
    if substance_keys is None:
        substance_keys = rsys.substances.keys()
    if axes is None:
        _fig, axes = plt.subplots(len(substance_keys))
    rates = rate_exprs_cb(varied, concs)

    if combine_equilibria:
        eqk1, eqk2 = zip(*rsys.identify_equilibria())
        eqs = [Equilibrium(rsys.rxns[ri].reac,
                           rsys.rxns[ri].prod,
                           inact_reac=rsys.rxns[ri].inact_reac,
                           inact_prod=rsys.rxns[ri].inact_prod) for ri in eqk1]
    else:
        eqk1, eqk2, eqs = [], [], []

    for sk, ax in zip(substance_keys, axes):
        si = rsys.as_substance_index(sk)
        reaction_effects = rsys.per_reaction_effect_on_substance(sk)
        for ri, n in reaction_effects.items():
            if ri in eqk1:
                otheri = eqk2[eqk1.index(ri)]
                y = n*rates[:, ri] + reaction_effects[otheri]*rates[:, otheri]
                lbl = r'$\mathrm{%s}$' % eqs[eqk1.index(ri)].latex(rsys.substances)
            elif ri in eqk2:
                continue
            else:
                y = n*rates[:, ri]
                lbl = r'$\mathrm{%s}$' % rsys.rxns[ri].latex(rsys.substances)

            if relative:
                y /= concs[:, si]
            if np.all(np.abs(y) < linthreshy):
                continue
            ax.plot(varied, y, label=lbl)
        if rsys.substances[sk].latex_name is None:
            ttl = rsys.substances[sk].name
            ttl_template = '%s'
        else:
            ttl = rsys.substances[sk].latex_name
            ttl_template = r'\mathrm{$%s$}'

        if yscale == 'symlog':
            #ax.hlines([linthreshy, -linthreshy], 0, 1, transform=ax.get_yaxis_transform(), linestyle='--', color='k')
            ax.axhline(linthreshy, linestyle='--', color='k')
            ax.axhline(-linthreshy, linestyle='--', color='k')
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
