# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from operator import mul, add
from functools import reduce

from ..util.pyutil import defaultnamedtuple


def _eval_k(k, state):
    is_func = callable(k) and not k.__class__.__name__ == 'Symbol'
    return k(state) if is_func else k


def law_of_mass_action_rates(conc, rsys, params=None, state=None):
    """ Returns a generator of reaction rate expressions

    Rates according to the law of mass action (:attr:`rsys.inact_reac` ignored)
    from a :class:`ReactionSystem`.

    Parameters
    ----------
    conc: array_like
        concentrations (floats or symbolic objects)
    rsys: ReactionSystem instance
        See :class:`ReactionSystem`
    params: array_like (optional)
        to override rate parameters of the reactions
    state: object (optional)
        argument for reaction parameters

    Examples
    --------
    >>> from chempy import ReactionSystem, Reaction
    >>> line, keys = 'H2O -> H+ + OH- ; 1e-4', 'H2O H+ OH-'
    >>> rsys = ReactionSystem([Reaction.from_string(line, keys)], keys)
    >>> next(law_of_mass_action_rates([55.4, 1e-7, 1e-7], rsys))
    0.00554

    """
    for rxn_idx, rxn in enumerate(rsys.rxns):
        rate = 1
        for substance_key, coeff in rxn.reac.items():
            s_idx = rsys.as_substance_index(substance_key)
            rate *= conc[s_idx]**coeff
        if params is None:
            yield rate * _eval_k(rxn.param, state)
        else:
            yield rate * _eval_k(params[rxn_idx], state)


class MassAction(defaultnamedtuple('MassAction', 'k ref', [None])):

    def expr(self, ri, rsys, concs, state=None, params=None):
        if params is not None:
            param = params[0]
        else:
            param = self.k
        return _eval_k(param, state) * reduce(
            mul, [concs[rsys.as_substance_index(k)]**v for
                  k, v in rsys.rxns[ri].reac.items()])


class GeneralPow(defaultnamedtuple('GeneralPow', 'k powers ref', [None])):
    def expr(self, ri, rsys, concs, state=None, params=None):
        if params is None:
            rateconst, powers = self.k, self.powers
        else:
            rateconst, powers = params
        return _eval_k(rateconst, state) * reduce(
            mul, [concs[rsys.as_substance_index(base_key)]**power for
                  base_key, power in powers.items()])


class Sum(defaultnamedtuple('Quotient', 'terms ref', [None])):
    def expr(self, ri, rsys, concs, state=None, params=None):
        if params is None:
            params = self.terms
        return reduce(add, [term.expr(ri, rsys, concs, state) for
                            term in params])


class Quotient(defaultnamedtuple('Quotient', 'numer denom ref', [None])):
    def expr(self, ri, rsys, concs, state=None, params=None):
        if params is None:
            numerator, denominator = self.numer, self.denom
        else:
            numerator, denominator = params
        return (numerator.expr(ri, rsys, concs, state) /
                denominator.expr(ri, rsys, concs, state))
