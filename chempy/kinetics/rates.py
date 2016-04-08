# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from functools import reduce
from operator import mul, add

import numpy as np

from ..util.pyutil import defaultnamedtuple, deprecated


def _eval_k(k, state):
    is_func = callable(k) and not k.__class__.__name__ == 'Symbol'
    return k(state) if is_func else k


@deprecated()
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


# Tree structure
class _RateExpr(object):

    nargs = 0  # override

    def __init__(self, args, hard_args=None, ref=None):
        if self.nargs == 0:
            raise ValueError("nargs need to be set")
        if self.nargs == -1:
            self.nargs = len(args)
        if len(args) != self.nargs:
            raise ValueError("Incorrect number of parameters")
        self.args = args
        self.hard_args = hard_args  # cannot be overrided during integration
        self.ref = ref  # arbitrary placeholder
        self._nscalars = self._get_nscalars()
        self._accum_nscalars = np.cumsum([0] + self._nscalars)

    def _get_nscalars(self):
        nscalars = []
        for arg in self.args:
            if isinstance(arg, _RateExpr):
                nscalars.append(sum(arg._nscalars))
            else:
                nscalars.append(1)
        return nscalars

    def eval_arg(self, rsys, ri, concs, i, params=None):
        arg = self.args[i]
        if isinstance(arg, _RateExpr):
            if params is None:
                p = None
            else:
                p = params[self._accum_nscalars[i]:self._accum_nscalars[i+1]]
            return arg.eval(rsys, ri, concs, p)
        else:
            if params is None:
                return arg
            else:
                if len(params) != 1:
                    raise ValueError("Incorrect length of params")
                return params[0]

    def eval(self, rsys, ri, concs, params=None):
        raise NotImplementedError


class MassAction(_RateExpr):

    nargs = 1

    def eval(self, rsys, ri, concs, params=None):
        return self.eval_arg(rsys, ri, concs, 0, params) * reduce(
            mul, [concs[rsys.as_substance_index(k)]**v for
                  k, v in rsys.rxns[ri].reac.items()])


class GeneralPow(_RateExpr):

    nargs = 1

    def eval(self, rsys, ri, concs, params=None):
        return self.eval_arg(rsys, ri, concs, 0, params) * reduce(
            mul, [concs[rsys.as_substance_index(base_key)]**power for
                  base_key, power in self.hard_args.items()])


class Sum(_RateExpr):

    nargs = -1

    def eval(self, rsys, ri, concs, params=None):
        return reduce(add, [self.eval_arg(rsys, ri, concs, idx, params) for
                            idx in range(self.nargs)])


class Quotient(_RateExpr):

    nargs = 2

    def eval(self, rsys, ri, concs, params=None):
        return (self.eval_arg(rsys, ri, concs, 0, params) /
                self.eval_arg(rsys, ri, concs, 1, params))


class _MassAction(defaultnamedtuple('MassAction', 'k ref', [None])):

    def expr(self, ri, rsys, concs, state=None, params=None):
        if params is not None:
            param = params[0]
        else:
            param = self.k
        if hasattr(param, 'expr'):
            param = param.expr(ri, rsys, concs, state, params)
        return _eval_k(param, state) * reduce(
            mul, [concs[rsys.as_substance_index(k)]**v for
                  k, v in rsys.rxns[ri].reac.items()])


class _GeneralPow(defaultnamedtuple('GeneralPow', 'k powers ref', [None])):
    def expr(self, ri, rsys, concs, state=None, params=None):
        if params is None:
            rateconst, powers = self.k, self.powers
        else:
            rateconst, powers = params
        return _eval_k(rateconst, state) * reduce(
            mul, [concs[rsys.as_substance_index(base_key)]**power for
                  base_key, power in powers.items()])


class _Sum(defaultnamedtuple('Quotient', 'terms ref', [None])):
    def expr(self, ri, rsys, concs, state=None, params=None):
        if params is None:
            params = self.terms
        return reduce(add, [term.expr(ri, rsys, concs, state) for
                            term in params])


class _Quotient(defaultnamedtuple('Quotient', 'numer denom ref', [None])):
    def expr(self, ri, rsys, concs, state=None, params=None):
        if params is None:
            numerator, denominator = self.numer, self.denom
        else:
            numerator, denominator = params
        return (numerator.expr(ri, rsys, concs, state) /
                denominator.expr(ri, rsys, concs, state))
