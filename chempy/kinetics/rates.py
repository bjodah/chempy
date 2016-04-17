# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from functools import reduce
import math
from operator import mul, add


def _accumulate(iterable):
    tot = 0
    for elem in iterable:
        tot += elem
        yield tot


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
        self._accum_nscalars = [0] + list(_accumulate(self._nscalars))

    def _sub_params(self, i, params):
        if params is None:
            return None
        else:
            return params[self._accum_nscalars[i]:self._accum_nscalars[i+1]]

    def rebuild(self, params=None):
        args = []
        for i, arg in enumerate(self.args):
            if isinstance(arg, _RateExpr):
                args.append(arg.rebuild(self._sub_params(i, params)))
            else:
                if params is None:
                    args.append(arg)
                else:
                    if len(params) != 1:
                        raise ValueError("Incorrect length of params")
                    args.append(params[0])

        return self.__class__(args, self.hard_args, self.ref)

    def _get_nscalars(self):  # how many scalars does this RateExpr consume
        nscalars = []
        for arg in self.args:
            if isinstance(arg, _RateExpr):
                nscalars.append(sum(arg._nscalars))
            else:
                nscalars.append(1)
        return nscalars

    def eval_arg(self, rsys, ri, concs, i, params=None, global_p=None,
                 backend=math):
        arg = self.args[i]
        if isinstance(arg, _RateExpr):
            return arg.eval(rsys, ri, concs, self._sub_params(i, params),
                            global_p, backend)
        else:
            if params is None:
                return arg
            else:
                if len(params) != 1:
                    raise ValueError("Incorrect length of params")
                return params[0]

    def eval(self, rsys, ri, concs, params=None, global_p=None, backend=math):
        raise NotImplementedError  # to be overrided

    def get_params(self):
        params = []
        for arg in self.args:
            if isinstance(arg, _RateExpr):
                params.extend(arg.get_params())
            else:
                params.append(arg)
        return params


class MassAction(_RateExpr):
    """ args[0] * prod(conc[i]**order(i), for all i in reactants) """
    nargs = 1

    def eval(self, rsys, ri, concs, params=None, global_p=None, backend=math):
        return (self.eval_arg(rsys, ri, concs, 0, params, global_p, backend) *
                reduce(mul, [concs[rsys.as_substance_index(k)]**v for
                             k, v in rsys.rxns[ri].reac.items()]))


class GeneralPow(_RateExpr):
    """ args[0] * prod(conc[idx(k)]**hard_args[k], all k in hard_args) """
    nargs = 1

    def eval(self, rsys, ri, concs, params=None, global_p=None, backend=math):
        # hard_args map substance to power
        return (self.eval_arg(rsys, ri, concs, 0, params, global_p, backend) *
                reduce(mul, [concs[rsys.as_substance_index(base_key)]**power
                             for base_key, power in self.hard_args.items()]))


class Sum(_RateExpr):
    """ sum(args)"""

    nargs = -1

    def eval(self, rsys, ri, concs, params=None, global_p=None, backend=math):
        return reduce(add, [
            self.eval_arg(rsys, ri, concs, idx, params, global_p, backend) for
            idx in range(self.nargs)])


class Quotient(_RateExpr):
    """ args[0] / args[1] """
    nargs = 2

    def eval(self, rsys, ri, concs, params=None, global_p=None, backend=math):
        return (self.eval_arg(rsys, ri, concs, 0, params, global_p, backend) /
                self.eval_arg(rsys, ri, concs, 1, params, global_p, backend))


class ExpReciprocalT(_RateExpr):
    """ args[0] * exp(args[1]/global_p[self.T]) """
    global_p_keys = ('T',)
    nargs = 2

    def eval(self, rsys, ri, concs, params=None, global_p=None, backend=math):
        return (self.eval_arg(rsys, ri, concs, 0, params, global_p, backend) *
                backend.exp(
                    self.eval_arg(rsys, ri, concs, 1, params, global_p,
                                  backend) / global_p[self.global_p_keys[0]]))
