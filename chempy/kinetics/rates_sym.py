# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from functools import reduce
import math
from operator import mul, add


class RateExpr(object):

    def __init__(self, args, rxn=None, arg_keys=None, ref=None):
        self.args = args
        self.rxn = rxn
        self.arg_keys = arg_keys
        self.ref = ref

    def __call__(self, variables, backend=math):
        raise NotImplementedError


class MassAction(RateExpr):

    state_keys = ()

    def __call__(self, variables, args=None, backend=math):
        args = args or self.args
        if self.arg_keys is None:
            result = args[0]
        else:
            result = variables.get(self.arg_keys[0], args[0])
        for k, v in self.rxn.reac.items():
            result *= variables[k]**v
        return result


class ArrheniusMassAction(MassAction):

    state_keys = ('temperature',)

    def __call__(self, variables, args=None, backend=math):
        A, Ea = args or self.args
        k = arrhenius_equation(A, Ea, T=variables['temperature'])
        return super(ArrheniusMassAction, self).__call__(
            variables, (k,), backend=backend)

if False:

    _fields = 'cb variable_names state_names parameter_names ref'
    Parameterisation = defaultnamedtuple('Parameterisation', _fields, [None])

    def _rate_param_counter():
        idx = 0
        while True:
            yield 'k%d' % idx
            idx += 1

    rate_parameter_dummies = _rate_param_counter()

    class MassAction(Parameterisation):

        def __new__(cls, const, reactants, ref=None):
            param_names = set([const]) if isinstance(const, str) else set()
            try:
                state_names = const.state_variables()
            except AttributeError:
                state_names = set()

            def cb(variables, backend=math):
                if isinstance(const, str):
                    result = variables[const]
                else:
                    try:
                        result = const(variables, backend=backend)
                    except TypeError:
                        result = const
                for k, v in reactants.items():
                    result *= variables[k]**v
                return result
            return super(MassAction, cls).__new__(
                cls, cb, reactants.keys(), state_names, param_names, ref)

    class ArrheniusMassAction(Parameterisation):

        def __new__(cls, A, Ea_over_R, reactants, ref=None):
            param_names =

    class _MassAction(object):

        def __init__(self, rate_const, rate_const_key=None, ref=None):
            if callable(rate_const):
                self.rate_const = rate_expr
            else:
                self.rate_const = lambda state: rate_const
            self.rate_const_key = rate_const_key
            self.ref = None

        def __call__(self, conc, state, rxn_params, backend=math):
            raise NotImplementedError

        def state_variables(self):
            try:
                return self.rate_const.state_variables()
            except AttributeError:
                return set()

        def from_rxn(self, rxn, rate_const_keys=None):
            rck = rate_const_keys or (self.rate_const_key or (None,))

            def cb(conc, state, rxn_params, backend=math):
                result = rxn_params.get(rck[0], self.rate_const(state))
                for k, v in rxn.reac.items():
                    result *= conc[k]**v
                return result

    def evaluate(rxn_param, rxn=None):
        try:
            return rxn_param.from_rxn(rnx, )

    def _accumulate(iterable):
        tot = 0
        for elem in iterable:
            tot += elem
            yield tot

    def _eval_k(k, state):
        is_func = callable(k) and not k.__class__.__name__ == 'Symbol'
        return k(state) if is_func else k

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
            self.hard_args = hard_args  # cannot be overrided during integrat.
            self.ref = ref  # arbitrary placeholder
            self._nscalars = self._get_nscalars()
            self._accum_nscalars = [0] + list(_accumulate(self._nscalars))

        def _sub_params(self, i, params):
            if params is None:
                return None
            else:
                return params[self._accum_nscalars[i]:
                              self._accum_nscalars[i+1]]

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

        def eval(self, rsys, ri, concs, params=None, global_p=None,
                 backend=math):
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

        def eval(self, rsys, ri, concs, params=None, global_p=None,
                 backend=math):
            return (self.eval_arg(rsys, ri, concs, 0, params, global_p,
                                  backend) *
                    reduce(mul, [concs[rsys.as_substance_index(k)]**v for
                                 k, v in rsys.rxns[ri].reac.items()]))

    class GeneralPow(_RateExpr):
        """ args[0] * prod(conc[idx(k)]**hard_args[k], all k in hard_args) """
        nargs = 1

        def eval(self, rsys, ri, concs, params=None, global_p=None,
                 backend=math):
            # hard_args map substance to power
            return (self.eval_arg(rsys, ri, concs, 0, params, global_p,
                                  backend) *
                    reduce(mul, [
                        concs[rsys.as_substance_index(base_key)]**power
                        for base_key, power in self.hard_args.items()]))

    class Sum(_RateExpr):
        """ sum(args)"""

        nargs = -1

        def eval(self, rsys, ri, concs, params=None, global_p=None,
                 backend=math):
            return reduce(add, [
                self.eval_arg(rsys, ri, concs, idx, params, global_p, backend)
                for idx in range(self.nargs)])

    class Quotient(_RateExpr):
        """ args[0] / args[1] """
        nargs = 2

        def eval(self, rsys, ri, concs, params=None, global_p=None,
                 backend=math):
            return (self.eval_arg(rsys, ri, concs, 0, params, global_p,
                                  backend) /
                    self.eval_arg(rsys, ri, concs, 1, params, global_p,
                                  backend))

    class ExpReciprocalT(_RateExpr):
        """ args[0] * exp(args[1]/global_p[self.T]) """
        global_p_keys = ('T',)
        nargs = 2

        def eval(self, rsys, ri, concs, params=None, global_p=None,
                 backend=math):
            return (self.eval_arg(rsys, ri, concs, 0, params, global_p,
                                  backend) *
                    backend.exp(
                        self.eval_arg(rsys, ri, concs, 1, params, global_p,
                                      backend) /
                        global_p[self.global_p_keys[0]]))
