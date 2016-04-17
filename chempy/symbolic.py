# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

from ._util import AttributeDict


def get_constant_symbols(Symbol=None):
    if Symbol is None:
        from sympy import Symbol
    consts = [
        ('Faraday_constant', 'F'),
        ('Avogadro_constant', 'N_A'),
        ('vacuum_permittivity', 'epsilon_0'),
        ('Boltzmann_constant', 'k_B'),
        ('pi', 'pi'),
        ('molar_gas_constant', 'R'),
    ]
    return AttributeDict({k: Symbol(v) for k, v in consts})


class SymbolRegistry(object):
    def __init__(self, backend='sympy', kw=None):
        try:
            self.backend = __import__(backend)
        except ImportError:
            self.backend = None
        self._reg = {}
        self._kw = kw or {}
        self.one = self.backend.S(1)

    def __getitem__(self, key):
        if key not in self._reg:
            self._reg[key] = self.backend.Symbol(key, **self._kw)
        return self._reg[key]

    def __contains__(self, key):
        return key in self._reg


concentrations = SymbolRegistry(kw=dict(positive=True))
state = SymbolRegistry()
parameters = SymbolRegistry()


def evaluate(expr, variables, registries=()):
    def _get(key):
        for reg in registries:
            if key in reg:
                return reg[key]
        raise KeyError("Unkown key: %s" % key)
    subsd = {_get(k): v for k, v in variables.items()}
    result = expr.subs(subsd)
    fs = result.free_symbols
    if len(fs) != 0:
        raise ValueError("%s missing in variables" % ', '.join(fs))
    return result
