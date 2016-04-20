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
