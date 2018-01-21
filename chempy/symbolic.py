# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

from .util.pyutil import AttributeContainer


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
    return AttributeContainer(**{k: Symbol(v) for k, v in consts})
