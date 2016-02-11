#!/usr/bin/env python
# -*- coding: utf-8 -*-

from chempy.chemistry import Substance
assert Substance.from_formula('H2O').composition == {1: 2, 8: 1}
