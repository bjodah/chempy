# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


def ionic_strength(molalities, charges):
    tot = 0
    for b, z in zip(molalities, charges):
        tot += b*z**2
    return tot/2


class ActivityProduct(object):

    def __init__(self, stoich, *args):
        self.stoich = stoich
        self.args = args

    def __call__(self, c):
        pass
