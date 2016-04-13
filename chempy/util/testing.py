# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import pytest


class requires(object):

    def __init__(self, *pymod):
        self.missing = []
        for name in pymod:
            try:
                __import__(name)
            except ImportError:
                self.missing.append(name)

    def __call__(self, cb):
        reason = "missing modules %s" % ', '.join(self.missing)
        return pytest.mark.skipif(self.missing, reason=reason)(cb)
