# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from pkg_resources import parse_requirements, parse_version

from operator import lt, le, eq, ne, ge, gt
import pytest

_relop = dict(zip('< <= == != >= >'.split(), (lt, le, eq, ne, ge, gt)))


class requires(object):
    """ Conditional skipping (on requirements) of tests in pytest

    Examples
    --------
    >>> @requires('numpy', 'scipy')
    ... def test_sqrt():
    ...     import numpy as np
    ...     assert np.sqrt(4) == 2
    ...     from scipy.special import zeta
    ...     assert zeta(2) < 2
    ...
    >>> @requires('numpy>=1.9.0')
    ... def test_nanmedian():
    ...     import numpy as np
    ...     a = np.array([[10.0, 7, 4], [3, 2, 1]])
    ...     a[0, 1] = np.nan
    ...     assert np.nanmedian(a) == 3
    ...

    """
    def __init__(self, *reqs):
        self.missing = []
        self.incomp = []
        self.requirements = list(parse_requirements(reqs))
        for req in self.requirements:
            try:
                mod = __import__(req.project_name)
            except ImportError:
                self.missing.append(req.project_name)
            else:
                try:
                    ver = parse_version(mod.__version__)
                except AttributeError:
                    pass
                else:
                    for rel, vstr in req.specs:
                        if not _relop[rel](ver, parse_version(vstr)):
                            self.incomp.append(str(req))

    def __call__(self, cb):
        r = 'Unfulfilled requirements.'
        if self.missing:
            r += " Missing modules: %s." % ', '.join(self.missing)
        if self.incomp:
            r += " Incomp versions: %s." % ', '.join(self.incomp)
        return pytest.mark.skipif(self.missing or self.incomp, reason=r)(cb)
