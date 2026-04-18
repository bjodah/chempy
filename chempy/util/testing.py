# -*- coding: utf-8 -*-

import os
from operator import lt, le, eq, ne, ge, gt
import pytest

_relop = dict(zip("<= == != >= > <".split(), (le, eq, ne, ge, gt, lt)))


def _parse_version(vs: str, /):
    return tuple(map(int, vs.split('.')))

class requires(object):
    """Conditional skipping (on requirements) of tests in pytest

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
        for req in reqs:
            for rs, ro in _relop.items():
                if rs in req:
                    name, version = req.split(rs)
                    version = _parse_version(version)
            else:
                name, version = req, None
            
            try:
                mod = __import__(name)
            except ImportError:
                self.missing.append(name)
            else:
                if version is not None:
                    found_version = _parse_version(mod.__version__)
                    if not (fulfilled := ro(version, found_version)):
                        self.incomp.append("%s %s %s" % (found_version, rs, version))

    def __call__(self, cb):
        r = "Unfulfilled requirements."
        if self.missing:
            r += " Missing modules: %s." % ", ".join(self.missing)
        if self.incomp:
            r += " Incomp versions: %s." % ", ".join(self.incomp)
        return skipif(self.missing or self.incomp, reason=r)(cb)


def skipif(predicate, *, reason):
    if os.environ.get("CHEMPY_SKIP_NO_TESTS", "0") == "1":
        return pytest.mark.skipif(False, reason=reason)
    else:
        return pytest.mark.skipif(predicate, reason=reason)
