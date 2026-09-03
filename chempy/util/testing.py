# -*- coding: utf-8 -*-

import os
from operator import lt, le, eq, ne, ge, gt
import re
import pytest

_relop = dict(zip("<= == != >= > <".split(), (le, eq, ne, ge, gt, lt)))


def _parse_version(vs):
    parts = []
    for part in re.split(r'[.\-+_]', vs):
        match = re.match(r'(\d+)', part)
        if match:
            parts.append(int(match.group(1)))
        elif part:
            break
    return tuple(parts)


def _parse_requirement(req):
    for rel in sorted(_relop, key=len, reverse=True):
        if rel in req:
            name, version = req.split(rel, 1)
            return name.strip(), rel, _parse_version(version.strip())
    return req.strip(), None, None


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
            name, rel, version = _parse_requirement(req)
            try:
                mod = __import__(name)
            except ImportError:
                self.missing.append(name)
            else:
                if version is not None:
                    try:
                        found_version = _parse_version(mod.__version__)
                    except AttributeError:
                        pass
                    else:
                        if not _relop[rel](found_version, version):
                            self.incomp.append(req)

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
