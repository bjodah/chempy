from __future__ import division

from chempy.util.testing import requires
from ..integrated import pseudo_irrev, pseudo_rev, binary_irrev, binary_rev

import pytest

try:
    import sympy
except ImportError:
    sympy = None
else:
    one = sympy.S(1)

    t, kf, kb, prod, major, minor = sympy.symbols(
        't kf kb prod major minor', negative=False, nonnegative=True, real=True)

    subsd = {t: one*2, kf: one*3, kb: one*7, major: one*11,
             minor: one*13, prod: one*0}


@requires('sympy')
def test_pseudo_irrev():
    f = pseudo_irrev(t, kf, prod, major, minor, backend=sympy)
    dfdt = f.diff(t)
    num_dfdt = dfdt.subs(subsd)
    assert (num_dfdt - (
        major*kf*(minor - f)
    ).subs(subsd)).simplify() == 0


@requires('sympy')
def test_pseudo_rev():
    f = pseudo_rev(t, kf, kb, prod, major, minor, backend=sympy)
    dfdt = f.diff(t)
    num_dfdt = dfdt.subs(subsd)
    assert (num_dfdt - (major*kf*(minor - f) - kb*f).subs(subsd)).simplify() == 0


@pytest.mark.slow
@requires('sympy')
def test_binary_irrev():
    f = binary_irrev(t, kf, prod, major, minor, backend=sympy)
    dfdt = f.diff(t)
    num_dfdt = dfdt.subs(subsd)
    assert (num_dfdt - (kf*(minor - f)*(major - f)).subs(subsd)).simplify() == 0


@pytest.mark.slow
@requires('sympy')
def test_binary_rev():
    f = binary_rev(t, kf, kb, prod, major, minor, backend=sympy)
    dfdt = f.diff(t)
    num_dfdt = dfdt.subs(subsd)
    ans = kf*(minor - f)*(major - f) - kb*f
    # symbolic susbsitution fails:
    assert abs(float(num_dfdt) - float(ans.subs(subsd))) < 2e-14
