from __future__ import division

from chempy.util.testing import requires
from ..integrated import pseudo_irrev, pseudo_rev, binary_irrev, binary_rev

try:
    import sympy
except ImportError:
    sympy = None
else:
    one = sympy.S(1)

    t, kf, P0, t0, excess_C, limiting_C, eps_l, beta = sympy.symbols(
        't k_f P0 t0 Y Z epsilon beta', negative=False)  # t0 => -t0

    subsd = {t: one*2, kf: one*3, P0: one*5, t0: one*7, excess_C: one*11,
             limiting_C: one*13, eps_l: one*17, beta: one*23}
    excl_params = {t0: one*0, P0: one, eps_l: one}


@requires('sympy')
def test_pseudo_irrev():
    f = pseudo_irrev(t, kf, P0, -t0, excess_C, limiting_C, eps_l,
                     backend=sympy).subs(excl_params)
    dfdt = f.diff(t)
    num_dfdt = dfdt.subs(subsd)
    assert (num_dfdt - (
        excess_C*kf*(limiting_C - f)
    ).subs(subsd)).simplify() == 0


@requires('sympy')
def test_pseudo_rev():
    f = pseudo_rev(t, kf, P0, -t0, excess_C, limiting_C, eps_l,
                   beta, backend=sympy).subs(excl_params)
    dfdt = f.diff(t)
    num_dfdt = dfdt.subs(subsd)
    assert (num_dfdt - (
        excess_C*kf*(limiting_C - f) - kf/beta*f
    ).subs(subsd)).simplify() == 0


@requires('sympy')
def test_binary_irrev():
    f = binary_irrev(t, kf, P0, -t0, excess_C, limiting_C, eps_l,
                     backend=sympy).subs(excl_params)
    dfdt = f.diff(t)
    num_dfdt = dfdt.subs(subsd)
    assert (num_dfdt - (
        kf*(limiting_C - f)*(excess_C - f)
    ).subs(subsd)).simplify() == 0


@requires('sympy')
def test_binary_rev():
    f = binary_rev(t, kf, P0, -t0, excess_C, limiting_C, eps_l,
                   beta, backend=sympy).subs(excl_params)
    dfdt = f.diff(t)
    num_dfdt = dfdt.subs(subsd)
    ans = kf*(limiting_C - f)*(excess_C - f) - kf/beta*f
    # symbolic susbsitution fails:
    assert abs(float(num_dfdt) - float(ans.subs(subsd))) < 2e-14
