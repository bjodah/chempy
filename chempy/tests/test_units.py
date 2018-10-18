# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from collections import defaultdict

try:
    import numpy as np
except ImportError:
    np = None
import pytest

from ..util.testing import requires
from ..units import (
    allclose, concatenate, get_derived_unit, is_unitless, linspace, logspace_from_lin,
    SI_base_registry, unitless_in_registry, format_string, get_physical_dimensionality,
    to_unitless, magnitude, default_unit_in_registry, Backend, latex_of_unit,
    unit_of, unit_registry_to_human_readable, units_library, simplified, uniform,
    unit_registry_from_human_readable, _sum, UncertainQuantity, compare_equality,
    default_units as u, patched_numpy as pnp, default_constants as dc
)


@requires(units_library)
def test_default_units():
    u.metre
    u.second
    u.hour
    u.decimetre
    u.mole
    u.kilogram
    u.ampere
    u.kelvin
    u.candela
    u.molar
    u.per100eV
    u.joule
    u.gray
    u.eV
    u.MeV
    u.metre
    u.decimetre
    u.centimetre
    u.micrometre
    u.nanometre
    u.gram
    u.molar
    u.hour
    u.perMolar_perSecond
    u.per100eV
    u.umol
    u.umol_per_J


@requires(units_library)
def test_allclose():
    assert allclose(42, 42)
    assert allclose(42*u.meter, 0.042*u.km)
    assert not allclose(42, 43)
    assert not allclose(42, 42*u.meter)
    assert not allclose(42, 43*u.meter)
    assert not allclose(42*u.meter, 42)

    a = np.linspace(2, 3)*u.second
    b = np.linspace(2/3600., 3/3600.)*u.hour
    assert allclose(a, b)
    assert allclose([3600*u.second, 2*u.metre/u.hour],
                    [1*u.hour, 2/3600*u.metre/u.second])
    c1 = [[3000, 4000], [3000, 4000]]*u.mol/u.metre**3
    c2 = [[3000, 4000], [436.2, 5281.89]]*u.mol/u.metre**3
    assert not allclose(c1, c2)
    assert allclose(0*u.second, 0*u.second)

    # Possibly allow comparison with scalars in future (broadcasting):
    # assert allclose(2, [2, 2])
    # assert allclose([2, 2], 2)

    # assert not allclose(2, [2, 3])
    # assert not allclose([2, 3], 2)

    # assert allclose(2*u.second, [2, 2]*u.second)
    # assert allclose([2, 2]*u.second, 2*u.second)

    # assert not allclose(2*u.second, [2, 3]*u.second)
    # assert not allclose([2, 3]*u.second, 2*u.second)


@requires(units_library)
def test_is_unitless():
    assert not is_unitless(1*u.second)
    assert is_unitless(1)
    assert is_unitless({'a': 1, 'b': 2.0})
    assert not is_unitless({'a': 2, 'b': 5.0*u.second, 'c': 3})
    assert is_unitless(7*u.molar/u.mole*u.dm3)
    assert is_unitless([2, 3, 4])
    assert not is_unitless([2*u.m, 3*u.m])
    assert not is_unitless([3, 4*u.m])


@requires(units_library)
def test_unit_of():
    assert compare_equality(unit_of(0.1*u.metre/u.second), u.metre/u.second)
    assert not compare_equality(unit_of(0.1*u.metre/u.second), u.kilometre/u.second)
    assert compare_equality(unit_of(7), 1)
    assert unit_of(u.gray).dimensionality == u.gray.dimensionality
    ref = (u.joule/u.kg).simplified.dimensionality
    assert unit_of(u.gray, simplified=True).dimensionality == ref

    assert compare_equality(unit_of(dict(foo=3*u.molar, bar=2*u.molar)), u.molar)
    assert not compare_equality(unit_of(dict(foo=3*u.molar, bar=2*u.molar)), u.second)
    with pytest.raises(Exception):
        unit_of(dict(foo=3*u.molar, bar=2*u.second))
    assert not compare_equality(unit_of(dict(foo=3*u.molar, bar=2*u.molar)), u.mol/u.metre**3)


@requires(units_library)
def test_to_unitless():
    dm = u.decimetre
    vals = [1.0*dm, 2.0*dm]
    result = to_unitless(vals, u.metre)
    assert result[0] == 0.1
    assert result[1] == 0.2
    with pytest.raises(ValueError):
        to_unitless([42, 43], u.metre)

    with pytest.raises(ValueError):
        to_unitless(np.array([42, 43]), u.metre)

    vals = [1.0, 2.0]*dm
    result = to_unitless(vals, u.metre)
    assert result[0] == 0.1
    assert result[1] == 0.2

    length_unit = 1000*u.metre
    result = to_unitless(1.0*u.metre, length_unit)
    assert abs(result - 1e-3) < 1e-12

    amount_unit = 1e-9  # nano
    assert abs(to_unitless(1.0, amount_unit) - 1e9) < 1e-6
    assert abs(to_unitless(3/(u.second*u.molar),
                           u.metre**3/u.mole/u.second) - 3e-3) < 1e-12
    assert abs(to_unitless(2*u.dm3, u.cm3) - 2000) < 1e-12
    assert abs(to_unitless(2*u.m3, u.dm3) - 2000) < 1e-12
    assert (float(to_unitless(UncertainQuantity(2, u.dm3, .3), u.cm3)) - 2000) < 1e-12

    g1 = UncertainQuantity(4.46, u.per100eV, 0)
    g_unit = get_derived_unit(SI_base_registry, 'radiolytic_yield')
    assert abs(to_unitless(g1, g_unit) - 4.46 * 1.036e-7) < 1e-9
    g2 = UncertainQuantity(-4.46, u.per100eV, 0)
    assert abs(to_unitless(-g2, g_unit) - 4.46 * 1.036e-7) < 1e-9

    vals = np.array([1.*dm, 2.*dm], dtype=object)
    result = to_unitless(vals, u.metre)
    assert result[0] == 0.1
    assert result[1] == 0.2

    one_billionth_molar_in_nanomolar = to_unitless(1e-9*u.molar, u.nanomolar)
    assert one_billionth_molar_in_nanomolar == 1


@requires(units_library)
def test_UncertainQuantity():
    a = UncertainQuantity([1, 2], u.m, [.1, .2])
    assert a[1] == [2.]*u.m
    assert (-a)[0] == [-1.]*u.m
    assert (-a).uncertainty[0] == [0.1]*u.m
    assert (-a)[0] == (a*-1)[0]
    assert (-a).uncertainty[0] == (a*-1).uncertainty[0]


@requires(units_library, 'sympy')
def test_to_unitless__sympy():
    import sympy as sp
    assert sp.cos(to_unitless(sp.pi)) == -1
    with pytest.raises(AttributeError):
        to_unitless(sp.pi, u.second)


@requires(units_library)
def test_linspace():
    ls = linspace(2*u.second, 3*u.second)
    assert abs(to_unitless(ls[0], u.hour) - 2/3600.) < 1e-15


@requires(units_library)
def test_logspace_from_lin():
    ls = logspace_from_lin(2*u.second, 3*u.second)
    assert abs(to_unitless(ls[0], u.hour) - 2/3600.) < 1e-15
    assert abs(to_unitless(ls[-1], u.hour) - 3/3600.) < 1e-15


@requires(units_library)
def test_get_derived_unit():
    registry = SI_base_registry.copy()
    registry['length'] = 1e-1*registry['length']
    conc_unit = get_derived_unit(registry, 'concentration')
    dm = u.decimetre
    assert abs(conc_unit - 1*u.mole/(dm**3)) < 1e-12*u.mole/(dm**3)

    registry = defaultdict(lambda: 1)
    registry['amount'] = 1e-9  # nano
    assert abs(to_unitless(1.0, get_derived_unit(
        registry, 'concentration')) - 1e9) < 1e-6


@requires(units_library)
def test_unit_registry_to_human_readable():
    # Not as much human readable as JSON serializable...
    d = defaultdict(lambda: 1)
    assert unit_registry_to_human_readable(d) == dict(
        (x, (1, 1)) for x in SI_base_registry.keys())

    ur = {
        'length': 1e3*u.metre,
        'mass': 1e-2*u.kilogram,
        'time': 1e4*u.second,
        'current': 1e-1*u.ampere,
        'temperature': 1e1*u.kelvin,
        'luminous_intensity': 1e-3*u.candela,
        'amount': 1e4*u.mole
    }
    assert unit_registry_to_human_readable(ur) == {
        'length': (1e3, 'm'),
        'mass': (1e-2, 'kg'),
        'time': (1e4, 's'),
        'current': (1e-1, 'A'),
        'temperature': (1e1, 'K'),
        'luminous_intensity': (1e-3, 'cd'),
        'amount': (1e4, 'mol')
    }
    assert unit_registry_to_human_readable(ur) != {
        'length': (1e2, 'm'),
        'mass': (1e-2, 'kg'),
        'time': (1e4, 's'),
        'current': (1e-1, 'A'),
        'temperature': (1e1, 'K'),
        'luminous_intensity': (1e-3, 'cd'),
        'amount': (1e4, 'mol')
    }


@requires(units_library)
def test_unit_registry_from_human_readable():
    hr = unit_registry_to_human_readable(defaultdict(lambda: 1))
    assert hr == dict((x, (1, 1)) for x in SI_base_registry.keys())
    ur = unit_registry_from_human_readable(hr)
    assert ur == dict((x, 1) for x in SI_base_registry.keys())

    hr = unit_registry_to_human_readable(SI_base_registry)
    assert hr == {
        'length': (1.0, 'm'),
        'mass': (1.0, 'kg'),
        'time': (1.0, 's'),
        'current': (1.0, 'A'),
        'temperature': (1.0, 'K'),
        'luminous_intensity': (1.0, 'cd'),
        'amount': (1.0, 'mol')
    }
    ur = unit_registry_from_human_readable(hr)
    assert ur == SI_base_registry

    ur = unit_registry_from_human_readable({
        'length': (1.0, 'm'),
        'mass': (1.0, 'kg'),
        'time': (1.0, 's'),
        'current': (1.0, 'A'),
        'temperature': (1.0, 'K'),
        'luminous_intensity': (1.0, 'cd'),
        'amount': (1.0, 'mol')
    })
    assert ur == {
        'length': u.metre,
        'mass': u.kilogram,
        'time': u.second,
        'current': u.ampere,
        'temperature': u.kelvin,
        'luminous_intensity': u.candela,
        'amount': u.mole
    }

    ur = unit_registry_from_human_readable({
        'length': (1e3, 'm'),
        'mass': (1e-2, 'kg'),
        'time': (1e4, 's'),
        'current': (1e-1, 'A'),
        'temperature': (1e1, 'K'),
        'luminous_intensity': (1e-3, 'cd'),
        'amount': (1e4, 'mol')
    })
    assert ur == {
        'length': 1e3*u.metre,
        'mass': 1e-2*u.kilogram,
        'time': 1e4*u.second,
        'current': 1e-1*u.ampere,
        'temperature': 1e1*u.kelvin,
        'luminous_intensity': 1e-3*u.candela,
        'amount': 1e4*u.mole
    }

    assert ur != {
        'length': 1e2*u.metre,
        'mass': 1e-3*u.kilogram,
        'time': 1e2*u.second,
        'current': 1e-2*u.ampere,
        'temperature': 1e0*u.kelvin,
        'luminous_intensity': 1e-2*u.candela,
        'amount': 1e3*u.mole
    }


@requires(units_library)
def test_unitless_in_registry():
    mag = magnitude(unitless_in_registry(3*u.per100eV, SI_base_registry))
    ref = 3*1.0364268834527753e-07
    assert abs(mag - ref) < 1e-14
    ul = unitless_in_registry([3*u.per100eV, 5*u.mol/u.J], SI_base_registry)
    assert allclose(ul, [ref, 5], rtol=1e-6)


@requires(units_library)
def test_compare_equality():
    assert compare_equality(3*u.m, 3*u.m)
    assert compare_equality(3*u.m, 3e-3*u.km)
    assert compare_equality(3e+3*u.mm, 3*u.m)
    assert not compare_equality(3*u.m, 2*u.m)
    assert not compare_equality(3*u.m, 3*u.s)
    assert not compare_equality(3*u.m, 3*u.m**2)
    assert not compare_equality(3*u.m, np.array(3))
    assert not compare_equality(np.array(3), 3*u.m)
    assert compare_equality([3, None], [3, None])
    assert not compare_equality([3, None, 3], [3, None, None])
    assert not compare_equality([None, None, 3], [None, None, 2])
    assert compare_equality([3*u.m, None], [3, None])
    assert not compare_equality([3*u.m, None], [3*u.km, None])


@requires(units_library)
def test_get_physical_dimensionality():
    assert get_physical_dimensionality(3*u.mole) == {'amount': 1}
    assert get_physical_dimensionality([3*u.mole]) == {'amount': 1}
    assert get_physical_dimensionality(42) == {}


@requires(units_library)
def test_default_unit_in_registry():
    mol_per_m3 = default_unit_in_registry(3*u.molar, SI_base_registry)
    assert magnitude(mol_per_m3) == 1
    assert mol_per_m3 == u.mole/u.metre**3

    assert default_unit_in_registry(3, SI_base_registry) == 1
    assert default_unit_in_registry(3.0, SI_base_registry) == 1


@requires(units_library)
def test__sum():
    # sum() does not work here...
    assert (_sum([0.1*u.metre, 1*u.decimetre]) - 2*u.decimetre)/u.metre == 0


@requires(units_library)
def test_Backend():
    b = Backend()
    with pytest.raises(ValueError):
        b.exp(-3*u.metre)
    assert abs(b.exp(1234*u.metre/u.kilometre) - b.exp(1.234)) < 1e-14


@requires(units_library, 'numpy')
def test_Backend__numpy():
    import numpy as np
    b = Backend(np)
    b.sum([1000*u.metre/u.kilometre, 1], axis=0) == 2.0

    with pytest.raises(AttributeError):
        b.Piecewise


@requires('sympy')
def test_Backend__sympy():
    b = Backend('sympy')
    b.sin(b.pi) == 0

    with pytest.raises(AttributeError):
        b.min


@requires(units_library)
def test_format_string():
    assert format_string(3*u.gram/u.metre**2) == ('3', 'g/m**2')
    assert format_string(3*u.gram/u.metre**2, tex=True) == ('3', r'\mathrm{\frac{g}{m^{2}}}')


@requires(units_library)
def test_joule_html():
    joule_htm = 'kg&sdot;m<sup>2</sup>/s<sup>2</sup>'
    joule = u.J.dimensionality.simplified
    assert joule.html == joule_htm


@requires(units_library)
def test_latex_of_unit():
    assert latex_of_unit(u.gram/u.metre**2) == r'\mathrm{\frac{g}{m^{2}}}'


@requires(units_library)
def test_concatenate():
    a = [1, 2]*u.metre
    b = [2, 3]*u.mm
    ref = [1, 2, 2e-3, 3e-3]*u.metre
    assert allclose(concatenate((a, b)), ref)


@requires(units_library)
def test_pow0():
    a = [1, 2]*u.metre
    b = a**0
    assert allclose(b, [1, 1])

    c = a**2
    assert allclose(c, [1, 4]*u.m**2)


@requires(units_library)
def test_tile():
    a = [2*u.m, 3*u.km]
    assert allclose(pnp.tile(a, 2), [2*u.m, 3000*u.m, 2e-3*u.km, 3*u.km])


@requires(units_library)
def test_simplified():
    assert allclose(simplified(dc.molar_gas_constant), 8.314*u.J/u.mol/u.K, rtol=2e-3)


@requires(units_library)
def test_polyfit_polyval():
    p1 = pnp.polyfit([0, 1, 2], [0, 1, 4], 2)
    assert allclose(p1, [1, 0, 0], atol=1e-14)
    assert allclose(pnp.polyval(p1, 3), 9)
    assert allclose(pnp.polyval(p1, [4, 5]), [16, 25])

    p2 = pnp.polyfit([0, 1, 2]*u.s, [0, 1, 4]*u.m, 2)
    for _p, _r, _a in zip(p2, [1*u.m/u.s**2, 0*u.m/u.s, 0*u.m],
                          [0*u.m/u.s**2, 1e-15*u.m/u.s, 1e-15*u.m]):
        assert allclose(_p, _r, atol=_a)
    assert allclose(pnp.polyval(p2, 3*u.s), 9*u.m)
    assert allclose(pnp.polyval(p2, [4, 5]*u.s), [16, 25]*u.m)


@requires(units_library)
def test_uniform():
    base = [3*u.km, 200*u.m]
    refs = [np.array([3000, 200]), np.array([3, 0.2])]

    def _check(case, ref):
        assert np.any(np.all(magnitude(uniform(case)) == ref, axis=1))
    _check(base, refs)
    _check(tuple(base), refs)
    keys = 'foo bar'.split()
    assert magnitude(uniform(dict(zip(keys, base)))) in [dict(zip(keys, r)) for r in refs]
