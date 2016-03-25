# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from collections import defaultdict
import numpy as np
from ..units import (
    allclose, default_units, get_derived_unit, is_unitless, linspace,
    SI_base_registry,
    to_unitless,
    unit_of,
    unit_registry_to_human_readable,
    unit_registry_from_human_readable,
)

metre = default_units.metre
second = default_units.second
hour = default_units.hour
dm = default_units.decimetre
mole = default_units.mole
kilogram = default_units.kilogram
ampere = default_units.ampere
kelvin = default_units.kelvin
candela = default_units.candela
molar = mole / dm**3


def test_default_units():
    # In addition to the above, we guarantee that theses are available
    default_units.joule
    default_units.gray
    default_units.eV
    default_units.MeV
    default_units.metre
    default_units.decimetre
    default_units.centimetre
    default_units.micrometre
    default_units.nanometre
    default_units.gram
    default_units.molar
    default_units.hour
    default_units.perMolar_perSecond
    default_units.per100eV
    default_units.umol
    default_units.umol_per_J


def test_allclose():
    a = np.linspace(2, 3)*second
    b = np.linspace(2/3600., 3/3600.)*hour
    assert allclose(a, b)
    assert allclose([3600*second, 2*metre/hour], [1*hour, 2/3600*metre/second])


def test_is_unitless():
    assert not is_unitless(1*second)
    assert is_unitless(1)


def test_unit_of():
    assert unit_of(0.1*metre/second) == metre/second
    assert unit_of(7) == 1


def test_to_unitless():
    vals = [1.0*dm, 2.0*dm]
    result = to_unitless(vals, metre)
    assert result[0] == 0.1
    assert result[1] == 0.2

    vals = [1.0, 2.0]*dm
    result = to_unitless(vals, metre)
    assert result[0] == 0.1
    assert result[1] == 0.2

    length_unit = 1000*metre
    result = to_unitless(1.0*metre, length_unit)
    assert abs(result - 1e-3) < 1e-12

    amount_unit = 1e-9  # nano
    assert abs(to_unitless(1.0, amount_unit) - 1e9) < 1e-6

    assert abs(to_unitless(3/(second*molar),
                           metre**3/mole/second) - 3e-3) < 1e-12


def test_linspace():
    ls = linspace(2*second, 3*second)
    assert abs(to_unitless(ls[0], hour) - 2/3600.) < 1e-15


def test_get_derived_unit():
    registry = SI_base_registry.copy()
    registry['length'] = 1e-1*registry['length']
    conc_unit = get_derived_unit(registry, 'concentration')
    assert abs(conc_unit - 1*mole/(dm**3)) < 1e-12*mole/(dm**3)

    registry = defaultdict(lambda: 1)
    registry['amount'] = 1e-9  # nano
    assert abs(to_unitless(1.0, get_derived_unit(
        registry, 'concentration')) - 1e9) < 1e-6


def test_unit_registry_to_human_readable():
    # Not as much human readable as JSON serializable...
    d = defaultdict(lambda: 1)
    assert unit_registry_to_human_readable(d) == dict(
        (x, (1, 1)) for x in SI_base_registry.keys())

    ur = {
        'length': 1e3*metre,
        'mass': 1e-2*kilogram,
        'time': 1e4*second,
        'current': 1e-1*ampere,
        'temperature': 1e1*kelvin,
        'luminous_intensity': 1e-3*candela,
        'amount': 1e4*mole
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
        'length': metre,
        'mass': kilogram,
        'time': second,
        'current': ampere,
        'temperature': kelvin,
        'luminous_intensity': candela,
        'amount': mole
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
        'length': 1e3*metre,
        'mass': 1e-2*kilogram,
        'time': 1e4*second,
        'current': 1e-1*ampere,
        'temperature': 1e1*kelvin,
        'luminous_intensity': 1e-3*candela,
        'amount': 1e4*mole
    }

    assert ur != {
        'length': 1e2*metre,
        'mass': 1e-3*kilogram,
        'time': 1e2*second,
        'current': 1e-2*ampere,
        'temperature': 1e0*kelvin,
        'luminous_intensity': 1e-2*candela,
        'amount': 1e3*mole
    }
