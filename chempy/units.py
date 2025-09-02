# -*- coding: utf-8 -*-
""" The units module provides the following attributes:

- ``chempy.units.default_units``
- ``chempy.units.default_constants``
- ``chempy.units.SI_base_registry``

together with some functions.

Currently `quantities <https://pypi.python.org/pypi/quantities>`_ is used as
the underlying package to handle units. If it is possible you should try to
only use the ``chempy.units`` module (since it is likely that ``ChemPy``
will change this backend at some point in the future). Therefore you should not
rely on any attributes of the ``Quantity`` instances (and rather use
getter & setter functions in `chempy.units`).

"""

from functools import reduce
from operator import mul
import sys
import warnings
from typing import Iterable

from .util.arithmeticdict import ArithmeticDict
from .util.pyutil import NameSpace, deprecated

units_library = "quantities"  # info used for selective testing.


try:
    pq = __import__(units_library)
except ImportError:
    UncertainQuantity = None
    default_constants = None
    default_units = None
    SI_base_registry = None
    np = None
else:
    import numpy as np
    from .util._quantities import _patch_quantities

    _patch_quantities(pq)
    UncertainQuantity = pq.UncertainQuantity
    # Let us extend the underlying pq namespace with some common units in
    # chemistry
    default_constants = NameSpace(pq.constants)

    default_units = NameSpace(pq)
    default_units.dm = default_units.decimetre = pq.UnitQuantity(
        "decimetre", default_units.m / 10.0, u_symbol="dm"
    )
    default_units.m3 = default_units.metre ** 3
    default_units.dm3 = default_units.decimetre ** 3
    default_units.cm3 = default_units.centimetre ** 3
    if not hasattr(default_units, "molar"):
        default_units.molar = pq.UnitQuantity(
            "M", 1e3 * default_units.mole / default_units.m3, u_symbol="M"
        )
    if not hasattr(default_units, "millimolar"):
        default_units.millimolar = pq.UnitQuantity(
            "mM", 1 * default_units.mole / default_units.m3, u_symbol="mM"
        )
    if not hasattr(default_units, "micromolar"):
        default_units.micromolar = pq.UnitQuantity(
            "uM", 1e-3 * default_units.mole / default_units.m3, u_symbol="μM"
        )
    if not hasattr(default_units, "nanomolar"):
        default_units.nanomolar = pq.UnitQuantity(
            "nM", 1e-6 * default_units.mole / default_units.m3, u_symbol="nM"
        )
    if not hasattr(default_units, "molal"):
        default_units.molal = pq.UnitQuantity(
            "molal", default_units.mole / default_units.kg, u_symbol="molal"
        )
    if not hasattr(default_units, "per100eV"):
        default_units.per100eV = pq.UnitQuantity(
            "per100eV",
            1 / (100 * default_units.eV * default_constants.Avogadro_constant),
            u_symbol="(100eV)**-1",
        )
    if not hasattr(default_units, "micromole"):
        default_units.micromole = pq.UnitQuantity(
            "micromole", pq.mole / 1e6, u_symbol="μmol"
        )
    if not hasattr(default_units, "nanomole"):
        default_units.nanomole = pq.UnitQuantity(
            "nanomole", pq.mole / 1e9, u_symbol="nmol"
        )
    if not hasattr(default_units, "kilojoule"):
        default_units.kilojoule = pq.UnitQuantity(
            "kilojoule", 1e3 * pq.joule, u_symbol="kJ"
        )
    if not hasattr(default_units, "kilogray"):
        default_units.kilogray = pq.UnitQuantity(
            "kilogray", 1e3 * pq.gray, u_symbol="kGy"
        )
    if not hasattr(default_units, "perMolar_perSecond"):
        default_units.perMolar_perSecond = 1 / default_units.molar / pq.s
    if not hasattr(default_units, "umol"):
        default_units.umol = default_units.micromole
    if not hasattr(default_units, "umol_per_J"):
        default_units.umol_per_J = default_units.umol / pq.joule

    # unit registry data and logic:

    SI_base_registry = {
        "length": default_units.metre,
        "mass": default_units.kilogram,
        "time": default_units.second,
        "current": default_units.ampere,
        "temperature": default_units.kelvin,
        "luminous_intensity": default_units.candela,
        "amount": default_units.mole,
    }


def magnitude(value):
    try:
        return value.magnitude
    except AttributeError:
        return value


def uncertainty(uval):
    uncert = uval.uncertainty
    if not is_quantity(uncert):
        warnings.warn(f"Handling unexpected type: {type(uval)}")
    return uncert


def simplified(value):
    if hasattr(value, "simplified"):
        return value.simplified
    else:
        return to_unitless(value)


def is_quantity(arg):
    if arg.__class__.__name__ == "Quantity":
        return True  # this checks works even if quantities is not installed.
    else:
        return False


# SI Base Quantities:
time = ArithmeticDict(int, {"time": 1})
length = ArithmeticDict(int, {"length": 1})
mass = ArithmeticDict(int, {"mass": 1})
current = ArithmeticDict(int, {"current": 1})
temperature = ArithmeticDict(int, {"temperature": 1})
amount = ArithmeticDict(int, {"amount": 1})
# intensity = ArithmeticDict(int, {'intensity': 1}) what's wrong with photon flux? (human eyes, bah!)

energy = ArithmeticDict(int, {"mass": 1, "length": 2, "time": -2})
volume = ArithmeticDict(int, {"length": 3})
concentration = {"amount": 1} - volume


def get_derived_unit(registry, key):
    """Get the unit of a physical quantity in a provided unit system.

    Parameters
    ----------
    registry: dict (str: unit)
        mapping 'length', 'mass', 'time', 'current', 'temperature',
        'luminous_intensity', 'amount'. If registry is ``None`` the
        function returns 1.0 unconditionally.
    key: str
        one of the registry keys or one of: 'diffusivity', 'electricalmobility',
        'permittivity', 'charge', 'energy', 'concentration', 'density',
        'radiolytic_yield'.

    Examples
    --------
    >>> m, s = default_units.meter, default_units.second
    >>> bool(get_derived_unit(SI_base_registry, 'diffusivity') == m**2/s)
    True

    """
    if registry is None:
        return 1.0
    derived = {
        "diffusivity": registry["length"] ** 2 / registry["time"],
        "electrical_mobility": (
            registry["current"] * registry["time"] ** 2 / registry["mass"]
        ),
        "permittivity": (
            registry["current"] ** 2
            * registry["time"] ** 4
            / (registry["length"] ** 3 * registry["mass"])
        ),
        "charge": registry["current"] * registry["time"],
        "energy": registry["mass"] * registry["length"] ** 2 / registry["time"] ** 2,
        "concentration": registry["amount"] / registry["length"] ** 3,
        "density": registry["mass"] / registry["length"] ** 3,
    }
    derived["diffusion"] = derived["diffusivity"]  # 'diffusion' is deprecated
    derived["radiolytic_yield"] = registry["amount"] / derived["energy"]
    derived["doserate"] = derived["energy"] / registry["mass"] / registry["time"]
    derived["linear_energy_transfer"] = derived["energy"] / registry["length"]

    try:
        return derived[key]
    except KeyError:
        return registry[key]


def unit_registry_to_human_readable(unit_registry):
    """Serialization of a unit registry."""
    if unit_registry is None:
        return None
    new_registry = {}
    integer_one = 1
    for k in SI_base_registry:
        if unit_registry[k] is integer_one:
            new_registry[k] = 1, 1
        else:
            dim_list = list(unit_registry[k].dimensionality)
            if len(dim_list) != 1:
                raise TypeError("Compound units not allowed: {}".format(dim_list))
            u_symbol = dim_list[0].u_symbol
            new_registry[k] = float(unit_registry[k]), u_symbol
    return new_registry


def _latex_from_dimensionality(dim):
    # see https://github.com/python-quantities/python-quantities/issues/148
    from quantities.markup import format_units_latex

    return format_units_latex(dim, mult=r"\\cdot")


def latex_of_unit(quant):
    """Returns LaTeX representation of the unit of a quantity

    Examples
    --------
    >>> print(latex_of_unit(1/default_units.kelvin))
    \\mathrm{\\frac{1}{K}}

    """
    return _latex_from_dimensionality(quant.dimensionality).strip("$")


def unicode_of_unit(quant):
    """Returns unicode representation of the unit of a quantity

    Examples
    --------
    >>> print(unicode_of_unit(1/default_units.kelvin))
    1/K

    """
    return quant.dimensionality.unicode


def html_of_unit(quant):
    """Returns HTML representation of the unit of a quantity

    Examples
    --------
    >>> print(html_of_unit(2*default_units.m**2))
    m<sup>2</sup>

    """
    return quant.dimensionality.html


def unit_registry_from_human_readable(unit_registry):
    """Deserialization of unit_registry."""
    if unit_registry is None:
        return None
    new_registry = {}
    for k in SI_base_registry:
        factor, u_symbol = unit_registry[k]
        if u_symbol == 1:
            unit_quants = [1]
        else:
            unit_quants = list(pq.Quantity(0, u_symbol).dimensionality.keys())

        if len(unit_quants) != 1:
            raise TypeError("Unknown UnitQuantity: {}".format(unit_registry[k]))
        else:
            new_registry[k] = factor * unit_quants[0]
    return new_registry


# Abstraction of underlying package providing units and dimensional analysis:


def is_unitless(expr):
    """Returns ``True`` if ``expr`` is unitless, otherwise ``False``

    Examples
    --------
    >>> is_unitless(42)
    True
    >>> is_unitless(42*default_units.kilogram)
    False

    """
    if hasattr(expr, "dimensionality"):
        if expr.dimensionality == pq.dimensionless:
            return True
        else:
            return expr.simplified.dimensionality == pq.dimensionless.dimensionality
    if isinstance(expr, dict):
        return all(is_unitless(_) for _ in expr.values())
    elif isinstance(expr, (tuple, list)):
        return all(is_unitless(_) for _ in expr)
    return True


def unit_of(expr, simplified=False):
    """Returns the unit of a quantity

    Examples
    --------
    >>> bool(unit_of(42*pq.second) == unit_of(12*pq.second))
    True
    >>> unit_of(42)
    1

    """
    if isinstance(expr, (tuple, list)):
        return unit_of(uniform(expr)[0], simplified)
    elif isinstance(expr, dict):
        return unit_of(list(uniform(expr).values())[0], simplified)

    try:
        if simplified:
            return expr.units.simplified
        else:
            return expr.units
    except AttributeError:
        return 1


def rescale(value, unit):
    try:
        return value.rescale(unit)
    except AttributeError:
        if unit == 1:
            return value
        else:
            raise


def to_unitless(value, new_unit=None):
    """ Nondimensionalization of a quantity.

    Parameters
    ----------
    value: quantity
    new_unit: unit

    Examples
    --------
    >>> '%.1g' % to_unitless(1*default_units.metre, default_units.nm)
    '1e+09'
    >>> '%.1g %.1g' % tuple(to_unitless([1*default_units.m, 1*default_units.mm], default_units.nm))
    '1e+09 1e+06'

    """
    integer_one = 1
    if new_unit is None:
        new_unit = pq.dimensionless

    if isinstance(value, np.ndarray) and not hasattr(value, 'rescale'):
        if is_unitless(new_unit) and new_unit == 1:
            if value.ndim == 0:
                return value.item()
            else:
                return value
        return np.array([to_unitless(elem, new_unit) for elem in value])
    elif isinstance(value, dict):
        new_value = dict(value.items())  # value.copy()
        for k in value:
            new_value[k] = to_unitless(value[k], new_unit)
        return new_value
    elif isinstance(value, (int, float)) and new_unit is integer_one or new_unit is None:
        return value
    elif isinstance(value, Iterable) and getattr(value, 'ndim', -1) != 0:
        if len(value) == 0:
            return to_unitless(value[0], new_unit)
        else:
            return np.array([to_unitless(elem, new_unit) for elem in value])
    elif isinstance(value, str):
        raise ValueError("str not supported")

    ori_unit = unit_of(value)
    ori_mag = magnitude(value)
    scale_factor = simplified(ori_unit/new_unit)
    if is_unitless(scale_factor):
        sc = magnitude(scale_factor)
        if getattr(ori_mag, 'ndim', -1) == 0 or not isinstance(ori_mag, Iterable):
            return ori_mag*sc
        else:
            return np.array(ori_mag*sc)
    else:
        raise ValueError("Incompatible units")


def uniform(container):
    """Turns a list, tuple or dict with mixed units into one with uniform units.

    Parameters
    ----------
    container : tuple, list or dict

    Examples
    --------
    >>> km, m = default_units.kilometre, default_units.metre
    >>> uniform(dict(a=3*km, b=200*m))  # doctest: +SKIP
    {'b': array(200.0) * m, 'a': array(3000.0) * m}

    """
    if isinstance(container, (tuple, list)):
        unit = unit_of(container[0])
    elif isinstance(container, dict):
        unit = unit_of(list(container.values())[0])
        return container.__class__(
            [(k, to_unitless(v, unit) * unit) for k, v in container.items()]
        )
    else:
        return container
    return to_unitless(container, unit) * unit


def get_physical_dimensionality(value):
    if is_unitless(value):
        return {}
    _quantities_mapping = {
        pq.UnitLength: "length",
        pq.UnitMass: "mass",
        pq.UnitTime: "time",
        pq.UnitCurrent: "current",
        pq.UnitTemperature: "temperature",
        pq.UnitLuminousIntensity: "luminous_intensity",
        pq.UnitSubstance: "amount",
    }
    return {
        _quantities_mapping[k.__class__]: v
        for k, v in uniform(value).simplified.dimensionality.items()
    }


@deprecated(use_instead=get_physical_dimensionality, will_be_missing_in="0.8.0")
def get_physical_quantity(value):
    return get_physical_dimensionality(value)


def _get_unit_from_registry(dimensionality, registry):
    return reduce(mul, [registry[k] ** v for k, v in dimensionality.items()])


def default_unit_in_registry(value, registry):
    _dimensionality = get_physical_dimensionality(value)
    if _dimensionality == {}:
        return 1
    return _get_unit_from_registry(_dimensionality, registry)


def unitless_in_registry(value, registry):
    _default_unit = default_unit_in_registry(value, registry)
    return to_unitless(value, _default_unit)


# NumPy like functions for compatibility:


def compare_equality(a, b):
    """Returns True if two arguments are equal.

    Both arguments need to have the same dimensionality.

    Parameters
    ----------
    a : quantity
    b : quantity

    Examples
    --------
    >>> km, m = default_units.kilometre, default_units.metre
    >>> compare_equality(3*km, 3)
    False
    >>> bool(compare_equality(3*km, 3000*m))
    True

    """
    # Work around for https://github.com/python-quantities/python-quantities/issues/146
    try:
        a + b
    except TypeError:
        # We might be dealing with e.g. None (None + None raises TypeError)
        try:
            len(a)
        except TypeError:
            # Assumed scalar
            return a == b
        else:
            if len(a) != len(b):
                return False
            return all(compare_equality(_a, _b) for _a, _b in zip(a, b))
    except ValueError:
        return False
    else:
        return a == b


def allclose(a, b, rtol=1e-8, atol=None):
    """Analogous to ``numpy.allclose``."""
    if a.__class__.__name__ == "UncertainQuantity":
        return allclose(pq.Quantity(a), b, rtol=rtol, atol=atol)
    if b.__class__.__name__ == "UncertainQuantity":
        return allclose(a, pq.Quantity(b), rtol=rtol, atol=atol)

    try:
        d = abs(a - b)
    except Exception:
        try:
            if len(a) == len(b):
                return all(allclose(_a, _b, rtol, atol) for _a, _b in zip(a, b))
            else:
                return False
        except Exception:
            return False
    lim = abs(a) * rtol
    if atol is not None:
        lim += atol

    try:
        len(d)
    except TypeError:
        return d <= lim
    else:
        try:
            len(lim)
        except TypeError:
            return np.all([_d <= lim for _d in d])
        else:
            return np.all([_d <= _lim for _d, _lim in zip(d, lim)])


def linspace(start, stop, num=50):
    """Analogous to ``numpy.linspace``.

    Examples
    --------
    >>> bool(abs(linspace(2, 8, num=3)[1] - 5) < 1e-15)
    True

    """

    # work around for quantities v0.10.1 and NumPy
    unit = unit_of(start)
    start_ = to_unitless(start, unit)
    stop_ = to_unitless(stop, unit)
    return np.linspace(start_, stop_, num) * unit


def logspace_from_lin(start, stop, num=50):
    """Logarithmically spaced data points

    Examples
    --------
    >>> bool(abs(logspace_from_lin(2, 8, num=3)[1] - 4) < 1e-15)
    True

    """
    unit = unit_of(start)
    start_ = np.log2(to_unitless(start, unit))
    stop_ = np.log2(to_unitless(stop, unit))
    return np.exp2(np.linspace(start_, stop_, num)) * unit


def _sum(iterable):
    try:
        result = next(iterable)
    except TypeError:
        result = iterable[0]
        for elem in iterable[1:]:
            result += elem
        return result
    else:
        try:
            while True:
                result += next(iterable)
        except StopIteration:
            return result
        else:
            raise ValueError("Not sure how this point was reached")


class Backend(object):
    """Wrapper around modules such as numpy and math

    Instances of Backend wraps a module, e.g. `numpy` and ensures that
    arguments passed on are unitless, i.e. it raises an error if a
    transcendental function is used with quantities with units.

    Parameters
    ----------
    underlying_backend : module, str or tuple of str
        e.g. 'numpy' or ('sympy', 'math')

    Examples
    --------
    >>> import math
    >>> km, m = default_units.kilometre, default_units.metre
    >>> math.exp(3*km) == math.exp(3*m)
    True
    >>> be = Backend('math')
    >>> be.exp(3*km)  # DOCTEST: +ELLIPSIS
    Traceback (most recent call last):
        ...
    ValueError: Incompatible units
    >>> import numpy as np
    >>> np.sum([1000*pq.metre/pq.kilometre, 1]).item()
    1001.0
    >>> be_np = Backend(np)
    >>> be_np.sum([[1000*pq.metre/pq.kilometre, 1], [3, 4]], axis=1)
    array([2., 7.])

    """

    def __init__(self, underlying_backend=("numpy", "math")):
        if isinstance(underlying_backend, tuple):
            for name in underlying_backend:
                try:
                    self.be = __import__(name)
                except ImportError:
                    continue
                else:
                    break
            else:
                raise ValueError("Could not import any of %s" % str(underlying_backend))
        elif isinstance(underlying_backend, str):
            self.be = __import__(underlying_backend)
        else:
            self.be = underlying_backend

    def __getattr__(self, attr):
        be_attr = getattr(self.be, attr)
        if callable(be_attr):
            return lambda *args, **kwargs: be_attr(*map(to_unitless, args), **kwargs)
        else:
            return be_attr


# TODO: decide whether to deprecate in favor of "number_to_scientific_latex"?
def format_string(value, precision="%.5g", tex=False):
    """Formats a scalar with unit as two strings

    Parameters
    ----------
    value: float with unit
    precision: str
    tex: bool
       LaTeX formatted or not? (no '$' signs)

    Examples
    --------
    >>> print(' '.join(format_string(0.42*default_units.mol/default_units.decimetre**3)))
    0.42 mol/decimetre**3
    >>> print(' '.join(format_string(2/default_units.s, tex=True)))
    2 \\mathrm{\\frac{1}{s}}

    """
    if tex:
        unit_str = latex_of_unit(value)
    else:
        from quantities.markup import config

        attr = "unicode" if config.use_unicode else "string"
        unit_str = getattr(value.dimensionality, attr)
    return precision % float(value.magnitude), unit_str


def concatenate(arrays, **kwargs):
    """Patched version of numpy.concatenate

    Examples
    --------
    >>> from chempy.units import default_units as u
    >>> all(concatenate(([2, 3]*u.s, [4, 5]*u.s)) == [2, 3, 4, 5]*u.s)
    True

    """
    unit = unit_of(arrays[0])
    result = np.concatenate([to_unitless(arr, unit) for arr in arrays], **kwargs)
    return result * unit


def tile(array, *args, **kwargs):
    """Patched version of numpy.tile (with support for units)"""
    try:
        elem = array[0, ...]
    except TypeError:
        elem = array[0]

    unit = unit_of(elem)
    result = np.tile(to_unitless(array, unit), *args, **kwargs)
    return result * unit


def polyfit(x, y, deg, **kwargs):
    u_x = unit_of(x[0])
    u_y = unit_of(y[0])
    _x, _y = to_unitless(x, u_x), to_unitless(y, u_y)
    p = np.polyfit(_x, _y, deg)
    return [v * u_y * u_x ** (i - deg) for i, v in enumerate(p)]


def polyval(p, x):
    try:
        u_x = unit_of(x[0])
    except (TypeError, IndexError):
        u_x = unit_of(x)
    u_y = unit_of(p[-1])
    deg = len(p) - 1
    _p = [to_unitless(v, u_y * u_x ** (i - deg)) for i, v in enumerate(p)]
    _x = to_unitless(x, u_x)
    _y = np.polyval(_p, _x)
    return _y * u_y


def _wrap_numpy(k):
    numpy_func = getattr(np, k)
    if sys.version_info[0] > 2:
        from functools import wraps
    else:

        def wraps(_meta_fun):
            return lambda x: x  # py2: numpy.ufunc lacks "__module__"

    @wraps(numpy_func)
    def f(*args, **kwargs):
        return numpy_func(*map(to_unitless, args), **kwargs)

    return f


if np is None:
    patched_numpy = None
else:
    patched_numpy = NameSpace(np)
    patched_numpy.allclose = allclose
    patched_numpy.concatenate = concatenate
    patched_numpy.linspace = linspace
    patched_numpy.tile = tile
    patched_numpy.polyfit = polyfit
    patched_numpy.polyval = polyval
    for k in "log log10 log2 log1p exp expm1 logaddexp logaddexp2".split():
        setattr(patched_numpy, k, _wrap_numpy(k))


def fold_constants(arg):
    if hasattr(arg, "dimensionality"):
        m = arg.magnitude
        d = 1
        for k, v in arg.dimensionality.items():
            if isinstance(k, pq.UnitConstant):
                m = m * k.simplified ** v
            else:
                d = d * k ** v
        return m * d
    else:
        return arg
