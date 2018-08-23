# -*- coding: utf-8 -*-
"""
General utilities and exceptions.
"""
from __future__ import (absolute_import, division, print_function)

from collections import defaultdict, namedtuple, Mapping, OrderedDict
from functools import wraps
from itertools import product
import os
import warnings

from .. import __url__
from .deprecation import Deprecation


class NoConvergence(Exception):
    pass


class ChemPyDeprecationWarning(DeprecationWarning):
    pass


def deprecated(*args, **kwargs):
    """ Helper to :class:`Deprecation` for using ChemPyDeprecationWarning. """
    return Deprecation(
        *args, issues_url=lambda s: __url__ + '/issues/' +
        s.lstrip('gh-'), warning=ChemPyDeprecationWarning, **kwargs)

warnings.simplefilter(os.environ.get('CHEMPY_DEPRECATION_FILTER', 'once'),
                      ChemPyDeprecationWarning)


class DeferredImport(object):
    def __init__(self, modname, arg=None, decorators=None):
        self._modname = modname
        self._arg = arg
        self._decorators = decorators
        self._cache = None

    @property
    def cache(self):
        if self._cache is None:
            if self._arg is None:
                obj = __import__(self._modname)
            else:
                obj = getattr(__import__(self._modname, globals(), locals(), [self._arg]), self._arg)
            if self._decorators is not None:
                for deco in self._decorators:
                    obj = deco(obj)
            self._cache = obj
        return self._cache

    def __getattribute__(self, attr):
        if attr in ('_modname', '_arg', '_cache', 'cache', '_decorators'):
            return object.__getattribute__(self, attr)
        else:
            return getattr(self.cache, attr)

    def __call__(self, *args, **kwargs):
        return self.cache(*args, **kwargs)


class NameSpace:
    """ Used to wrap, e.g. modules.

    Parameters
    ----------
    default : module
        The underlying module. Acts as a fallback for attribute access.

    Examples
    --------
    >>> import numpy
    >>> my_numpy = NameSpace(numpy)
    >>> my_numpy.array = lambda *args, **kwargs: list(numpy.array(*args, **kwargs))
    >>> isinstance(my_numpy.array([2, 3]), list)
    True
    >>> isinstance(numpy.array([2, 3]), list)
    False

    """

    def __init__(self, default):
        self._NameSpace_default = default
        self._NameSpace_attr_store = {}

    def __getattr__(self, attr):
        if attr.startswith('_NameSpace_'):
            return self.__dict__[attr]
        else:
            try:
                return self._NameSpace_attr_store[attr]
            except KeyError:
                return getattr(self._NameSpace_default, attr)

    def __setattr__(self, attr, val):
        if attr.startswith('_NameSpace_'):
            self.__dict__[attr] = val
        else:
            self._NameSpace_attr_store[attr] = val

    def as_dict(self):
        items = self._NameSpace_default.__dict__.items()
        result = {k: v for k, v in items if not k.startswith('_')}
        result.update(self._NameSpace_attr_store)
        return result


class AttributeContainer(object):
    """ Used to turn e.g. a dictionary to a module-like object.

    Parameters
    ----------
    d : dict

    Examples
    --------
    >>> def RT(T, const):
    ...     return T*const.molar_gas_constant
    ...
    >>> from quantities import constants
    >>> RT(273.15, constants)
    array(273.15) * R
    >>> my_constants = AttributeContainer(molar_gas_constant=42)
    >>> RT(273.15, my_constants)
    11472.3

    """
    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

    def as_dict(self):
        return self.__dict__.copy()


class defaultkeydict(defaultdict):
    """ defaultdict where default_factory should have the signature key -> value

    Examples
    --------
    >>> d = defaultkeydict(lambda k: '[%s]' % k, {'a': '[a]', 'b': '[B]'})
    >>> d['a']
    '[a]'
    >>> d['b']
    '[B]'
    >>> d['c']
    '[c]'

    """

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError("Missing key: %s" % key)
        else:
            self[key] = self.default_factory(key)
        return self[key]


def defaultnamedtuple(typename, field_names, defaults=()):
    """ Generates a new subclass of tuple with default values.

    Parameters
    ----------
    typename : string
        The name of the class.
    field_names : str or iterable
        An iterable of splitable string.
    defaults : iterable
        Default values for ``field_names``, counting ``[-len(defaults):]``.

    Examples
    --------
    >>> Body = defaultnamedtuple('Body', 'x y z density', (1.0,))
    >>> Body.__doc__
    'Body(x, y, z, density)'
    >>> b = Body(10, z=3, y=5)
    >>> b._asdict()
    OrderedDict([('x', 10), ('y', 5), ('z', 3), ('density', 1.0)])

    Returns
    -------
    A new tuple subclass named ``typename``

    """
    Tuple = namedtuple(typename, field_names)
    Tuple.__new__.__defaults__ = (None,) * len(Tuple._fields)
    if isinstance(defaults, Mapping):
        Tuple.__new__.__defaults__ = tuple(Tuple(**defaults))
    else:
        nmissing = len(Tuple._fields) - len(defaults)
        defaults = (None,)*nmissing + tuple(defaults)
        Tuple.__new__.__defaults__ = tuple(Tuple(*defaults))
    return Tuple


def multi_indexed_cases(od):
    """ Returns a list of length-2 tuples

    Each tuple consist of a multi-index (tuple of integers) and a dictionary.

    Parameters
    ----------
    od : OrderedDict
        Maps each key to a number of values. ``list`` and ``tuple`` instances
        are converted to ``OrderedDict``.

    Examples
    --------
    >>> from chempy.util.pyutil import multi_indexed_cases
    >>> cases = multi_indexed_cases([('a', [1, 2, 3]), ('b', [False, True])])
    >>> len(cases)
    6
    >>> midxs, case_kws = zip(*cases)
    >>> midxs[0]
    (0, 0)
    >>> case_kws[0] == {'a': 1, 'b': False}
    True

    Returns
    -------
    List of length-2 tuples, each consisting of one tuple of inidices and one dictionary.

    """
    if isinstance(od, (list, tuple)):
        od = OrderedDict(od)
    if not isinstance(od, OrderedDict):
        raise NotImplementedError("Expected an OrderedDcit")
    keys, values = map(list, [od.keys(), od.values()])
    return [(mi, {k: v[i] for k, v, i in zip(keys, values, mi)}) for mi in product(*map(range, map(len, values)))]


def memoize(max_nargs=0):
    def decorator(func):
        @wraps(func)
        def wrapper(*args):
            if max_nargs is not None and len(args) > max_nargs:
                raise ValueError("memoization error")
            if args not in wrapper.results:
                wrapper.results[args] = func(*args)
            return wrapper.results[args]
        wrapper.results = {}
        return wrapper
    return decorator
