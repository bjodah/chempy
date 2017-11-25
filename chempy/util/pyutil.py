# -*- coding: utf-8 -*-
"""
General utilities and exceptions.
"""
from __future__ import (absolute_import, division, print_function)

from collections import defaultdict, namedtuple, Mapping
from functools import wraps
import os
import warnings

from .. import __url__
from .deprecation import Deprecation


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


class NoConvergence(Exception):
    pass  # why is this not in the Python standard library!?


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
