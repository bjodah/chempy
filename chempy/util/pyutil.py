# -*- coding: utf-8 -*-
"""
General utilities and exceptions.
"""
from __future__ import (absolute_import, division, print_function)

from collections import namedtuple, Mapping
import os
import warnings

from .. import __url__
from .deprecation import Deprecation


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
    return Deprecation(*args, issues_url=__url__ + '/issues/%s/',
                       warning=ChemPyDeprecationWarning, **kwargs)

warnings.simplefilter(os.environ.get('CHEMPY_DEPRECATION_FILTER', 'once'),
                      ChemPyDeprecationWarning)
