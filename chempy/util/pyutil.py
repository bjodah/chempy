# -*- coding: utf-8 -*-
"""
General utilities and exceptions.
"""
from __future__ import (absolute_import, division, print_function)

from collections import namedtuple, Mapping
import os
import warnings

from .. import __url__


def defaultnamedtuple(typename, field_names, defaults=()):
    """ Generates a new subclass of tuple with default values.

    Parameters
    ----------
    typename: string
         the name of the class
    field_names: str or iterable
        an iterable of splitable string
    defaults: iterable
        default values for field_names, counting [-len(defaults):]

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


class Deprecation(object):
    """ decorator for deprecating functions and classes

    Search the code base for more examples

    Parameters
    ----------
    last_supported_version: str
        version string, e.g. '0.2.1'
    will_be_missing_in: str
        version string, e.g. '0.3.0'
    use_instead: object or str
        function or class to use instead or descriptive string
    issue: str
        issue identifier, e.g. 'gh-15' where 15 is the issue number on github.

    Examples
    --------
    >>> @Deprecation(last_supported_version='0.4.0')
    ... def some_old_function(x):
    ...     return x*x - x
    ...


    Notes
    -----
    :class:`DeprecationWarning` is ignored by default. Use custom warning
    and filter appropriately. Alternatively, run python with ``-W`` flag or set
    the appropriate environment variable:

    ::

        $ python -c 'import warnings as w; w.warn("X", DeprecationWarning)'
        $ python -Wd -c 'import warnings as w; w.warn("X", DeprecationWarning)'
        -c:1: DeprecationWarning: X
        $ export PYTHONWARNINGS=d
        $ python -c 'import warnings as w; w.warn("X", DeprecationWarning)'
        -c:1: DeprecationWarning: X

    """

    def __init__(self, last_supported_version=None, will_be_missing_in=None,
                 use_instead=None, issue=None,
                 deprecation_warning=DeprecationWarning):
        self.last_supported_version = last_supported_version
        self.will_be_missing_in = will_be_missing_in
        self.use_instead = use_instead
        self.issue = issue
        self.deprecation_warning = deprecation_warning
        self.warning_message = self._warning_message_template()

    def _warning_message_template(self):
        msg = '%(func_name)s is deprecated'
        if self.last_supported_version is not None:
            msg += ' since (not including) % s' % self.last_supported_version
        if self.will_be_missing_in is not None:
            msg += ', it will be missing in %s' % self.will_be_missing_in
        if self.issue is not None:
            assert self.issue.startswith('gh-')  # currently tracked at github
            msg += ' (see %s)' % (__url__ + '/issues/' + self.issue[3:])
        if self.use_instead is not None:
            try:
                msg += '. Use %s instead' % self.use_instead.__name__
            except AttributeError:
                msg += '. Use %s instead' % self.use_instead
        return msg + '.'

    def __call__(self, wrapped):
        """ Decorates function to be deprecated """
        msg = self.warning_message % {'func_name': wrapped.__name__}
        wrapped_doc = wrapped.__doc__ or ''
        if hasattr(wrapped, '__mro__'):  # wrapped is a class
            class _Wrapper(wrapped):
                __doc__ = msg + '\n\n' + wrapped_doc

                def __init__(_self, *args, **kwargs):
                    warnings.warn(msg, self.deprecation_warning, stacklevel=3)
                    wrapped.__init__(_self, *args, **kwargs)

        else:  # wrapped is a function
            def _Wrapper(*args, **kwargs):
                warnings.warn(msg, self.deprecation_warning, stacklevel=3)
                return wrapped(*args, **kwargs)
            _Wrapper.__doc__ = msg + '\n\n' + wrapped_doc

        _Wrapper._deprecation = self
        _Wrapper.__name__ = wrapped.__name__
        _Wrapper.__module__ = wrapped.__module__
        return _Wrapper


def deprecated(*args, **kwargs):
    """ Helper to :class:`Deprecation` for using ChemPyDeprecationWarning. """
    return Deprecation(*args, deprecation_warning=ChemPyDeprecationWarning,
                       **kwargs)

warnings.simplefilter(os.environ.get('CHEMPY_DEPRECATION_FILTER', 'once'),
                      ChemPyDeprecationWarning)
