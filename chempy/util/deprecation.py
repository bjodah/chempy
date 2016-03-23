# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import warnings


class Deprecation(object):
    """ factory of decorators for deprecating functions or classes

    Search the code base for more examples. Note that the wrapper gets
    :attr:`_deprecation` attribute with the :class:`Deprecation` instance.

    Parameters
    ----------
    last_supported_version : str
        Version string, e.g. ``'0.2.1'``.
    will_be_missing_in : str
        Version string, e.g. ``'0.3.0'``.
    use_instead : object or str
        Function or class to use instead or descriptive string.
    issue : str
        Issue identifier, e.g. 'gh-15' where 15 is the issue number on github.
        Note that 'gh-' prefix is a special case and stripped away.
    issues_url : str
        Url temlpate, e.g. ``'https://github.com/user/repo/issues/%s/'``.
    warning: DeprecationWarning
        Any subclass of DeprecationWarning, tip: you may invoke:
        ``warnings.simplefilter('once', MyWarning)`` at module init.

    Examples
    --------
    >>> import warnings
    >>> warnings.simplefilter("error", DeprecationWarning)
    >>> @Deprecation()
    ... def f():
    ...     return 1
    ...
    >>> f()  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    DeprecationWarning: f is deprecated.
    >>> @Deprecation(last_supported_version='0.4.0')
    ... def some_old_function(x):
    ...     return x*x - x
    ...
    >>> some_old_function._deprecation.last_supported_version
    '0.4.0'
    >>> @Deprecation(will_be_missing_in='1.0')
    ... class ClumsyClass(object):
    ...     pass
    ...
    >>> ClumsyClass._deprecation.will_be_missing_in
    '1.0'
    >>> warnings.resetwarnings()

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
                 use_instead=None, issue=None, issues_url='%s',
                 warning=DeprecationWarning):
        self.last_supported_version = last_supported_version
        self.will_be_missing_in = will_be_missing_in
        self.use_instead = use_instead
        self.issue = issue
        self.issues_url = issues_url
        self.warning = warning
        self.warning_message = self._warning_message_template()

    def _warning_message_template(self):
        msg = '%(func_name)s is deprecated'
        if self.last_supported_version is not None:
            msg += ' since (not including) % s' % self.last_supported_version
        if self.will_be_missing_in is not None:
            msg += ', it will be missing in %s' % self.will_be_missing_in
        if self.issue is not None:
            if self.issue.startswith('gh-'):
                msg += ' (see %s)' % (self.issues_url % self.issue[3:])
            else:
                msg += ' (see %s)' % (self.issues_url % self.issue)
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
                    warnings.warn(msg, self.warning, stacklevel=3)
                    wrapped.__init__(_self, *args, **kwargs)

        else:  # wrapped is a function
            def _Wrapper(*args, **kwargs):
                warnings.warn(msg, self.warning, stacklevel=3)
                return wrapped(*args, **kwargs)
            _Wrapper.__doc__ = msg + '\n\n' + wrapped_doc

        _Wrapper._deprecation = self
        _Wrapper.__name__ = wrapped.__name__
        _Wrapper.__module__ = wrapped.__module__
        return _Wrapper
