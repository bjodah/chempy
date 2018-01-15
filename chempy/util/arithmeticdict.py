# -*- coding: utf-8 -*-
from __future__ import division

from collections import defaultdict
from itertools import chain


def _imul(d1, d2):
    if hasattr(d2, 'keys'):
        for k in set(chain(d1.keys(), d2.keys())):
            d1[k] = d1[k]*d2[k]
    else:
        for k in d1:
            d1[k] *= d2


def _itruediv(d1, d2):
    if hasattr(d2, 'keys'):
        for k in set(chain(d1.keys(), d2.keys())):
            d1[k] = d1[k]/d2[k]
    else:
        for k in d1:
            d1[k] /= d2


class ArithmeticDict(defaultdict):
    """ A dictionary which supports arithmetics

    Subclassed from defaultdict, with support for addition, subtraction,
    multiplication and division. If other term/factor has a :meth:`keys` method
    the arithmetics are performed on a key per key basis. If :meth:`keys` is
    missing, the operation is broadcasted onto all values.
    Nonexisting keys are interpreted to signal a zero.

    Notes
    -----
    ``__eq__`` ignores values equal to ``self.default_factory()``

    Examples
    --------
    >>> d1 = ArithmeticDict(float, {'a': 2.0, 'b': 3.0})
    >>> d2 = ArithmeticDict(float, {'b': 5.0, 'c': 7.0})
    >>> (d1 + d2) == {'a': 2., 'b': 8., 'c': 7., 'd': 0.}
    True
    >>> (d1 * d1) == {'a': 4.0, 'b': 9.0, 'z': 0}
    True
    >>> (d1 * d2) == {'b': 15}
    True
    >>> d1*2 == {'a': 4, 'b': 6}
    True
    >>> (d1 / {'a': 2, 'b': 11})['b'] == 3./11
    True
    >>> d2/3 == {'b': 5./3, 'c': 7./3}
    True

    """

    def copy(self):
        return self.__class__(self.default_factory, self.items())

    def __iadd__(self, other):
        try:
            for k, v in other.items():
                self[k] += v
        except AttributeError:
            for k in self:
                self[k] += other
        return self

    def __isub__(self, other):
        try:
            for k, v in other.items():
                self[k] -= v
        except AttributeError:
            for k in self:
                self[k] -= other
        return self

    def __add__(self, other):
        a = self.copy()
        a += other
        return a

    def __sub__(self, other):
        a = self.copy()
        a -= other
        return a

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        return -1*self + other

    def __imul__(self, other):
        _imul(self, other)
        return self

    def __mul__(self, other):
        a = self.copy()
        a *= other
        return a

    def __rmul__(self, other):
        return self * other

    def __itruediv__(self, other):
        _itruediv(self, other)
        return self

    def __truediv__(self, other):
        a = self.copy()
        a /= other
        return a

    def __rtruediv__(self, other):
        """ other / self """
        return self.__class__(self.default_factory,
                              {k: other/v for k, v in self.items()})

    def __ifloordiv__(self, other):
        if hasattr(other, 'keys'):
            for k in set(chain(self.keys(), other.keys())):
                self[k] = self[k]//other[k]
        else:
            for k in self:
                self[k] //= other
        return self

    def __floordiv__(self, other):
        a = self.copy()
        a //= other
        return a

    def __rfloordiv__(self, other):
        """ other // self """
        return self.__class__(self.default_factory,
                              {k: other//v for k, v in self.items()})

    __div__ = __truediv__  # Py2 compatibility (or: import division from __future__)
    __idiv__ = __itruediv__  # Py2 compatibility
    __rdiv__ = __rtruediv__  # Py2 compatibility

    def __repr__(self):
        return "{}({}, {})".format(self.__class__.__name__,
                                   repr(self.default_factory),
                                   dict(self))

    def _element_eq(self, a, b):
        return a == b

    def _discrepancy(self, other, cb):
        default = self.default_factory()
        _self = self.copy()  # getitem is not idempotent on defaultdict
        _other = other.copy()
        try:
            for k in set(chain(_self.keys(), _other.keys())):
                if not cb(_self[k], _other.get(k, default)):
                    return False
            return True
        except TypeError:
            return False

    def __eq__(self, other):
        return self._discrepancy(other, self._element_eq)

    def isclose(self, other, rtol=1e-12, atol=None):
        def _isclose(a, b):
            lim = abs(rtol*b)
            if atol is not None:
                lim += atol
            return abs((a-b)) <= lim
        return self._discrepancy(other, _isclose)

    def all_non_negative(self):
        for v in self.values():
            if v < v*0:
                return False
        return True
