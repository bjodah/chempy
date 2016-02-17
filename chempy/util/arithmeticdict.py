# -*- coding: utf-8 -*-
from __future__ import division

from collections import defaultdict
from itertools import chain


class ArithmeticDict(defaultdict):
    """ A dictionary which supports arithmetics

    Subclassed from defaultdict, with support for addition, subtraction,
    multiplication and division. If other term/factor has a :meth:`keys` method
    the arithmetics are performed on a key per key basis. If :meth:`keys` is
    missing, the operation is broadcasted onto all values.
    Nonexisting keys are interpreted to signal a zero

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
        if hasattr(other, 'keys'):
            for k in set(chain(self.keys(), other.keys())):
                self[k] = self[k]*other[k]
        else:
            for k in self:
                self[k] *= other
        return self

    def __mul__(self, other):
        a = self.copy()
        a *= other
        return a

    def __rmul__(self, other):
        return self * other

    def __itruediv__(self, other):
        if hasattr(other, 'keys'):
            for k in set(chain(self.keys(), other.keys())):
                self[k] = self[k]/other[k]
        else:
            for k in self:
                self[k] /= other
        return self

    def __truediv__(self, other):
        a = self.copy()
        a /= other
        return a

    def __rtruediv__(self, other):
        """ other / self """
        return self.__class__(self.default_factory,
                              {k: other/v for k, v in self.items()})

    __div__ = __truediv__  # Py2 compability
    __idiv__ = __itruediv__  # Py2 compability
    __rdiv__ = __rtruediv__  # Py2 compability

    def __repr__(self):
        return "{}({}, {})".format(self.__class__.__name__,
                                   repr(self.default_factory),
                                   dict(self))

    def __eq__(self, other):
        default = self.default_factory()
        try:
            for k, v in self.items():
                if v == default:
                    continue
                if v != other[k]:
                    return False
            for k, v in other.items():
                if v == default:
                    continue
                if k in self:
                    continue
                return False
            return True
        except TypeError:
            return False

    def __hash__(self):
        default = self.default_factory()
        l = [self.default_factory]
        for k, v in self.items():
            if v != default:
                l.append((k, v))
        return hash(tuple(l))
