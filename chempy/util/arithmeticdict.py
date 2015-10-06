from __future__ import division

from collections import defaultdict


class ArithmeticDict(defaultdict):
    """
    A dictionary which supports:

    addition, subtraction, multiplication and division

    subclassed from defaultdict, if other term/factor
    has items() the arithmetics is performed on a key
    per key basis. If AttributeError is raised when trying
    to access items() the operation is broadcasted onto all values.

    Nonexisting keys are interpreted to signal a zero

    __eq__ ignores values equal to ``self.default_factory()``
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
        try:
            for k, v in other.items():
                self[k] *= v
        except AttributeError:
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
        try:
            for k, v in other.items():
                self[k] /= v
        except AttributeError:
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
