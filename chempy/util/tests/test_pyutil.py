# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from ..pyutil import defaultkeydict, defaultnamedtuple


def test_defaultnamedtuple():
    Point2 = defaultnamedtuple('Point2', 'x y', [10])
    p = Point2(3)
    assert p.x == 3 and p.y == 10

    Point3 = defaultnamedtuple('Point3', 'x y z', [10, 20])
    p = Point3(3)
    assert p.x == 3 and p.y == 10 and p.z == 20

    p = Point3(3, z=30)
    assert p.x == 3 and p.y == 10 and p.z == 30

    p = Point3(3, 4, 5)
    assert p.x == 3 and p.y == 4 and p.z == 5

    _default_y = Point2.__new__.__defaults__[-1]

    class MySubclass(Point2):

        def __new__(cls, x2, y2=_default_y):
            return super(MySubclass, cls).__new__(cls, x2**.5, y2**.5)

    p2 = MySubclass(9, 4)
    assert isinstance(p2, tuple)
    assert isinstance(p2, Point2)
    assert isinstance(p2, MySubclass)
    assert not isinstance(p, MySubclass)
    assert p2.x == 3
    assert p2.y == 2

    p3 = MySubclass(9)
    assert p3.x == 3
    assert p3.y == 10**.5


def test_defaultkeydict():
    d = defaultkeydict(lambda k: k*2)
    assert d['as'] == 'asas'
