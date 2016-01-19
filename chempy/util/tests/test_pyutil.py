# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from ..pyutil import defaultnamedtuple


def test_defaultnamedtuple():
    Point2 = defaultnamedtuple('Point2', 'x y', [10])
    p = Point2(3)
    assert p.x == 3 and p.y == 10

    Point3 = defaultnamedtuple('Point2', 'x y z', [10, 20])
    p = Point3(3)
    assert p.x == 3 and p.y == 10 and p.z == 20

    p = Point3(3, z=30)
    assert p.x == 3 and p.y == 10 and p.z == 30

    p = Point3(3, 4, 5)
    assert p.x == 3 and p.y == 4 and p.z == 5
