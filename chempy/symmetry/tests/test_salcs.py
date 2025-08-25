# -*- coding: utf-8 -*-

import sympy
from ..salcs import (
    calc_salcs_projection,
    calc_salcs_func,
    _expand_irreducible,
    _angles_to_vectors,

)


def test_calc_salcs_projection():
    # ammonia hydrogens
    a, b, c = sympy.symbols('a b c')
    assert (calc_salcs_projection([a, b, c, a, b, c], 'c3v') ==
            [2*a + 2*b + 2*c, 0, 2*a - b - c])

    # trigonal bipyramidal
    a1, a2, e1, e2, e3 = sympy.symbols('a1, a2, e1, e2, e3')
    assert (calc_salcs_projection([e1, e2, e3, -e1, -e2, -e3, -e1,
                                   -e2, -e3, e1, e2, e3], 'd3h') ==
            [0, 0, 0, 0, 4*e1 + 4*e2 + 4*e3, 4*e1 - 2*e2 - 2*e3])
    assert (calc_salcs_projection([a1, a1, a1, -a2, -a2, -a2, -a2,
                                   -a2, -a2, a1, a1, a1], 'd3h') ==
            [6*a1 - 6*a2, 0, 0, 0, 6*a1 + 6*a2, 0])


def test_calc_salcs_func():
    # square planar
    a, b, c, d = sympy.symbols('a b c d')
    salc_true1 = [a + b + c + d, 0, a - b + c - d, 0, 0, 0, 0, 0, 0,
                 [a - c, b - d]]
    assert (calc_salcs_func([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
                            'd4h', [a, b, c, d], mode='vector') == salc_true1)

    # trigonal bipyramidal
    a1, a2, e1, e2, e3 = sympy.symbols('a1, a2, e1, e2, e3')
    salc_true2 = [[1.0*e1 + 1.0*e2 + 1.0*e3, 1.0*a1 + 1.0*a2], 0,
                 [1.0*e1 - 0.5*e2 - 0.5*e3, 1.0*e2 - 1.0*e3,
                  1.0*e1 - 0.5*e2 - 0.5*e3, 1.0*e2 - 1.0*e3], 0,
                 1.0*a1 - 1.0*a2, 0]
    angles = [[0, 0], [120, 0], [240, 0], [0, 90], [0, -90]]
    assert(calc_salcs_func(angles, 'd3h', [e1, e2, e3, a1, a2], mode='angle')
           == salc_true2)

    # seesaw - such as SF4
    salc_true3 =  [a1 + a2, 0, a1 - a2, 0]
    a1, a2, e1, e2 = sympy.symbols('a1 a2 e1 e2')
    assert(calc_salcs_func([[0, 0], [-180, 0], [90, -30], [-90, -30]],
                           'c2v', [a1, a2, e1, e2], mode='angle') == salc_true3)

    # octahedral
    a, b, c, d, e, f = sympy.symbols('a b c d e f')
    salc_true4 = [a + b + c + d + e + f, 0,
                  [-a - b - c - d + 2*e + 2*f, a - b + c - d],
                  0, 0, 0, 0, 0, [a - c, b - d, e - f], 0]
    oh_angle = calc_salcs_func([[0, 0], [90, 0], [180, 0], [270, 0], [0, 90],
                                [0, -90]], 'oh', [a, b, c, d, e, f], mode='angle')
    oh_vector = calc_salcs_func([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0],
                                 [0, 0, 1], [0, 0, -1]], 'oh', [a, b, c, d, e, f])
    assert(oh_angle == salc_true4)
    assert(oh_vector == salc_true4)
