# -*- coding: utf-8 -*-

import sympy
from ..salcs import (
    calc_salcs_projection,
    calc_salcs_basis,
    _expand_irreducible,
    _angles_to_vectors,
    
)
from ..tables import (
    tables, 
    headers, 
    mulliken, 
    rot_trans_modes, 
    IR_active, 
    Raman_active, 
    masks, 
    atom_contribution
)

def test_calc_salcs_projection():
    a, b, c = sympy.symbols('a b c')
    
    assert calc_salcs_projection([a, b, c, a, b, c], 'c3v') == [2*a + 2*b + 2*c, 0, 2*a - b - c]
    
    
    a1, a2, e1, e2, e3 = sympy.symbols('a1, a2, e1, e2, e3')
    assert (calc_salcs_projection([e1, e2, e3, -e1, -e2, -e3, -e1, -e2, -e3, e1, e2, e3], 'd3h') == 
            [0, 0, 0, 0, 4*e1 + 4*e2 + 4*e3, 4*e1 - 2*e2 - 2*e3])
    assert (calc_salcs_projection([a1, a1, a1, -a2, -a2, -a2, -a2, -a2, -a2, a1, a1, a1], 'd3h') == 
            [6*a1 - 6*a2, 0, 0, 0, 6*a1 + 6*a2, 0]) 
    
    
    
def test_calc_salcs_basis():
    assert calc_salcs_basis([[0,0], [120,0], [240,0]], 'c3v', mode='angle') == [[1.0, 1.0, 1.0],
     0,
     [[1.0, -0.5, -0.5], 
      [0.0, 1.0, -1.0], 
      [1.0, -0.5, -0.5], 
      [0.0, -1.0, 1.0]]]
    
    assert calc_salcs_basis([[1,0,0], [0,1,0], [-1,0,0], [0,-1,0]], 'd4h', mode='vector') == [
             [1, 1, 1, 1],
              0,
             [1, -1, 1, -1],
             0, 0, 0, 0, 0, 0,
             [[1, 0, -1, 0], 
              [0, 1, 0, -1]]]