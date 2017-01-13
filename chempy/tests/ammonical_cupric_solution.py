# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function)

from chempy.chemistry import Species, Equilibrium
from chempy.equilibria import EqSystem


def get_ammonical_cupric_eqsys():
    NH3_complexes = ['CuNH3+2', 'Cu(NH3)2+2', 'Cu(NH3)3+2', 'Cu(NH3)4+2',
                     'Cu(NH3)5+2']
    OH_complexes = ['Cu2(OH)2+2', 'Cu(OH)3-', 'Cu(OH)4-2']
    substances = [
        Species.from_formula(n) for n in [
            'H+', 'OH-', 'NH4+', 'NH3', 'H2O', 'Cu+2'] +
        NH3_complexes + OH_complexes + ['Cu(OH)2(s)']
    ]

    (Hp, OHm, NH4p, NH3, H2O, Cupp, CuNH31pp, CuNH32pp,
     CuNH33pp, CuNH34pp, CuNH35pp, Cu2OH2pp, CuOH3m,
     CuOH4mm, CuOH2) = [s.name for s in substances]
    substances[-1].phase_idx = 1
    init_conc = {Hp: 1e-7, OHm: 1e-7, NH4p: 0, NH3: 1.0, Cupp: 1e-2,
                 CuNH31pp: 0, CuNH32pp: 0, CuNH33pp: 0, CuNH34pp: 0,
                 CuNH35pp: 0, H2O: 55.5, Cu2OH2pp: 0, CuOH2: 0, CuOH3m: 0,
                 CuOH4mm: 0}
    H2O_c = init_conc[H2O]
    w_autop = Equilibrium({H2O: 1}, {Hp: 1, OHm: 1}, 10**-14/H2O_c)
    NH4p_pr = Equilibrium({NH4p: 1}, {Hp: 1, NH3: 1}, 10**-9.26)
    CuOH2_s = Equilibrium({CuOH2: 1}, {Cupp: 1, OHm: 2}, 10**-18.8)
    CuOH_B3 = Equilibrium({CuOH2: 1, OHm: 1}, {CuOH3m: 1}, 10**-3.6)
    CuOH_B4 = Equilibrium({CuOH2: 1, OHm: 2}, {CuOH4mm: 1}, 10**-2.7)
    Cu2OH2 = Equilibrium({Cupp: 2, H2O: 2}, {Cu2OH2pp: 1, Hp: 2},
                         10**-10.6 / H2O_c**2)
    CuNH3_B1 = Equilibrium({CuNH31pp: 1}, {Cupp: 1, NH3: 1}, 10**-4.3)
    CuNH3_B2 = Equilibrium({CuNH32pp: 1}, {Cupp: 1, NH3: 2}, 10**-7.9)
    CuNH3_B3 = Equilibrium({CuNH33pp: 1}, {Cupp: 1, NH3: 3}, 10**-10.8)
    CuNH3_B4 = Equilibrium({CuNH34pp: 1}, {Cupp: 1, NH3: 4}, 10**-13.0)
    CuNH3_B5 = Equilibrium({CuNH35pp: 1}, {Cupp: 1, NH3: 5}, 10**-12.4)
    equilibria = (w_autop, NH4p_pr, CuNH3_B1, CuNH3_B2, CuNH3_B3, CuNH3_B4,
                  CuNH3_B5, Cu2OH2, CuOH_B3, CuOH_B4, CuOH2_s)

    new_eqs = CuOH2_s - CuOH_B3, CuOH2_s - CuOH_B4
    skip_subs, skip_eq = (1, 3)
    simpl_subs = substances[:-skip_subs]
    simpl_eq = equilibria[:-skip_eq] + new_eqs
    simpl_c0 = {k.name: init_conc[k.name] for k in substances[:-skip_subs]}
    return EqSystem(simpl_eq, simpl_subs), simpl_c0
