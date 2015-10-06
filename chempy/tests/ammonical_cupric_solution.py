# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function)

from chempy.chemistry import Solute
from chempy.equilibria import Equilibrium, EqSystemLog
import periodictable


def get_ammonical_cupric_eqsys():
    substances = (Hp, OHm, NH4p, NH3, H2O, Cupp, CuNH31pp, CuNH32pp, CuNH33pp,
                  CuNH34pp, CuNH35pp, Cu2OH2pp, CuOH3m, CuOH4mm, CuOH2) = [
                      Solute(n, latex_name=l, formula=periodictable.formula(n))
                      for n, l in [
                              ('H{+}', 'H^+'), ('HO{-}', 'OH^-'),
                              ('NH3 + H{+}', 'NH_4^+'),
                              ('NH3', 'NH_3'), ('H2O', 'H_2O'),
                              ('Cu{2+}', 'Cu^{2+}'),
                              ('Cu{2+}NH3', 'Cu(NH_3)^{2+}'),
                              ('Cu{2+}(NH3)2', 'Cu(NH_3)_2^{2+}'),
                              ('Cu{2+}(NH3)3', 'Cu(NH_3)_3^{2+}'),
                              ('Cu{2+}(NH3)4', 'Cu(NH_3)_4^{2+}'),
                              ('Cu{2+}(NH3)5', 'Cu(NH_3)_5^{2+}'),
                              ('2Cu{2+} + 2HO{-}', 'Cu_2(OH)_2^{2+}'),
                              ('Cu{2+} + 3HO{-}', 'Cu(OH)_3^-'),
                              ('Cu{2+} + 4HO{-}', 'Cu(OH)_4^{2-}'),
                              ('Cu{2+} + 2HO{-}', 'Cu(OH_2)(s)'),
                      ]]
    CuOH2.solid = True
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
    simpl_c0 = {k: init_conc[k] for k in simpl_subs}
    return EqSystemLog(simpl_eq, simpl_subs), simpl_c0
