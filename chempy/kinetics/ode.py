from __future__ import (absolute_import, division, print_function)


def law_of_mass_action_rates(y, rsys, k=None):
    for idx, rxn in enumerate(rsys.rxns):
        if k is None:
            rate = rxn.param
        else:
            rate = k[idx]
        for substance, coeff in rxn.reac.items():
            idx = rsys.as_substance_index(substance)
            rate *= y[idx]**coeff
        yield rate


def dCdt(rsys, rates):
    f = [0]*rsys.ns
    net_stoichs = rsys.net_stoichs()
    for idx_s, sbstnc in enumerate(rsys.substances):
        for idx_r, rxn in enumerate(rsys.rxns):
            f[idx_s] += net_stoichs[idx_r, idx_s]*rates[idx_r]
    return f


def get_odesys(rsys, include_params=False, SymbolicSys=None):
    if SymbolicSys is None:
        from pyodesys.symbolic import SymbolicSys

    def dydt(t, y, p):
        rates = list(law_of_mass_action_rates(
            y, rsys, rsys.params() if include_params else p))
        return dCdt(rsys, rates)

    return SymbolicSys.from_callback(dydt, rsys.ns,
                                     0 if include_params else rsys.nr)
