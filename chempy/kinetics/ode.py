from __future__ import (absolute_import, division, print_function)


def law_of_mass_action_rates(y, rsys):
    for rxn in rsys.rxns:
        rate = rxn.params
        for substance, coeff in rxn.reac.items():
            idx = rsys.as_substance_index(substance)
            rate *= y[idx]**coeff
        yield rate


def dCdt(rsys, rates):
    f = [0]*rsys.ns
    for idx_s, sbstnc in enumerate(rsys.substances):
        for idx_r, rxn in enumerate(rsys.rxns):
            f[idx_s] += rsys.stoichs[idx_r, idx_s]*rates[idx_r]
    return f


def get_odesys(rsys, SymbolicSys=None):
    if SymbolicSys is None:
        from pyodesys.symbolic import SymbolicSys
    y = SymbolicSys.symarray('y', rsys.ns)
    t = SymbolicSys.Symbol('t')
    rates = list(law_of_mass_action_rates(y, rsys))
    return SymbolicSys(zip(y, dCdt(rsys, rates)), t,
                       names=[s.name for s in rsys.substances])


def plot():
    """ Convenience function for ploting results from OdeSys.integrate() """
    pass
