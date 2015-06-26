# -*- coding: utf-8 -*-

from aqchem.equilibria import equilibrium_quotient

class Substance(object):
    """
    Parameters
    ----------
    name: str
    z: int
    mass: float
    latex_name: str
    """
    def __init__(self, name=None, z=0, mass=None, latex_name=None):
        self.name = name
        self.z = z
        self.mass = mass
        self.latex_name = latex_name


class Reaction(object):
    """
    Parameters
    ----------
    reac: dict
    prod: dict
    params: float or callable

    inact_reac: dict (optional)
    inact_prod: dict (optional)
    """

    def __init__(self, reac, prod, params=None, inact_reac=None, inact_prod=None):
        self.reac = reac
        self.prod = prod
        self.params = params
        self.inact_reac = inact_reac
        self.inact_prod = inact_prod

    def active_net_stoich(self, substances):
        return tuple(self.prod[k] - self.reac[k] for k in substances)


class Equilibrium(Reaction):

    def K(self, T=None):
        if callable(self.params):
            if T is None:
                raise ValueError("No T provided")
            return self.params(T)
        else:
            if T is not None:
                raise ValueError("T provided but params not callble")
            return self.params


class ReactionSystem(object):

    def __init__(self, rxns, substances):
        self.rxns = rxns
        self.substances = substances


class EqSystem(ReactionSystem):
    pass
