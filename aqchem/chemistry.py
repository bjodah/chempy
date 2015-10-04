# -*- coding: utf-8 -*-

from __future__ import division

from operator import itemgetter
from collections import defaultdict


def elements(formula):
    # see https://github.com/pkienzle/periodictable/issues/14
    d = defaultdict(int)
    for atm, n in formula.atoms.items():
        try:
            d[atm.element] += n
        except AttributeError:
            d[atm] += n
    return d


class Substance(object):
    """
    Parameters
    ----------
    name: str
    charge: int
    mass: float
    latex_name: str
    formula: carries dict attribute "atoms"
    composition: dict (int -> int) e.g. {atomic number: count}
    """
    def __init__(self, name=None, charge=None, mass=None,
                 latex_name=None, formula=None,
                 composition=None):
        self.name = name
        self.charge = charge
        self.mass = mass
        self.latex_name = latex_name
        self.formula = formula
        self.composition = composition

        if charge is None:
            try:
                self.charge = formula.charge
            except AttributeError:
                pass
        if mass is None:
            try:
                self.mass = formula.mass
            except AttributeError:
                pass
        if composition is None and formula is not None:
            self.composition = {
                k.number: v for k, v in elements(formula).items()}

    def __repr__(self):
        kw = ['name=' + self.name + ', ...']  # Too verbose
        return "{}({})".format(self.__class__.__name__, ','.join(kw))

    def __str__(self):
        return str(self.name)


class Solute(Substance):

    def __init__(self, *args, **kwargs):
        self.solid = kwargs.pop('solid', False)
        super(self.__class__, self).__init__(*args, **kwargs)


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

    latex_arrow = '\\rightarrow'

    def __init__(self, reac, prod, params=None, inact_reac=None,
                 inact_prod=None):
        self.reac = reac
        self.prod = prod
        self.params = params
        self.inact_reac = inact_reac
        self.inact_prod = inact_prod

    def net_stoich(self, substances):
        return tuple(self.prod.get(k, 0) - self.reac.get(k, 0)
                     for k in substances)

    def latex(self):
        reac, prod = [[
            ((str(v)+' ') if v > 1 else '') + getattr(k, 'latex_name')
            for k, v in filter(itemgetter(1), d.items())
        ] for d in (self.reac, self.prod)]
        fmtstr = "{} " + self.latex_arrow + " {}"
        return fmtstr.format(" + ".join(reac),
                             " + ".join(prod))


class ReactionSystem(object):

    def __init__(self, rxns, substances):
        self.rxns = rxns
        self.substances = substances

    @property
    def nr(self):
        return len(self.rxns)

    @property
    def ns(self):
        return len(self.substances)
