# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np

from operator import itemgetter
from collections import defaultdict


def elements(formula):
    """
    Returns a dict mapping {periodictable.core.Element: int}

    Parameters
    ----------
    formula: periodictable.formulas.Formula

    """
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
    composition: dict or None (default)
        dict (int -> int) e.g. {atomic number: count}
    excess_electrons: bool (default: True)
        when this is set to True and composition is None:
        composition[0] will be number of excess electrons
        (-1 for H+ for example).
    """
    def __init__(self, name=None, charge=None, mass=None,
                 latex_name=None, formula=None,
                 composition=None, excess_electrons=True):
        self.name = name
        self.latex_name = latex_name

        if isinstance(formula, str):
            import periodictable
            formula = periodictable.formula(formula)
        self.formula = formula

        self.charge = charge
        if charge is None:
            try:
                self.charge = formula.charge
            except AttributeError:
                pass

        self.mass = mass
        if mass is None:
            try:
                self.mass = formula.mass
            except AttributeError:
                pass

        self.composition = composition
        if composition is None and formula is not None:
            self.composition = {
                k.number: v for k, v in elements(formula).items()}
            if excess_electrons:
                self.composition[0] = -self.charge

    def __repr__(self):
        kw = ['name=' + self.name + ', ...']  # Too verbose to print all
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

    str_arrow = '->'
    latex_arrow = '\\rightarrow'

    def __init__(self, reac, prod, params=None, inact_reac=None,
                 inact_prod=None):
        self.reac = reac
        self.prod = prod
        self.params = params
        self.inact_reac = inact_reac
        self.inact_prod = inact_prod

    def __repr__(self):
        try:
            s = ' K=%.2g' % self.params
        except:
            s = ''
        return self._get_str('name', 'str_arrow') + s

    def __eq__(lhs, rhs):
        if not isinstance(lhs, Reaction) or not isinstance(rhs, Reaction):
            return NotImplemented
        for attr in ['reac', 'prod', 'params', 'inact_reac', 'inact_prod']:
            if getattr(lhs, attr) != getattr(rhs, attr):
                return False
        return True

    def net_stoich(self, substances):
        return tuple(self.prod.get(k, 0) - self.reac.get(k, 0)
                     for k in substances)

    def _get_str(self, name_attr, arrow_attr):
        reac, prod = [[
            ((str(v)+' ') if v > 1 else '') + getattr(k, name_attr)
            for k, v in filter(itemgetter(1), d.items())
        ] for d in (self.reac, self.prod)]
        fmtstr = "{} " + getattr(self, arrow_attr) + " {}"
        return fmtstr.format(" + ".join(reac),
                             " + ".join(prod))

    def latex(self):
        return self._get_str('latex_name', 'latex_arrow')


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

    def as_per_substance_array(self, cont):
        """ Turns e.g. a dict into an ordered array """
        if isinstance(cont, np.ndarray):
            pass
        elif isinstance(cont, dict):
            cont = [cont[k] for k in self.substances]
        cont = np.asarray(cont)
        if cont.shape[-1] != self.ns:
            raise ValueError("Incorrect size")
        return cont

    def as_substance_index(self, sbstnc):
        """ Returns the index of a Substance in the system"""
        if isinstance(sbstnc, int):
            return sbstnc
        else:
            return self.substances.index(sbstnc)

    @property
    def stoichs(self):
        return np.array([eq.net_stoich(self.substances)
                         for eq in self.rxns])
