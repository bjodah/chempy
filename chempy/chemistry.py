# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np

from itertools import chain
from operator import itemgetter
from collections import defaultdict

from .util.arithmeticdict import ArithmeticDict


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
    name: str (optional)
    """

    str_arrow = '->'
    latex_arrow = '\\rightarrow'
    param_char = 'k'  # convention

    def __init__(self, reac, prod, params=None, inact_reac=None, name=None):
        self.reac = reac
        self.prod = prod
        self.params = params
        self.inact_reac = inact_reac
        self.name = name

    def __eq__(lhs, rhs):
        if not isinstance(lhs, Reaction) or not isinstance(rhs, Reaction):
            return NotImplemented
        for attr in ['reac', 'prod', 'params', 'inact_reac']:
            if getattr(lhs, attr) != getattr(rhs, attr):
                return False
        return True

    def net_stoich(self, substances):
        return tuple(self.prod.get(k, 0) - self.reac.get(k, 0) - (
            0 if self.inact_reac is None else self.inact_reac.get(k, 0)
        ) for k in substances)

    def _xsolid_stoich(self, substances, xor):
        return tuple((
            0 if xor ^ k.solid else
            self.prod.get(k, 0) - self.reac.get(k, 0) - (
                0 if self.inact_reac is None else self.inact_reac.get(k, 0)
            )) for k in substances)

    def solid_stoich(self, substances):
        """ Only stoichiometry of solids """
        net = self._xsolid_stoich(substances, True)
        found1 = -1
        for idx in range(len(net)):
            if net[idx] != 0:
                if found1 == -1:
                    found1 = idx
                else:
                    raise NotImplementedError("Only one solid assumed.")
        return net, net[idx], idx

    def non_solid_stoich(self, substances):
        """ Only stoichiometry of non-solids """
        return self._xsolid_stoich(substances, False)

    def has_solids(self):
        for subst in chain(self.reac.keys(), self.prod.keys(),
                           (self.inact_reac or {}).keys()):
            if subst.solid:
                return True
        return False

    def __str__(self):
        try:
            s = ' %s=%.2g' % (self.param_char, self.params)
        except:
            s = ''
        return self._get_str('name', 'str_arrow') + s

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


def equilibrium_quotient(concs, stoich):
    if not hasattr(concs, 'ndim') or concs.ndim == 1:
        tot = 1
    else:
        tot = np.ones(concs.shape[0])
        concs = concs.T

    for nr, conc in zip(stoich, concs):
        tot *= conc**nr
    return tot


class Equilibrium(Reaction):
    """
    Represents equilibrium reaction

    See :py:class:`Reaction` for parameters
    """

    str_arrow = '<->'
    latex_arrow = '\\rightleftharpoons'
    param_char = 'k'  # convention

    def __init__(self, reac, prod, params, *args):
        if not all(arg is None for arg in args):
            raise NotImplementedError("Inactive reac/prod not implemented")
        return super(Equilibrium, self).__init__(reac, prod, params)

    def K(self, state=None):
        """ Return equilibrium quotient (possibly state dependent) """
        if callable(self.params):
            if state is None:
                raise ValueError("No state provided")
            return self.params(state)
        else:
            if state is not None:
                raise ValueError("state provided but params not callable")
            return self.params

    def Q(self, substances, concs):
        stoich = self.non_solid_stoich(substances)
        return equilibrium_quotient(concs, stoich)

    def solid_factor(self, substances, sc_concs):
        factor = 1
        for r, n in self.reac.items():
            if r.solid:
                factor *= sc_concs[substances.index(r)]**-n
        for p, n in self.prod.items():
            if p.solid:
                factor *= sc_concs[substances.index(p)]**n
        return factor

    def dimensionality(self):
        result = 0
        for r, n in self.reac.items():
            if r.solid:
                continue
            result -= n
        for p, n in self.prod.items():
            if p.solid:
                continue
            result += n
        return result

    def __rmul__(lhs, rhs):  # This works on both Py2 and Py3
        if not isinstance(rhs, int) or not isinstance(lhs, Equilibrium):
            return NotImplemented
        return Equilibrium(dict(rhs*ArithmeticDict(int, lhs.reac)),
                           dict(rhs*ArithmeticDict(int, lhs.prod)),
                           lhs.params**rhs)

    def __add__(lhs, rhs):
        keys = set()
        for key in chain(lhs.reac.keys(), lhs.prod.keys(),
                         rhs.reac.keys(), rhs.prod.keys()):
            keys.add(key)
        reac, prod = {}, {}
        for key in keys:
            n = (lhs.prod.get(key, 0) - lhs.reac.get(key, 0) +
                 rhs.prod.get(key, 0) - rhs.reac.get(key, 0))
            if n < 0:
                reac[key] = -n
            elif n > 0:
                prod[key] = n
            else:
                pass  # n == 0
        return Equilibrium(reac, prod, lhs.params * rhs.params)

    def __sub__(lhs, rhs):
        return lhs + -1*rhs


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

    def net_stoichs(self):
        return np.array([(eq.net_stoich(self.substances)) for
                         idx, eq in enumerate(self.rxns)], dtype=np.int)

    def stoichs(self, non_precip_rids=()):
        return np.array([(
            -np.array(eq.solid_stoich(self.substances)[0]) if idx
            in non_precip_rids else
            eq.non_solid_stoich(self.substances)
        ) for idx, eq in enumerate(self.rxns)], dtype=np.int)
