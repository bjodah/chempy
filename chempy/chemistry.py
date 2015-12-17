# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np

from itertools import chain
from operator import itemgetter
from collections import defaultdict, OrderedDict
import warnings

from .util.arithmeticdict import ArithmeticDict


def elements(formula):
    """
    Returns a dict mapping {periodictable.core.Element: int}

    Parameters
    ----------
    formula: periodictable.formulas.Formula

    Returns
    -------
    defaultdict(int) with atom numbers as keys and respective number
        of occuring atoms in formula.
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
    charge: int (optional, default: None)
        will be stored in composition[0].
    latex_name: str
    formula: carries dict attribute "atoms"
    composition: dict or None (default)
        dict (int -> number) e.g. {atomic number: count}, zero has special meaning (net charge)
    other_properties: dict
        free form dictionary. Could be simple such as ``{'mp': 0, 'bp': 100}``
        or considerably more involved, e.g.: ``{'diffusion_coefficient': {
            'water': lambda T: 2.1*m**2/s/K*(T - 273.15*K)}}``
    """
    attrs = ('name', 'mass', 'latex_name', 'formula', 'composition',
             'other_properties')

    @property
    def charge(self):
        return self.composition[0]  # electron deficiency

    @property
    def mass(self):
        try:
            return self.formula.mass
        except:
            return None  # we could use atomic masses here

    def __init__(self, name=None, charge=None, latex_name=None, formula=None,
                 composition=None, other_properties=None):
        self.name = name
        self.latex_name = latex_name

        if isinstance(formula, str):
            import periodictable
            formula = periodictable.formula(formula)
        self.formula = formula

        if composition is None and formula is not None:
            composition = {
                k.number: v for k, v in elements(formula).items()}
        self.composition = composition or {}

        if 0 in self.composition:
            if charge is not None:
                raise KeyError("Cannot give both charge and composition[0]")
        else:
            if charge is None:
                try:
                    charge = self.formula.charge
                except AttributeError:
                    pass
            if charge is not None:
                self.composition[0] = charge
        self.other_properties = other_properties

    def __repr__(self):
        kw = ['name=' + self.name + ', ...']  # Too verbose to print all
        return "<{}({})>".format(self.__class__.__name__, ','.join(kw))

    def __str__(self):
        return str(self.name)

    @staticmethod
    def composition_keys(substance_iter):
        keys = set()
        for s in substance_iter:
            for k in s.composition.keys():
                keys.add(k)
        return sorted(keys)


class Solute(Substance):

    def __init__(self, *args, **kwargs):
        self.solid = kwargs.pop('solid', False)
        super(self.__class__, self).__init__(*args, **kwargs)


class Reaction(object):
    """
    Parameters
    ----------
    reac: dict (str -> int)
    prod: dict (str -> int)
    param: float or callable
    inact_reac: dict (optional)
    name: str (optional)
    k: deprecated (alias for param)
    """

    str_arrow = '->'
    latex_arrow = '\\rightarrow'
    param_char = 'k'  # convention

    def __init__(self, reac, prod, param=None, inact_reac=None, name=None,
                 k=None):
        self.reac = reac
        self.prod = prod
        if k is not None:
            if param is not None:
                raise ValueError("Got both param and k")
            param = k
            warnings.warn("Use param instead", DeprecationWarning)
        self.param = param
        self.inact_reac = inact_reac
        self.name = name

    def __eq__(lhs, rhs):
        if not isinstance(lhs, Reaction) or not isinstance(rhs, Reaction):
            return NotImplemented
        for attr in ['reac', 'prod', 'param', 'inact_reac']:
            if getattr(lhs, attr) != getattr(rhs, attr):
                return False
        return True

    def net_stoich(self, substances):
        return tuple(self.prod.get(k, 0) - self.reac.get(k, 0) - (
            0 if self.inact_reac is None else self.inact_reac.get(k, 0)
        ) for k in substances)

    def _xsolid_stoich(self, substances, xor):
        return tuple((
            0 if xor ^ v.solid else
            self.prod.get(k, 0) - self.reac.get(k, 0) - (
                0 if self.inact_reac is None else self.inact_reac.get(k, 0)
            )) for k, v in substances.items())

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
            s = ' %s=%.2g' % (self.param_char, self.param)
        except:
            s = ''
        return self._get_str('name', 'str_arrow', {
            k: k for k in self.reac.keys() + self.prod.keys()}) + s

    def _get_str(self, name_attr, arrow_attr, substances):
        reac, prod = [[
            ((str(v)+' ') if v > 1 else '') + getattr(substances[k], name_attr)
            for k, v in filter(itemgetter(1), d.items())
        ] for d in (self.reac, self.prod)]
        fmtstr = "{} " + getattr(self, arrow_attr) + " {}"
        return fmtstr.format(" + ".join(reac),
                             " + ".join(prod))

    def latex(self, substances):
        return self._get_str('latex_name', 'latex_arrow', substances)

    def _violation(self, substances, attr):
        net = 0.0
        for substance, coeff in zip(substances.values(), self.net_stoich(substances.keys())):
            net += getattr(substance, attr) * coeff
        return net

    def mass_balance_violation(self, substances):
        return self._violation(substances, 'mass')

    def charge_neutrality_violation(self, substances):
        return self._violation(substances, 'charge')

    def composition_violation(self, substances, composition_keys=None):
        if composition_keys is None:
            composition_keys = Substance.composition_keys(substances.values())
        net = [0]*len(composition_keys)
        for substance, coeff in zip(substances.values(), self.net_stoich(substances.keys())):
            for idx, key in enumerate(composition_keys):
                net[idx] += substance.composition.get(key, 0) * coeff
        return net



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

    def __init__(self, reac, prod, param, *args):
        if not all(arg is None for arg in args):
            raise NotImplementedError("Inactive reac/prod not implemented")
        return super(Equilibrium, self).__init__(reac, prod, param)

    def K(self, state=None):
        """ Return equilibrium quotient (possibly state dependent) """
        if callable(self.param):
            if state is None:
                raise ValueError("No state provided")
            return self.param(state)
        else:
            if state is not None:
                raise ValueError("state provided but param not callable")
            return self.param

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
                           lhs.param**rhs)

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
        return Equilibrium(reac, prod, lhs.param * rhs.param)

    def __sub__(lhs, rhs):
        return lhs + -1*rhs


class ReactionSystem(object):
    """
    Collection of reactions forming a system (model).

    Parameters
    ----------
    rxns: sequence
         sequence of :py:class:`Reaction` instances
    substances: OrderedDict

    substance_names: iterable of strings (optional)
        if not provided substacens.keys() will be used.
    name: string (optional)
         Name of ReactionSystem (e.g. model name / citation key)

    Attributes
    ----------
    rxns: list of objects
        sequence of :class:`Reaction` instances
    species_names: set of strings
        names of occurring species
    ns: int
        number of species
    nr: int
        number of reactions

    """

    def __init__(self, rxns, substances, name=None):
        self.rxns = rxns
        if isinstance(substances, OrderedDict):
            self.substances = substances
        else:
            try:
                self.substances = OrderedDict(substances)
            except:
                self.substances = OrderedDict([(s.name, s) for s in substances])
        self.name = name

    @property
    def nr(self):
        """ Number of reactions """
        return len(self.rxns)

    @property
    def ns(self):
        """ Number of substances """
        return len(self.substances)

    def params(self):
        return [rxn.param for rxn in self.rxns]

    def as_per_substance_array(self, cont):
        """ Turns e.g. a dict into an ordered array """
        if isinstance(cont, np.ndarray):
            pass
        elif isinstance(cont, dict):
            cont = [cont[k] for k in self.substances.keys()]
        cont = np.asarray(cont)
        if cont.shape[-1] != self.ns:
            raise ValueError("Incorrect size")
        return cont

    def as_substance_index(self, sbstnc):
        """ Returns the index of a Substance in the system"""
        if isinstance(sbstnc, int):
            return sbstnc
        else:
            return self.substances.keys().index(sbstnc)

    def net_stoichs(self):
        return np.array([(eq.net_stoich(self.substances.keys())) for
                         idx, eq in enumerate(self.rxns)], dtype=np.int)

    def stoichs(self, non_precip_rids=()):
        return np.array([(
            -np.array(eq.solid_stoich(self.substances)[0]) if idx
            in non_precip_rids else
            eq.non_solid_stoich(self.substances)
        ) for idx, eq in enumerate(self.rxns)], dtype=np.int)

    def obeys_mass_balance(self):
        """ Returns True if all reactions obeys mass balance, else False. """
        for rxn in self.rxns:
            if rxn.mass_balance_violation(self.substances) != 0:
                return False
        return True

    def obeys_charge_neutrality(self):
        """ Returns False if any reaction violate charge neutrality. """
        for rxn in self.rxns:
            if rxn.charge_neutrality_violation(self.substances) != 0:
                return False
        return True
