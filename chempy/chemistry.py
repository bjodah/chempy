# -*- coding: utf-8 -*-

from __future__ import division

from itertools import chain
from operator import itemgetter
from collections import defaultdict, OrderedDict
import warnings

import numpy as np

from .arrhenius import arrhenius_equation
from .util.arithmeticdict import ArithmeticDict
from .util.pyutil import defaultnamedtuple
from .units import to_unitless, default_constants, default_units


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
        dict (int -> number) e.g. {atomic number: count}, zero has special
        meaning (net charge)
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
        self.other_properties = other_properties or {}

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
        self.precipitate = kwargs.pop('precipitate', False)
        super(self.__class__, self).__init__(*args, **kwargs)


class Reaction(object):
    """ Class representing a chemical reaction

    Consider for example:

        A + R --> A + P; r = k*A*R

    this would be represented as ``Reaction({'A': 1, 'R': 1},
    {'A': 1, 'P': 1}, param=k)``. Some reactions have a larger
    stoichiometric coefficient than what appears in the rate
    expression, e.g.:

        5*C1 + C2 --> B; r = k*C1*C2

    this can be represented as ``Reaction({'C1': 1, 'C2': 1},
    {'B': 1}, inact_reac={'C1': 4}, param=k)``.

    The rate constant information in ``param`` may be callable (with a single
    argument representing the state, e.g. temperature)

    Additional data may be stored in the ``other_properties`` dict.


    Parameters
    ----------
    reac: dict (str -> int)
    prod: dict (str -> int)
    param: float or callable
    inact_reac: dict (optional)
    inact_prod: dict (optional)
    name: str (optional)
    k: deprecated (alias for param)
    ref: object
        reference
    other_properties: dict (optional)
    """

    str_arrow = '->'
    latex_arrow = r'\rightarrow'
    param_char = 'k'  # convention

    def __init__(self, reac, prod, param=None, inact_reac=None,
                 inact_prod=None, name=None, k=None, ref=None,
                 other_properties=None):
        self.reac = reac
        self.prod = prod
        if k is not None:
            if param is not None:
                raise ValueError("Got both param and k")
            param = k
            warnings.warn("Use param instead", DeprecationWarning)
        self.param = param
        self.inact_reac = inact_reac or {}
        self.inact_prod = inact_prod or {}
        self.name = name
        self.ref = ref
        self.other_properties = other_properties or {}

    def __eq__(lhs, rhs):
        if not isinstance(lhs, Reaction) or not isinstance(rhs, Reaction):
            return NotImplemented
        for attr in ['reac', 'prod', 'param', 'inact_reac', 'inact_prod']:
            if getattr(lhs, attr) != getattr(rhs, attr):
                return False
        return True

    def order(self):
        return sum(self.reac.values())

    def net_stoich(self, substance_keys):
        return tuple(self.prod.get(k, 0) -
                     self.reac.get(k, 0) +
                     self.inact_prod.get(k, 0) -
                     self.inact_reac.get(k, 0) for k in substance_keys)

    def all_reac_stoich(self, substances):
        return tuple(self.reac.get(k, 0) + self.inact_reac.get(k, 0)
                     for k in substances)

    def all_prod_stoich(self, substances):
        return tuple(self.prod.get(k, 0) + self.inact_prod.get(k, 0)
                     for k in substances)

    def _xprecipitate_stoich(self, substances, xor):
        return tuple((
            0 if xor ^ v.precipitate else
            self.prod.get(k, 0) + self.inact_prod.get(k, 0) -
            self.reac.get(k, 0) - self.inact_reac.get(k, 0)
        ) for k, v in substances.items())

    def precipitate_stoich(self, substances):
        """ Only stoichiometry of precipitates """
        net = self._xprecipitate_stoich(substances, True)
        found1 = -1
        for idx in range(len(net)):
            if net[idx] != 0:
                if found1 == -1:
                    found1 = idx
                else:
                    raise NotImplementedError("Only one precipitate assumed.")
        return net, net[idx], idx

    def non_precipitate_stoich(self, substances):
        """ Only stoichiometry of non-precipitates """
        return self._xprecipitate_stoich(substances, False)

    def has_precipitates(self, substances):
        for s_name in chain(self.reac.keys(), self.prod.keys(),
                            self.inact_reac.keys(),
                            self.inact_prod.keys()):
            if substances[s_name].precipitate:
                return True
        return False

    def _get_str_parts(self, name_attr, arrow_attr, substances):
        reac, prod, i_reac, i_prod = [[
            ((str(v)+' ') if v > 1 else '') + str(getattr(
                substances[k], name_attr, k))
            for k, v in filter(itemgetter(1), d.items())
        ] for d in (self.reac, self.prod, self.inact_reac,
                    self.inact_prod)]
        r_str = " + ".join(reac)
        ir_str = ' (+ ' + " + ".join(i_reac) + ')' if len(i_reac) > 0 else ''
        arrow_str = getattr(self, arrow_attr)
        p_str = " + ".join(prod)
        ip_str = ' (+ ' + " + ".join(i_prod) + ')' if len(i_prod) > 0 else ''
        return r_str, ir_str, arrow_str, p_str, ip_str

    def _get_str(self, *args):
        return "{}{} {} {}{}".format(*self._get_str_parts(*args))

    def __str__(self):
        try:
            str_param = '%.3g %s' % (self.param, self.param.dimensionality)
        except AttributeError:
            try:
                str_param = '%.3g' % self.param
            except TypeError:
                str_param = str(self.param)
        s = '; ' + self.param_char + '=' + str_param
        return self._get_str('name', 'str_arrow', {
            k: k for k in chain(self.reac.keys(), self.prod.keys(),
                                self.inact_reac.keys(),
                                self.inact_prod.keys())
        }) + s

    def latex(self, substances):
        return self._get_str('latex_name', 'latex_arrow', substances)

    def _violation(self, substances, attr):
        net = 0.0
        for substance, coeff in zip(substances.values(),
                                    self.net_stoich(substances.keys())):
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
        for substance, coeff in zip(substances.values(),
                                    self.net_stoich(substances.keys())):
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


class ArrheniusRate(defaultnamedtuple('ArrheniusRate', 'A Ea ref', [None])):
    """
    Ea: float
        activation energy
    A: float
        preexponential prefactor (Arrhenius type eq.)
    ref: object (default: None)
        arbitrary reference (e.g. string representing citation key)
    """
    def __call__(self, T, constants=None, units=None, exp=None):
        """ See :py:func`chempy.arrhenius.arrhenius_equation`. """
        return arrhenius_equation(self.A, self.Ea, T, constants=constants,
                                  units=units, exp=exp)

    @staticmethod
    def _fmt(arg, precision, tex):
        if tex:
            unit_str = arg.dimensionality.latex.strip('$')
        else:
            from quantities.markup import config
            attr = 'unicode' if config.use_unicode else 'string'
            unit_str = getattr(arg.dimensionality, attr)
        return precision.format(float(arg.magnitude)) + " " + unit_str

    def format(self, precision, tex=False):
        try:
            str_A = self._fmt(self.A, precision, tex)
            str_Ea = self._fmt(self.Ea, precision, tex)
        except:
            str_A = precision.format(self.A)
            str_Ea = precision.format(self.Ea)
        return str_A, str_Ea

    def equation_as_string(self, precision, tex=False):
        if tex:
            return r"{}\exp \left(\frac{{{}}}{{RT}} \right)".format(*self.format(precision, tex))
        else:
            return "{}*exp({}/(R*T))".format(*self.format(precision, tex))

    def __str__(self):
        return self.equation_as_string('{0:.5g}')


class ArrheniusRateWithUnits(ArrheniusRate):
    def __call__(self, T, constants=default_constants, units=default_units,
                 exp=None):
        return super(ArrheniusRateWithUnits, self).__call__(
            T, constants, units, exp)


class Equilibrium(Reaction):
    """
    Represents equilibrium reaction

    See :py:class:`Reaction` for parameters
    """

    str_arrow = '<->'
    latex_arrow = r'\rightleftharpoons'
    param_char = 'K'  # convention

    def __init__(self, reac, prod, param, *args, **kwargs):
        ref = kwargs.pop('ref', None)
        if not all(arg is None for arg in args) or len(kwargs) > 0:
            raise NotImplementedError("Inactive reac/prod not implemented")
        return super(Equilibrium, self).__init__(reac, prod, param, ref=ref)

    def as_reactions(self, state=None, kf=None, kb=None, units=None):
        """ Creates a pair of Reaction instances """
        nb = sum(self.prod.values())
        nf = sum(self.reac.values())
        if units is None:
            if hasattr(kf, 'units') or hasattr(kb, 'units'):
                raise ValueError("units missing")
            c0 = 1
        else:
            c0 = 1*units.molar  # standard concentration IUPAC

        if kf is None:
            if kb is None:
                raise ValueError("Exactly one rate needs to be provided")
            kf = kb * self.K(state) * c0**(nb - nf)
        elif kb is None:
            kb = kf / (self.K(state) * c0**(nb - nf))
        else:
            raise ValueError("Exactly one rate needs to be provided")
        return (
            Reaction(self.reac, self.prod, kf, self.inact_reac,
                     self.inact_prod, ref=self.ref),
            Reaction(self.prod, self.reac, kb, self.inact_prod,
                     self.inact_reac, ref=self.ref)
        )

    def K(self, state=None):
        """ Return equilibrium quotient (possibly state dependent) """
        if callable(self.param):
            # For some cases a default might be obvious:
            # if state is None:
            #     raise ValueError("No state provided")
            return self.param(state)
        else:
            if state is not None:
                raise ValueError("state provided but param not callable")
            return self.param

    def Q(self, substances, concs):
        stoich = self.non_precipitate_stoich(substances)
        return equilibrium_quotient(concs, stoich)

    def precipitate_factor(self, substances, sc_concs):
        factor = 1
        for r, n in self.reac.items():
            if r.precipitate:
                factor *= sc_concs[substances.index(r)]**-n
        for p, n in self.prod.items():
            if p.precipitate:
                factor *= sc_concs[substances.index(p)]**n
        return factor

    def dimensionality(self):
        result = 0
        for r, n in self.reac.items():
            if r.precipitate:
                continue
            result -= n
        for p, n in self.prod.items():
            if p.precipitate:
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
         mapping str -> Substance instances
    name: string (optional)
         Name of ReactionSystem (e.g. model name / citation key)

    Attributes
    ----------
    rxns: list of objects
        sequence of :class:`Reaction` instances
    substances: OrderedDict
        mapping substance name to substance index
    ns: int
        number of substances
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
                self.substances = OrderedDict([(s.name, s) for
                                               s in substances])
        self._sanity_check()
        self.name = name

    def _sanity_check(self):
        for rxn in self.rxns:
            for key in chain(rxn.reac, rxn.prod, rxn.inact_reac,
                             rxn.inact_prod):
                if key not in self.substances:
                    raise ValueError("Unkown substance: %s" % key)

    def substance_names(self):
        return tuple(substance.name for substance in self.substances.values())

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

    def as_per_substance_array(self, cont, dtype=np.float64, unit=None):
        """ Turns e.g. a dict into an ordered array """
        if unit is not None:
            cont = to_unitless(cont, unit)
        if isinstance(cont, np.ndarray):
            pass
        elif isinstance(cont, dict):
            cont = [cont[k] for k in self.substances.keys()]
        cont = np.asarray(cont, dtype=dtype)
        if cont.shape[-1] != self.ns:
            raise ValueError("Incorrect size")
        return cont*(unit if unit is not None else 1)

    def as_substance_index(self, sbstnc):
        """ Returns the index of a Substance in the system"""
        if isinstance(sbstnc, int):
            return sbstnc
        else:
            return list(self.substances.keys()).index(sbstnc)

    def _stoichs(self, attr):
        # dtype: see https://github.com/sympy/sympy/issues/10295
        return np.array([(getattr(eq, attr)(self.substances.keys())) for
                         eq in self.rxns],
                        dtype=object)

    def net_stoichs(self):
        return self._stoichs('net_stoich')

    def all_reac_stoichs(self):
        return self._stoichs('all_reac_stoich')

    def all_prod_stoichs(self):
        return self._stoichs('all_prod_stoich')

    def stoichs(self, non_precip_rids=()):
        # dtype: see https://github.com/sympy/sympy/issues/10295
        return np.array([(
            -np.array(eq.precipitate_stoich(self.substances)[0]) if idx
            in non_precip_rids else
            eq.non_precipitate_stoich(self.substances)
        ) for idx, eq in enumerate(self.rxns)],
                        dtype=object)

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
