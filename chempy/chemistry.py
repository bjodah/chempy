# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from functools import reduce
from itertools import chain
from operator import itemgetter, mul
from collections import OrderedDict, defaultdict
import sys

from .util.arithmeticdict import ArithmeticDict
from .util.parsing import (
    formula_to_composition, mass_from_composition, to_reaction,
    formula_to_latex, formula_to_unicode, formula_to_html
)
from .units import to_unitless
from ._util import intdiv
from .util.pyutil import deprecated


class Substance(object):
    """ Class representing a chemical substance

    Parameters
    ----------
    name: str
    charge: int (optional, default: None)
        will be stored in composition[0], prefer composition when possible
    latex_name: str
    unicode_name: str
    html_name: str
    composition: dict or None (default)
        dict (int -> number) e.g. {atomic number: count}, zero has special
        meaning (net charge)
    other_properties: dict
        free form dictionary. Could be simple such as ``{'mp': 0, 'bp': 100}``
        or considerably more involved, e.g.: ``{'diffusion_coefficient': {\
 'water': lambda T: 2.1*m**2/s/K*(T - 273.15*K)}}``

    Attributes
    ----------
    mass
        maps to other_properties, and when unavailable looks for formula.mass
    attrs
        a tuple of attribute names for serialization

    Examples
    --------
    >>> ammonium = Substance('NH4+', 1, 'NH_4^+', composition={7: 1, 1: 4},
    ...     other_properties={'mass': 18.0385, 'pKa': 9.24})
    >>> ammonium.name
    'NH4+'
    >>> ammonium.composition  # note that charge was inserted as composition[0]
    {0: 1, 1: 4, 7: 1}
    >>> ammonium.other_properties['mass']
    18.0385
    >>> ammonium.other_properties['pKa']
    9.24
    >>> ammonium.mass  # mass is a special case (also attribute)
    18.0385
    >>> ammonium.pKa
    Traceback (most recent call last):
        ...
    AttributeError: 'Substance' object has no attribute 'pKa'
    >>> nh4p = Substance.from_formula('NH4+')  # simpler
    >>> nh4p.composition == {7: 1, 1: 4, 0: 1}
    True
    >>> nh4p.latex_name
    'NH_{4}^{+}'

    """

    attrs = (
        'name', 'latex_name', 'unicode_name', 'html_name',
        'composition', 'other_properties'
    )

    @property
    def charge(self):
        """ Convenience property for accessing ``composition[0]`` """
        return self.composition.get(0, 0)  # electron (net) deficiency

    @property
    def mass(self):
        """ Convenience property for accessing ``other_properties['mass']``

        when ``other_properties['mass']`` is missing the mass is calculated
        from the :attr:`composition` using
        :func:`chempy.util.parsing.mass_from_composition`.
        """
        try:
            return self.other_properties['mass']
        except KeyError:
            if self.composition is not None:
                return mass_from_composition(self.composition)

    def molar_mass(self, units):
        """ Returns the molar mass (with units) of the substance

        Examples
        --------
        >>> nh4p = Substance.from_formula('NH4+')  # simpler
        >>> from chempy.units import default_units as u
        >>> nh4p.molar_mass(u)
        array(18.0384511) * g/mol

        """
        return self.mass*units.g/units.mol

    def __init__(self, name=None, charge=None, latex_name=None, unicode_name=None,
                 html_name=None, composition=None, other_properties=None):
        self.name = name
        self.latex_name = latex_name
        self.unicode_name = unicode_name
        self.html_name = html_name
        self.composition = composition or {}

        if 0 in self.composition:
            if charge is not None:
                raise KeyError("Cannot give both charge and composition[0]")
        else:
            if charge is not None:
                self.composition[0] = charge
        self.other_properties = other_properties or {}

    @classmethod
    def from_formula(cls, formula, **kwargs):
        """ Creates a :class:`Substance` instance from its formula

        Parameters
        ----------
        formula: str
            e.g. 'Na+', 'H2O', 'Fe(CN)6-4'
        \*\*kwargs:
            keyword arguments passed on to `.Substance`

        Examples
        --------
        >>> NH3 = Substance.from_formula('NH3')
        >>> NH3.composition == {1: 3, 7: 1}
        True
        >>> '%.2f' % NH3.mass
        '17.03'
        >>> NH3.charge
        0
        >>> NH3.latex_name
        'NH_{3}'

        """

        return cls(formula, latex_name=formula_to_latex(formula),
                   unicode_name=formula_to_unicode(formula),
                   html_name=formula_to_html(formula),
                   composition=formula_to_composition(formula),
                   **kwargs)

    def __repr__(self):
        kw = ['name=' + self.name + ', ...']  # Too verbose to print all
        return "<{}({})>".format(self.__class__.__name__, ','.join(kw))

    def __str__(self):
        return str(self.name)

    def _repr_html_(self):
        return self.html_name

    @staticmethod
    def composition_keys(substance_iter):
        """ Occuring :attr:`composition` keys among a series of substances """
        keys = set()
        for s in substance_iter:
            for k in s.composition.keys():
                keys.add(k)
        return sorted(keys)


class Species(Substance):
    """ Substance belonging to a phase

    Species extends :class:`Substance` with the new attribute :attr:`phase_idx`

    Attributes
    ----------
    phase_idx: int
        Index of the phase (default is 0)
    """
    def __init__(self, *args, **kwargs):
        phase_idx = kwargs.pop('phase_idx', 0)
        super(Species, self).__init__(*args, **kwargs)
        self.phase_idx = phase_idx

    @property
    @deprecated(last_supported_version='0.3.0', will_be_missing_in='0.5.0')
    def precipitate(self):
        """ deprecated attribute, provided for compatibility for now """
        return self.phase_idx > 0

    @classmethod
    def from_formula(cls, formula, phases=('(s)', '(l)', '(g)'),
                     default_phase_idx=0, **kwargs):
        """ Create a :class:`Species` instance from its formula

        Analogous to :meth:`Substance.from_formula` but with the addition that
        phase_idx is determined from the formula (and a mapping provided by
        ``phases``)

        Parameters
        ----------
        formula: str
            e.g. 'H2O', 'NaCl(s)', 'CO2(aq)', 'CO2(g)'
        phases: iterable of str or dict mapping str -> int
            if not in \*\*kwargs, ``phase_idx`` is determined from the suffix
            of ``formula`` where the suffixes is mapped from phases:
                if ``phases`` is a dictionary:
                    ``phase_idx = phases[suffix]``
                else:
                    ``phase_idx = phases.index(suffix) + 1``
            and if suffixes is missing in phases phase_idx is taken to be 0
        default_phase_idx: int or None (default: 0)
            If ``default_phase_idx`` is ``None``, ``ValueError`` is raised for
                unkown suffixes.
            Else ``default_phase_idx`` is used as ``phase_idx`` in those cases.
        \*\*kwargs:
            Keyword arguments passed on.

        Examples
        --------
        >>> water = Species.from_formula('H2O')
        >>> water.phase_idx
        0
        >>> NaCl = Species.from_formula('NaCl(s)')
        >>> NaCl.phase_idx
        1
        >>> Hg_l = Species.from_formula('Hg(l)')
        >>> Hg_l.phase_idx
        2
        >>> CO2g = Species.from_formula('CO2(g)')
        >>> CO2g.phase_idx
        3
        >>> CO2aq = Species.from_formula('CO2(aq)', default_phase_idx=None)
        Traceback (most recent call last):
            ...
        ValueError: Could not determine phase_idx
        >>> CO2aq = Species.from_formula('CO2(aq)')
        >>> CO2aq.phase_idx
        0
        >>> CO2aq = Species.from_formula('CO2(aq)', ['(aq)'],
        ...     default_phase_idx=None)
        >>> CO2aq.phase_idx
        1
        >>> Species.from_formula('CO2(aq)', {'(aq)': 0}, None).phase_idx
        0


        Raises
        ------
        ValueError:
            if ``default_phase_idx`` is ``None`` and no suffix found in phases

        """
        if 'phase_idx' in kwargs:
            p_i = kwargs.pop('phase_idx')
        else:
            p_i = None
            if isinstance(phases, dict):
                for k, v in phases.items():
                    if formula.endswith(k):
                        p_i = v
                        break
            else:
                for idx, phase in enumerate(phases):
                    if formula.endswith(phase):
                        p_i = idx + 1
                        break
            if p_i is None:
                if default_phase_idx is None:
                    raise ValueError("Could not determine phase_idx")
                else:
                    p_i = default_phase_idx
        return super(Species, cls).from_formula(
            formula, phase_idx=p_i, **kwargs)


@deprecated(last_supported_version='0.3.0',
            will_be_missing_in='0.5.0', use_instead=Species)
class Solute(Substance):
    """ [DEPRECATED] Use `.Species` instead

    Counter-intuitive to its name Solute has an additional
    property 'precipitate'

    """

    def __init__(self, *args, **kwargs):
        precipitate = kwargs.pop('precipitate', False)
        Substance.__init__(self, *args, **kwargs)
        self.precipitate = precipitate

    @classmethod
    def from_formula(cls, formula, **kwargs):
        if formula.endswith('(s)'):
            kwargs['precipitate'] = True
        return cls(formula, latex_name=formula_to_latex(formula),
                   unicode_name=formula_to_unicode(formula),
                   html_name=formula_to_html(formula),
                   composition=formula_to_composition(formula),
                   **kwargs)


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
    reac : dict (str -> int)
    prod : dict (str -> int)
    param : float or callable
    inact_reac : dict (optional)
    inact_prod : dict (optional)
    name : str (optional)
    k : deprecated (alias for param)
    ref : object
        reference
    other_properties : dict (optional)
    checks : iterable of str
        raises value error if any method check_%s returns False
        for all %s in checks.

    Examples
    --------
    >>> r = Reaction({'H2': 2, 'O2': 1}, {'H2O': 2})
    >>> r.keys() == {'H2', 'O2', 'H2O'}
    True
    >>> r.order()
    3
    >>> r.net_stoich(['H2', 'H2O', 'O2'])
    (-2, 2, -1)
    >>> print(r)
    2 H2 + O2 -> 2 H2O; None

    """

    str_arrow = '->'
    latex_arrow = r'\rightarrow'
    unicode_arrow = u'→'
    html_arrow = '&rarr;'
    param_char = 'k'  # convention

    def __init__(
            self, reac, prod, param=None, inact_reac=None, inact_prod=None,
            name=None, ref=None, other_properties=None,
            checks=('any_effect', 'all_positive', 'all_integral')):
        self.reac = reac
        self.prod = prod
        self.param = param
        self.inact_reac = inact_reac or {}
        self.inact_prod = inact_prod or {}
        self.name = name
        self.ref = ref
        self.other_properties = other_properties or {}

        from .kinetics.rates import RateExpr
        if isinstance(self.param, RateExpr) and self.param.rxn is None:
            self.param.rxn = self  # register instance in rate expression

        for check in checks:
            if not getattr(self, 'check_'+check)():
                raise ValueError("Check failed %s" % check)

    @classmethod
    def from_string(cls, string, substance_keys, globals_=None):
        """ Parses a string into an instance

        Parameters
        ----------
        string : str
            string representation of the reaction
        substance_keys : iterable of strings or string
        globals_ : dict (optional)
            dict for eval for (default: None -> {'chempy': chempy})

        Examples
        --------
        >>> r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", 'H2O H+ OH-')
        >>> r.reac == {'H2O': 1} and r.prod == {'H+': 1, 'OH-': 1}
        True
        >>> r2 = Reaction.from_string("2 H2O -> 2 H2 + O2", 'H2O H2 O2')
        >>> r2.reac == {'H2O': 2} and r2.prod == {'H2': 2, 'O2': 1}
        True

        Notes
        -----
        :func:`chempy.util.parsing.to_reaction` is used which in turn calls
        eval which is a severe security concern for untrusted input.

        """
        if isinstance(substance_keys, str):
            if ' ' in substance_keys:
                substance_keys = substance_keys.split()
        return to_reaction(string, substance_keys, cls.str_arrow, cls, globals_)

    def check_any_effect(self):
        """ Checks if the reaction has any effect """
        if not any(self.net_stoich(self.keys())):
            return False
        return True

    def check_all_positive(self):
        """ Checks if all stoichiometric coefficients are positive """
        for cont in (self.reac, self.prod, self.inact_reac, self.inact_prod):
            for v in cont.values():
                if v < 0:
                    return False
            return True

    def check_all_integral(self):
        """ Checks if all stoichiometric coefficents are integers """
        for cont in (self.reac, self.prod, self.inact_reac, self.inact_prod):
            for v in cont.values():
                if v != type(v)(int(v)):
                    return False
        return True

    def __eq__(lhs, rhs):
        if not isinstance(lhs, Reaction) or not isinstance(rhs, Reaction):
            return NotImplemented
        for attr in ['reac', 'prod', 'param', 'inact_reac', 'inact_prod']:
            if getattr(lhs, attr) != getattr(rhs, attr):
                return False
        return True

    def order(self):
        """ Sum of (active) reactant stoichiometries """
        return sum(self.reac.values())

    def keys(self):
        return set(chain(self.reac.keys(), self.prod.keys(),
                         self.inact_reac.keys(), self.inact_prod.keys()))

    def net_stoich(self, substance_keys):
        """ Per substance net stoichiometry tuple (active & inactive) """
        return tuple(self.prod.get(k, 0) -
                     self.reac.get(k, 0) +
                     self.inact_prod.get(k, 0) -
                     self.inact_reac.get(k, 0) for k in substance_keys)

    def all_reac_stoich(self, substances):
        """ Per substance reactant stoichiometry tuple (active & inactive) """
        return tuple(self.reac.get(k, 0) + self.inact_reac.get(k, 0)
                     for k in substances)

    def all_prod_stoich(self, substances):
        """ Per substance product stoichiometry tuple (active & inactive) """
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

    def _get_str_parts(self, name_attr, arrow_attr, substances, _str=str):
        nullstr, space = _str(''), _str(' ')
        reac, prod, i_reac, i_prod = [[
            ((_str(v)+space) if v > 1 else nullstr) + _str(getattr(
                substances[k], name_attr, k))
            for k, v in filter(itemgetter(1), d.items())
        ] for d in (self.reac, self.prod, self.inact_reac,
                    self.inact_prod)]
        r_str = _str(" + ").join(sorted(reac))
        ir_str = (_str(' (+ ') + _str(" + ").join(sorted(i_reac)) + _str(')')
                  if len(i_reac) > 0 else nullstr)
        arrow_str = getattr(self, arrow_attr)
        p_str = _str(" + ").join(sorted(prod))
        ip_str = (_str(' (+ ') + _str(" + ").join(sorted(i_prod)) + _str(')')
                  if len(i_prod) > 0 else nullstr)
        return r_str, ir_str, arrow_str, p_str, ip_str

    def _get_str(self, *args, **kwargs):
        _str = kwargs.get('_str', str)
        return _str("{}{} {} {}{}").format(*self._get_str_parts(*args, **kwargs))

    def _str_param(self, fmt='%.3g'):
        try:
            return (fmt + ' %s') % (self.param, self.param.dimensionality)
        except AttributeError:
            try:
                return fmt % self.param
            except TypeError:
                return str(self.param)

    def __str__(self):
        s = '; ' + self._str_param()
        return self._get_str('name', 'str_arrow', {
            k: k for k in chain(self.reac.keys(), self.prod.keys(),
                                self.inact_reac.keys(),
                                self.inact_prod.keys())
        }) + s

    def latex(self, substances):
        """ Returns a LaTeX representation of the reaction

        Examples
        --------
        >>> keys = 'H2O H+ OH-'.split()
        >>> subst = {k: Substance.from_formula(k) for k in keys}
        >>> r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", subst)
        >>> r.latex(subst) == r'H_{2}O \\rightarrow H^{+} + OH^{-}'
        True

        """
        return self._get_str('latex_name', 'latex_arrow', substances)

    def unicode(self, substances):
        u""" Returns a unicode string representation of the reaction

        Examples
        --------
        >>> keys = 'H2O H+ OH-'.split()
        >>> subst = {k: Substance.from_formula(k) for k in keys}
        >>> r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", subst)
        >>> r.unicode(subst) == u'H₂O → H⁺ + OH⁻'
        True

        """
        return self._get_str('unicode_name', 'unicode_arrow', substances,
                             _str=str if sys.version_info[0] > 2 else unicode)

    def html(self, substances):
        """ Returns a HTML representation of the reaction

        Examples
        --------
        >>> keys = 'H2O H+ OH-'.split()
        >>> subst = {k: Substance.from_formula(k) for k in keys}
        >>> r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", subst)
        >>> r.html(subst)
        'H<sub>2</sub>O &rarr; H<sup>+</sup> + OH<sup>-</sup>'

        """
        return self._get_str('html_name', 'html_arrow', substances)

    def _violation(self, substances, attr):
        net = 0.0
        for substance, coeff in zip(substances.values(),
                                    self.net_stoich(substances.keys())):
            net += getattr(substance, attr) * coeff
        return net

    def mass_balance_violation(self, substances):
        """ Net amount of mass produced

        Parameters
        ----------
        substances: dict

        Returns
        -------
        float: amount of net mass produced/consumed

        """
        return self._violation(substances, 'mass')

    def charge_neutrality_violation(self, substances):
        """ Net amount of charge produced

        Parameters
        ----------
        substances: dict

        Returns
        -------
        float: amount of net charge produced/consumed

        """
        return self._violation(substances, 'charge')

    def composition_violation(self, substances, composition_keys=None):
        """ Net amount of constituent produced

        If composition keys correspond to conserved entities e.g. atoms
        in chemical reactions, this function should return a list of zeros.
        """
        if composition_keys is None:
            composition_keys = Substance.composition_keys(substances.values())
        net = [0]*len(composition_keys)
        for substance, coeff in zip(substances.values(),
                                    self.net_stoich(substances.keys())):
            for idx, key in enumerate(composition_keys):
                net[idx] += substance.composition.get(key, 0) * coeff
        return net


def equilibrium_quotient(concs, stoich):
    """ Calculates the equilibrium quotient of an equilbrium

    Parameters
    ----------
    concs: array_like
        per substance concentration
    stoich: iterable of integers
        per substance stoichiometric coefficient

    Examples
    --------
    >>> '%.12g' % equilibrium_quotient([1.0, 1e-7, 1e-7], [-1, 1, 1])
    '1e-14'

    """
    import numpy as np

    if not hasattr(concs, 'ndim') or concs.ndim == 1:
        tot = 1
    else:
        tot = np.ones(concs.shape[0])
        concs = concs.T

    for nr, conc in zip(stoich, concs):
        tot *= conc**nr
    return tot


class Equilibrium(Reaction):
    """ Represents an equilibrium reaction

    See :class:`Reaction` for parameters

    """
    str_arrow = '='
    latex_arrow = r'\rightleftharpoons'
    unicode_arrow = '⇌'
    html_arrow = '&harr;'
    param_char = 'K'  # convention

    def as_reactions(self, state=None, kf=None, kb=None, units=None):
        """ Creates a forward and backward :class:`Reaction` pair """
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
        """ Calculates the equilibrium qoutient """
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

    def __rmul__(self, other):  # This works on both Py2 and Py3
        try:
            other_is_int = other.is_integer
        except AttributeError:
            other_is_int = isinstance(other, int)
        if not other_is_int or not isinstance(self, Equilibrium):
            return NotImplemented
        param = None if self.param is None else self.param**other
        if other < 0:
            other *= -1
            flip = True
        else:
            flip = False
        reac = dict(other*ArithmeticDict(int, self.reac))
        prod = dict(other*ArithmeticDict(int, self.prod))
        inact_reac = dict(other*ArithmeticDict(int, self.inact_reac))
        inact_prod = dict(other*ArithmeticDict(int, self.inact_prod))
        if flip:
            reac, prod = prod, reac
            inact_reac, inact_prod = inact_prod, inact_reac
        return Equilibrium(reac, prod, param,
                           inact_reac=inact_reac, inact_prod=inact_prod)

    def __neg__(self):
        return -1*self

    def __mul__(self, other):
        return other*self

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
        if (lhs.param, rhs.param) == (None, None):
            param = None
        else:
            param = lhs.param * rhs.param
        return Equilibrium(reac, prod, param)

    def __sub__(lhs, rhs):
        return lhs + -1*rhs

    @staticmethod
    def eliminate(rxns, wrt):
        """ Linear combination coefficients for elimination of a substance

        Parameters
        ----------
        rxns : iterable of Equilibrium instances
        wrt : str (substance key)

        Examples
        --------
        >>> e1 = Equilibrium({'Cd+2': 4, 'H2O': 4}, {'Cd4(OH)4+4': 1, 'H+': 4}, 10**-32.5)
        >>> e2 = Equilibrium({'Cd(OH)2(s)': 1}, {'Cd+2': 1, 'OH-': 2}, 10**-14.4)
        >>> Equilibrium.eliminate([e1, e2], 'Cd+2')
        [1, 4]
        >>> print(1*e1 + 4*e2)
        4 Cd(OH)2(s) + 4 H2O = 4 H+ + 8 OH- + Cd4(OH)4+4; 7.94e-91

        """
        import sympy
        viol = [r.net_stoich([wrt])[0] for r in rxns]
        factors = defaultdict(int)
        for v in viol:
            for f in sympy.primefactors(v):
                factors[f] = max(factors[f], sympy.Abs(v//f))
        rcd = reduce(mul, (k**v for k, v in factors.items()))
        viol[0] *= -1
        return [rcd//v for v in viol]

    def cancel(self, rxn):
        """ multiplier of how many times rxn can be added/subtracted

        Parameters
        ----------
        rxn : Equilibrium

        Examples
        --------
        >>> e1 = Equilibrium({'Cd(OH)2(s)': 4, 'H2O': 4},
        ...                  {'Cd4(OH)4+4': 1, 'H+': 4, 'OH-': 8}, 7.94e-91)
        >>> e2 = Equilibrium({'H2O': 1}, {'H+': 1, 'OH-': 1}, 10**-14)
        >>> e1.cancel(e2)
        -4
        >>> print(e1 - 4*e2)
        4 Cd(OH)2(s) = 4 OH- + Cd4(OH)4+4; 7.94e-35

        """
        keys = rxn.keys()
        s1 = self.net_stoich(keys)
        s2 = rxn.net_stoich(keys)
        candidate = float('inf')
        for v1, v2 in zip(s1, s2):
            r = intdiv(-v1, v2)
            candidate = min(candidate, r, key=abs)
        return candidate


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
    check_balance: bool (default: None)
        if None => True if all substances has composition attribute.
    substance_factory: callback
        e.g. :meth:`Substance.from_formula`

    Attributes
    ----------
    rxns: list of objects
        sequence of :class:`Reaction` instances
    substances: OrderedDict or string or iterable of strings/Substance
        mapping substance name to substance index
    ns: int
        number of substances
    nr: int
        number of reactions

    Examples
    --------
    >>> from chempy import Reaction
    >>> r1 = Reaction({'R1': 1}, {'P1': 1}, 42.0)
    >>> rsys = ReactionSystem([r1], 'R1 P1')
    >>> rsys.as_per_substance_array({'R1': 2, 'P1': 3})
    array([ 2.,  3.])

    Raises
    ------
    ValueError
        When any reaction occurs more than once

    """

    def __init__(self, rxns, substances, name=None, check_balance=None,
                 substance_factory=Substance):
        self.rxns = rxns
        if substances is None:
            substances = set.union(*[set(rxn.keys()) for rxn in self.rxns])

        if isinstance(substances, OrderedDict):
            self.substances = substances
        elif isinstance(substances, str):
            if ' ' in substances:
                substances = substances.split()
            self.substances = OrderedDict([
                (s, substance_factory(s)) for s in substances])
        else:
            try:
                self.substances = OrderedDict([(s.name, s) for s in substances])
            except:
                try:
                    self.substances = OrderedDict(substances)
                except ValueError:
                    self.substances = OrderedDict((k, substance_factory(k)) for k in substances)
        self._sanity_check()
        self.name = name
        if check_balance is None:
            for subst in self.substances.values():
                if subst.composition is None:
                    check_balance = False
                    break
            else:
                check_balance = True
        if check_balance:
            self._balance_check()
        self._duplicate_check()

    def _repr_html_(self):
        def _format(r):
            return r.html(self.substances) + '; ' + r._str_param()
        return '<br>'.join(map(_format, self.rxns))

    def _duplicate_check(self):
        for i1, rxn1 in enumerate(self.rxns):
            for i2, rxn2 in enumerate(self.rxns[i1+1:]):
                if rxn1 == rxn2:
                    raise ValueError("Duplicate reactions %d & %d" % (i1, i2))

    def _balance_check(self):
        for rxn in self.rxns:
            for net in rxn.composition_violation(self.substances):
                if net != 0:
                    raise ValueError("Reaction not balanced: %s" % str(rxn))

    def _sanity_check(self):
        for rxn in self.rxns:
            for key in chain(rxn.reac, rxn.prod, rxn.inact_reac,
                             rxn.inact_prod):
                if key not in self.substances:
                    raise ValueError("Unkown substance: %s" % key)

    def substance_names(self):
        """ Returns a tuple of the substances' names """
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
        """ Returns list of per reaction ``param`` value """
        return [rxn.param for rxn in self.rxns]

    def as_per_substance_array(self, cont, dtype='float64', unit=None):
        """ Turns a dict into an ordered array """
        import numpy as np
        if unit is not None:
            cont = to_unitless(cont, unit)
        if isinstance(cont, np.ndarray):
            pass
        elif isinstance(cont, dict):
            substance_keys = self.substances.keys()
            for k in cont:
                if k not in substance_keys:
                    raise KeyError("Unkown substance key: %s" % k)
            cont = [cont[k] for k in substance_keys]

        cont = np.asarray(cont, dtype=dtype)
        if cont.shape[-1] != self.ns:
            raise ValueError("Incorrect size")
        return cont*(unit if unit is not None else 1)

    def as_per_substance_dict(self, arr):
        return dict(zip(self.substances.keys(), arr))

    def as_substance_index(self, sbstnc):
        """ Returns the index of a Substance in the system"""
        if isinstance(sbstnc, int):
            return sbstnc
        else:
            return list(self.substances.keys()).index(sbstnc)

    def _stoichs(self, attr, keys=None):
        import numpy as np
        if keys is None:
            keys = self.substances.keys()
        # dtype: see https://github.com/sympy/sympy/issues/10295
        return np.array([(getattr(eq, attr)(keys)) for eq in self.rxns],
                        dtype=object)

    def net_stoichs(self, keys=None):
        return self._stoichs('net_stoich', keys)

    def all_reac_stoichs(self, keys=None):
        return self._stoichs('all_reac_stoich', keys)

    def all_prod_stoichs(self, keys=None):
        return self._stoichs('all_prod_stoich', keys)

    def stoichs(self, non_precip_rids=()):  # TODO: rename to cond_stoichs
        """ Conditional stoichiometries depending on precipitation status """
        # dtype: see https://github.com/sympy/sympy/issues/10295
        import numpy as np
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
