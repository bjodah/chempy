# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict, defaultdict
from functools import reduce
from itertools import chain
from operator import mul, add
import math
import warnings

from .util.arithmeticdict import ArithmeticDict
from .util._expr import Expr
from .util.periodic import mass_from_composition
from .util.parsing import (
    formula_to_composition, to_reaction,
    formula_to_latex, formula_to_unicode, formula_to_html
)

from .units import default_units, is_quantity, unit_of, to_unitless
from ._util import intdiv
from .util.pyutil import deprecated, DeferredImport, ChemPyDeprecationWarning


ReactionSystem = DeferredImport('chempy.reactionsystem', 'ReactionSystem',
                                [deprecated(use_instead='chempy.ReactionSystem')])


class Substance(object):
    """ Class representing a chemical substance

    Parameters
    ----------
    name : str
    charge : int (optional, default: None)
        Will be stored in composition[0], prefer composition when possible.
    latex_name : str
    unicode_name : str
    html_name : str
    composition : dict or None (default)
        Dictionary (int -> number) e.g. {atomic number: count}, zero has special
        meaning (net charge). Avoid using the key 0 unless you specifically mean
        net charge. The motivation behind this is that it is easier to track a
        net-charge of e.g. 6 for U(VI) than it is to remember that uranium has 92
        electrons and use 86 as the value).
    data : dict
        Free form dictionary. Could be simple such as ``{'mp': 0, 'bp': 100}``
        or considerably more involved, e.g.: ``{'diffusion_coefficient': {\
 'water': lambda T: 2.1*m**2/s/K*(T - 273.15*K)}}``.

    Attributes
    ----------
    mass
        Maps to data['mass'], and when unavailable looks for ``formula.mass``.
    attrs
        A tuple of attribute names for serialization.
    composition : dict or None
        Dictionary mapping fragment key (str) to amount (int).
    data
        Free form dictionary.

    Examples
    --------
    >>> ammonium = Substance('NH4+', 1, 'NH_4^+', composition={7: 1, 1: 4},
    ...     data={'mass': 18.0385, 'pKa': 9.24})
    >>> ammonium.name
    'NH4+'
    >>> ammonium.composition == {0: 1, 1: 4, 7: 1}  # charge represented by key '0'
    True
    >>> ammonium.data['mass']
    18.0385
    >>> ammonium.data['pKa']
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
        'composition', 'data'
    )

    def __eq__(self, other):
        for attr in self.attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False
        return True

    @property
    def charge(self):
        """ Convenience property for accessing ``composition[0]`` """
        return self.composition.get(0, 0)  # electron (net) deficiency

    @property
    def mass(self):
        """ Convenience property for accessing ``data['mass']``

        when ``data['mass']`` is missing the mass is calculated
        from the :attr:`composition` using
        :func:`chempy.util.parsing.mass_from_composition`.
        """
        try:
            return self.data['mass']
        except KeyError:
            if self.composition is not None:
                return mass_from_composition(self.composition)

    @mass.setter
    def mass(self, value):
        self.data['mass'] = value

    def molar_mass(self, units=None):
        """ Returns the molar mass (with units) of the substance

        Examples
        --------
        >>> nh4p = Substance.from_formula('NH4+')  # simpler
        >>> from chempy.units import default_units as u
        >>> nh4p.molar_mass(u)
        array(18.0384511...) * g/mol

        """
        if units is None:
            units = default_units
        return self.mass*units.g/units.mol

    def __init__(self, name=None, charge=None, latex_name=None, unicode_name=None,
                 html_name=None, composition=None, data=None):
        self.name = name
        self.latex_name = latex_name
        self.unicode_name = unicode_name
        self.html_name = html_name
        self.composition = composition

        if self.composition is not None and 0 in self.composition:
            if charge is not None:
                raise KeyError("Cannot give both charge and composition[0]")
        else:
            if charge is not None and composition is not None:
                self.composition[0] = charge
        self.data = data or {}

    @classmethod
    def from_formula(cls, formula, **kwargs):
        """ Creates a :class:`Substance` instance from its formula

        Parameters
        ----------
        formula: str
            e.g. 'Na+', 'H2O', 'Fe(CN)6-4'
        \\*\\*kwargs:
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
    def composition_keys(substance_iter, skip_keys=()):
        """ Occuring :attr:`composition` keys among a series of substances """
        keys = set()
        for s in substance_iter:
            if s.composition is None:
                continue
            for k in s.composition.keys():
                if k in skip_keys:
                    continue
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
    @deprecated(last_supported_version='0.3.0', will_be_missing_in='0.8.0')
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
            if not in \\*\\*kwargs, ``phase_idx`` is determined from the suffix
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
        \\*\\*kwargs:
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
        return cls(
            formula,
            latex_name=formula_to_latex(formula, suffixes=phases),
            unicode_name=formula_to_unicode(formula, suffixes=phases),
            html_name=formula_to_html(formula, suffixes=phases),
            composition=formula_to_composition(formula, suffixes=phases),
            phase_idx=p_i, **kwargs
        )


@deprecated(last_supported_version='0.3.0',
            will_be_missing_in='0.8.0', use_instead=Species)
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

        2 R --> A + P; r = k*A*R*R

    this would be represented as ``Reaction({'A': 1, 'R': 2},
    {'A': 2, 'P': 1}, param=k)``. Some reactions have a larger
    stoichiometric coefficient than what appears in the rate
    expression, e.g.:

        5 A + B --> C; r = k*A*B

    this can be represented as ``Reaction({'C1': 1, 'C2': 1},
    {'B': 1}, inact_reac={'C1': 4}, param=k)``.

    The rate constant information in ``param`` may be a subclass of
    :class:`chempy.kinetics.rates.RateExpr` or carry a :meth:`as_RateExpr`,
    if neither: `param` will be assumed to be a rate constant for a mass-action
    type of kinetic expression.

    Additional data may be stored in the ``data`` dict.


    Parameters
    ----------
    reac : dict (str -> int)
        If ``reac`` is a ``set``, then multiplicities are assumed to be 1.
    prod : dict (str -> int)
        If ``prod`` is a ``set``, then multiplicities are assumed to be 1.
    param : float or callable
        Special case (side-effect): if param is a subclass of
        :class:`.kinetics.rates.RateExpr` and its :attr:`rxn`
        is `None` it will be set to `self`.
    inact_reac : dict (optional)
    inact_prod : dict (optional)
    name : str (optional)
    k : deprecated (alias for param)
    ref : object
        Reference (e.g. a string containing doi number).
    data : dict (optional)
    checks : iterable of str
        Raises ``ValueError`` if any method ``check_%s`` returns False
        for all ``%s`` in ``checks``. Default: ``Reaction.default_checks``.

    Attributes
    ----------
    reac : OrderedDict
    prod : OrderedDict
    param : object
    inact_reac : OrderedDict
    inact_prod : OrderedDict
    name : str
    ref : str
    data : dict

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
    2 H2 + O2 -> 2 H2O

    """

    _cmp_attr = ('reac', 'prod', 'param', 'inact_reac', 'inact_prod')
    _all_attr = _cmp_attr + ('name', 'ref', 'data')
    _str_arrow = '->'

    param_char = 'k'  # convention
    default_checks = ('any_effect', 'all_positive', 'all_integral', 'consistent_units')

    @staticmethod
    def _init_stoich(container):
        if isinstance(container, set):
            container = {k: 1 for k in container}
        container = container or {}
        if type(container) == dict:  # we don't want isinstance here in case of OrderedDict
            container = OrderedDict(sorted(container.items(), key=lambda kv: kv[0]))
        return container

    def __init__(
            self, reac, prod, param=None, inact_reac=None, inact_prod=None,
            name=None, ref=None, data=None, checks=None):
        self.reac = self._init_stoich(reac)
        self.inact_reac = self._init_stoich(inact_reac)
        self.prod = self._init_stoich(prod)
        self.inact_prod = self._init_stoich(inact_prod)
        self.param = param
        self.name = name
        self.ref = ref
        self.data = data or {}
        for check in (self.default_checks if checks is None else checks):
            if not getattr(self, 'check_'+check)():
                raise ValueError("Check failed: '%s'" % check)

    @classmethod
    def from_string(cls, string, substance_keys=None, globals_=None, **kwargs):
        """ Parses a string into a Reaction instance

        Parameters
        ----------
        string : str
            String representation of the reaction.
        substance_keys : convertible to iterable of strings or string or None
            Used prevent e.g. misspelling.
            if str: split is invoked, if None: no checking done.
        globals_ : dict (optional)
            Dictionary for eval for (default: None -> {'chempy': chempy})
            If ``False``: no eval will be called (useful for web-apps).
        \\*\\*kwargs :
            Passed on to constructor.

        Examples
        --------
        >>> r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", 'H2O H+ OH-')
        >>> r.reac == {'H2O': 1} and r.prod == {'H+': 1, 'OH-': 1}
        True
        >>> r2 = Reaction.from_string("2 H2O -> 2 H2 + O2", 'H2O H2 O2')
        >>> r2.reac == {'H2O': 2} and r2.prod == {'H2': 2, 'O2': 1}
        True
        >>> r3 = Reaction.from_string("A -> B; 1/second", 'A B')
        >>> from chempy.units import to_unitless, default_units as u
        >>> to_unitless(r3.param, u.hour**-1)
        3600.0
        >>> r4 = Reaction.from_string("A -> B; 'k'", 'A B')
        >>> r4.param.unique_keys
        ('k',)
        >>> r5 = Reaction.from_string("A -> B; 1/molar/second", 'A B')
        Traceback (most recent call last):
            ...
        ValueError: Check failed: 'consistent_units'


        Notes
        -----
        :func:`chempy.util.parsing.to_reaction` is used which in turn calls
        :func:`eval` which is a severe security concern for untrusted input.

        """
        if isinstance(substance_keys, str):
            if ' ' in substance_keys:
                substance_keys = substance_keys.split()
        return to_reaction(string, substance_keys, cls._str_arrow, cls, globals_, **kwargs)

    def copy(self, **kwargs):
        return self.__class__(**{k: kwargs.get(k, getattr(self, k)) for k in self._all_attr})

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

    def check_consistent_units(self):
        if is_quantity(self.param):  # This will assume mass action
            try:
                to_unitless(self.param/(
                    default_units.molar**(1-self.order())/default_units.s))
            except Exception:
                return False
            else:
                return True
        else:
            return True  # the user might not be using ``chempy.units``

    def __eq__(lhs, rhs):
        if lhs is rhs:
            return True
        if not isinstance(lhs, Reaction) or not isinstance(rhs, Reaction):
            return NotImplemented
        for attr in lhs._cmp_attr:
            if getattr(lhs, attr) != getattr(rhs, attr):
                return False
        return True

    def __hash__(self):
        return sum(map(hash, (getattr(self, k) for k in ['reac', 'prod', 'param', 'inact_reac', 'inact_prod'])))

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
        return tuple(self.reac.get(k, 0) + self.inact_reac.get(k, 0) for k in substances)

    def active_reac_stoich(self, substances):
        """ Per substance reactant stoichiometry tuple (active) """
        return tuple(self.reac.get(k, 0) for k in substances)

    def all_prod_stoich(self, substances):
        """ Per substance product stoichiometry tuple (active & inactive) """
        return tuple(self.prod.get(k, 0) + self.inact_prod.get(k, 0) for k in substances)

    def active_prod_stoich(self, substances):
        """ Per substance product stoichiometry tuple (active) """
        return tuple(self.prod.get(k, 0) for k in substances)

    def _xprecipitate_stoich(self, substances, xor):
        return tuple((
            0 if xor ^ (getattr(v, 'phase_idx', 0) > 0) else
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
        return net, net[found1], found1

    def non_precipitate_stoich(self, substances):
        """ Only stoichiometry of non-precipitates """
        return self._xprecipitate_stoich(substances, False)

    def has_precipitates(self, substances):
        for s_name in chain(self.reac.keys(), self.prod.keys(), self.inact_reac.keys(), self.inact_prod.keys()):
            if getattr(substances[s_name], 'phase_idx', 0) > 0:
                return True
        return False

    def string(self, substances=None, with_param=False, **kwargs):
        """ Returns a string representation of the reaction

        Parameters
        ----------
        substances: dict
            mapping substance keys to Substance instances
        with_param: bool
            whether to print the parameter (default: False)

        Examples
        --------
        >>> r = Reaction({'H+': 1, 'Cl-': 1}, {'HCl': 1}, 1e10)
        >>> r.string(with_param=False)
        'Cl- + H+ -> HCl'

        """
        from .printing import str_
        return str_(self, substances=substances, with_param=with_param, **kwargs)

    def __str__(self):
        return self.string(with_param=True)

    def latex(self, substances, with_param=False, **kwargs):
        r""" Returns a LaTeX representation of the reaction

        Parameters
        ----------
        substances: dict
            mapping substance keys to Substance instances
        with_param: bool
            whether to print the parameter (default: False)

        Examples
        --------
        >>> keys = 'H2O H+ OH-'.split()
        >>> subst = {k: Substance.from_formula(k) for k in keys}
        >>> r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", subst)
        >>> r.latex(subst) == r'H_{2}O \rightarrow H^{+} + OH^{-}'
        True
        >>> r2 = Reaction.from_string("H+ + OH- -> H2O; 1e8/molar/second", subst)
        >>> ref = r'H^{+} + OH^{-} \rightarrow H_{2}O; 10^{8} $\mathrm{\frac{1}{(s{\cdot}M)}}$'
        >>> r2.latex(subst, with_param=True) == ref
        True

        """
        from .printing import latex
        return latex(self, substances=substances, with_param=with_param, **kwargs)

    def unicode(self, substances, with_param=False, **kwargs):
        u""" Returns a unicode string representation of the reaction

        Examples
        --------
        >>> keys = 'H2O H+ OH-'.split()
        >>> subst = {k: Substance.from_formula(k) for k in keys}
        >>> r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", subst)
        >>> r.unicode(subst) == u'H₂O → H⁺ + OH⁻'
        True
        >>> r2 = Reaction.from_string("H+ + OH- -> H2O; 1e8/molar/second", subst)
        >>> r2.unicode(subst, with_param=True) == u'H⁺ + OH⁻ → H₂O; 10⁸ 1/(s·M)'
        True

        """
        from .printing import unicode_
        return unicode_(self, substances=substances, with_param=with_param, **kwargs)

    def html(self, substances, with_param=False, **kwargs):
        """ Returns a HTML representation of the reaction

        Examples
        --------
        >>> keys = 'H2O H+ OH-'.split()
        >>> subst = {k: Substance.from_formula(k) for k in keys}
        >>> r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", subst)
        >>> r.html(subst)
        'H<sub>2</sub>O &rarr; H<sup>+</sup> + OH<sup>-</sup>'
        >>> r2 = Reaction.from_string("H+ + OH- -> H2O; 1e8/molar/second", subst)
        >>> r2.html(subst, with_param=True)
        'H<sup>+</sup> + OH<sup>-</sup> &rarr; H<sub>2</sub>O&#59; 10<sup>8</sup> 1/(s*M)'

        """
        from .printing import html
        return html(self, with_param=with_param, substances=substances, **kwargs)

    def _repr_html_(self):
        return self.html({k: k for k in self.keys()})

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

        Parameters
        ----------
        substances : dict
        composition_keys : iterable of str, ``None`` or ``True``
            When ``None`` or True: composition keys are taken from substances.
            When ``True`` the keys are also return as an extra return value

        Returns
        -------
        - If ``composition_keys == True``: a tuple: (violations, composition_keys)
        - Otherwise: violations (list of coefficients)

        """
        keys, values = zip(*substances.items())
        ret_comp_keys = composition_keys is True
        if composition_keys in (None, True):
            composition_keys = Substance.composition_keys(values)
        net = [0]*len(composition_keys)
        for substance, coeff in zip(values, self.net_stoich(keys)):
            for idx, key in enumerate(composition_keys):
                net[idx] += substance.composition.get(key, 0) * coeff
        if ret_comp_keys:
            return net, composition_keys
        else:
            return net

    def rate_expr(self):
        """ Turns self.param into a RateExpr instance (if not already)

        Default is to create a ``MassAction`` instance. The parameter will
        be used as single instance in ``unique_keys`` if it is a string,
        otherwise it will be used as ``args``.

        Examples
        --------
        >>> r = Reaction.from_string('2 A + B -> 3 C; 7')
        >>> ratex = r.rate_expr()
        >>> ratex.args[0] == 7
        True

        """
        from .util._expr import Expr
        from .kinetics import MassAction
        if isinstance(self.param, Expr):
            return self.param
        else:
            try:
                convertible = self.param.as_RateExpr
            except AttributeError:
                if isinstance(self.param, str):
                    return MassAction.fk(self.param)
                else:
                    return MassAction([self.param])
            else:
                return convertible()

    def rate(self, variables=None, backend=math, substance_keys=None, ratex=None):
        """ Evaluate the rate of a reaction

        Parameters
        ----------
        variables : dict
        backend : module, optional
        substance_keys : iterable of str, optional
        ratex : RateExpr

        Returns
        -------
        Dictionary mapping substance keys to the reactions contribution to overall rates.

        Examples
        --------
        >>> rxn1 = Reaction.from_string('2 H2 + O2 -> 2 H2O; 3')
        >>> ref1 = 3*5*5*7
        >>> rxn1.rate({'H2': 5, 'O2': 7}) == {'H2': -2*ref1, 'O2': -ref1, 'H2O': 2*ref1}
        True
        >>> from sympy import Symbol
        >>> k = Symbol('k')
        >>> rxn2 = Reaction(rxn1.reac, rxn1.prod, k)
        >>> concentrations = {key: Symbol(key) for key in set.union(set(rxn1.reac), set(rxn1.prod))}
        >>> import pprint
        >>> pprint.pprint(rxn2.rate(concentrations))
        {'H2': -2*H2**2*O2*k, 'H2O': 2*H2**2*O2*k, 'O2': -H2**2*O2*k}

        """
        if variables is None:
            variables = {}
        if substance_keys is None:
            substance_keys = self.keys()
        if ratex is None:
            ratex = self.rate_expr()

        if isinstance(ratex, Expr):
            srat = ratex(variables, backend=backend, reaction=self)
        else:
            srat = ratex
        return {k: srat*v for k, v in zip(substance_keys, self.net_stoich(substance_keys))}


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
    _str_arrow = '='
    param_char = 'K'  # convention

    def check_consistent_units(self):
        if is_quantity(self.param):  # This will assume mass action
            exponent = sum(self.prod.values()) - sum(self.reac.values())
            return unit_of(self.param, simplified=True) == unit_of(
                default_units.molar**exponent, simplified=True)
        else:
            return True  # the user might not be using ``chempy.units``

    def as_reactions(self, kf=None, kb=None, units=None, variables=None, backend=math, new_name=None):
        """ Creates a forward and backward :class:`Reaction` pair

        Parameters
        ----------
        kf : float or RateExpr
        kb : float or RateExpr
        units : module
        variables : dict, optional
        backend : module

        """
        nb = sum(self.prod.values())
        nf = sum(self.reac.values())
        if units is None:
            if hasattr(kf, 'units') or hasattr(kb, 'units'):
                raise ValueError("units missing")
            c0 = 1
        else:
            c0 = 1*units.molar  # standard concentration IUPAC

        if kf is None:
            fw_name = self.name
            bw_name = new_name
            if kb is None:
                try:
                    kf, kb = self.param
                except TypeError:
                    raise ValueError("Exactly one rate needs to be provided")
            else:
                kf = kb * self.param * c0**(nb - nf)
        elif kb is None:
            kb = kf / (self.param * c0**(nb - nf))
            fw_name = new_name
            bw_name = self.name
        else:
            raise ValueError("Exactly one rate needs to be provided")

        return (
            Reaction(self.reac, self.prod, kf, self.inact_reac,
                     self.inact_prod, ref=self.ref, name=fw_name),
            Reaction(self.prod, self.reac, kb, self.inact_prod,
                     self.inact_reac, ref=self.ref, name=bw_name)
        )

    def equilibrium_expr(self):
        """ Turns self.param into a :class:`EqExpr` instance (if not already)

        Examples
        --------
        >>> r = Equilibrium.from_string('2 A + B = 3 C; 7')
        >>> eqex = r.equilibrium_expr()
        >>> eqex.args[0] == 7
        True

        """
        from .util._expr import Expr
        from .thermodynamics import MassActionEq
        if isinstance(self.param, Expr):
            return self.param
        else:
            try:
                convertible = self.param.as_EqExpr
            except AttributeError:
                return MassActionEq([self.param])
            else:
                return convertible()

    def equilibrium_constant(self, variables=None, backend=math):
        """ Return equilibrium constant

        Parameters
        ----------
        variables : dict, optional
        backend : module, optional

        """
        return self.equilibrium_expr().eq_const(variables, backend=backend)

    def equilibrium_equation(self, variables, backend=None, **kwargs):
        return self.equilibrium_expr().equilibrium_equation(
            variables, equilibrium=self, backend=backend, **kwargs)

    @deprecated(use_instead=equilibrium_constant)
    def K(self, *args, **kwargs):
        return self.equilibrium_constant(*args, **kwargs)

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

    def dimensionality(self, substances):
        result = 0
        for r, n in self.reac.items():
            if getattr(substances[r], 'phase_idx', 0) > 0:
                continue
            result -= n
        for p, n in self.prod.items():
            if getattr(substances[p], 'phase_idx', 0) > 0:
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

    def __add__(self, other):
        keys = set()
        for key in chain(self.reac.keys(), self.prod.keys(),
                         other.reac.keys(), other.prod.keys()):
            keys.add(key)
        reac, prod = {}, {}
        for key in keys:
            n = (self.prod.get(key, 0) - self.reac.get(key, 0) +
                 other.prod.get(key, 0) - other.reac.get(key, 0))
            if n < 0:
                reac[key] = -n
            elif n > 0:
                prod[key] = n
            else:
                pass  # n == 0
        if (self.param, other.param) == (None, None):
            param = None
        else:
            param = self.param * other.param
        return Equilibrium(reac, prod, param)

    def __sub__(self, other):
        return self + -1*other

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
        4 Cd(OH)2(s) + 4 H2O = Cd4(OH)4+4 + 4 H+ + 8 OH-; 7.94e-91

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
        """ Multiplier of how many times rxn can be added/subtracted.

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
        4 Cd(OH)2(s) = Cd4(OH)4+4 + 4 OH-; 7.94e-35

        """
        keys = rxn.keys()
        s1 = self.net_stoich(keys)
        s2 = rxn.net_stoich(keys)
        candidate = float('inf')
        for v1, v2 in zip(s1, s2):
            r = intdiv(-v1, v2)
            candidate = min(candidate, r, key=abs)
        return candidate


def _solve_balancing_ilp_pulp(A):
    import pulp
    x = [pulp.LpVariable('x%d' % i, lowBound=1, cat='Integer') for i in range(A.shape[1])]
    prob = pulp.LpProblem("chempy balancing problem", pulp.LpMinimize)
    prob += reduce(add, x)
    for expr in [pulp.lpSum([x[i]*e for i, e in enumerate(row)]) for row in A.tolist()]:
        prob += expr == 0
    prob.solve()
    return [pulp.value(_) for _ in x]


def balance_stoichiometry(reactants, products, substances=None,
                          substance_factory=Substance.from_formula,
                          parametric_symbols=None, underdetermined=True):
    """ Balances stoichiometric coefficients of a reaction

    Parameters
    ----------
    reactants : iterable of reactant keys
    products : iterable of product keys
    substances : OrderedDict or string or None
        Mapping reactant/product keys to instances of :class:`Substance`.
    substance_factory : callback
    parametric_symbols : generator of symbols
        Used to generate symbols for parametric solution for
        under-determined system of equations. Default is numbered "x-symbols" starting
        from 1.
    underdetermined : bool
        Allows to find a non-unique solution (in addition to a constant factor
        across all terms). Set to ``False`` to disallow (raise ValueError) on
        e.g. "C + O2 -> CO + CO2". Set to ``None`` if you want the symbols replaced
        so that the coefficients are the smallest possible positive (non-zero) integers.

    Examples
    --------
    >>> ref = {'C2H2': 2, 'O2': 3}, {'CO': 4, 'H2O': 2}
    >>> balance_stoichiometry({'C2H2', 'O2'}, {'CO', 'H2O'}) == ref
    True
    >>> ref2 = {'H2': 1, 'O2': 1}, {'H2O2': 1}
    >>> balance_stoichiometry('H2 O2'.split(), ['H2O2'], 'H2 O2 H2O2') == ref2
    True
    >>> reac, prod = 'CuSCN KIO3 HCl'.split(), 'CuSO4 KCl HCN ICl H2O'.split()
    >>> Reaction(*balance_stoichiometry(reac, prod)).string()
    '4 CuSCN + 7 KIO3 + 14 HCl -> 4 CuSO4 + 7 KCl + 4 HCN + 7 ICl + 5 H2O'
    >>> balance_stoichiometry({'Fe', 'O2'}, {'FeO', 'Fe2O3'}, underdetermined=False)
    Traceback (most recent call last):
        ...
    ValueError: The system was under-determined
    >>> r, p = balance_stoichiometry({'Fe', 'O2'}, {'FeO', 'Fe2O3'})
    >>> list(set.union(*[v.free_symbols for v in r.values()]))
    [x1]
    >>> b = balance_stoichiometry({'Fe', 'O2'}, {'FeO', 'Fe2O3'}, underdetermined=None)
    >>> b == ({'Fe': 3, 'O2': 2}, {'FeO': 1, 'Fe2O3': 1})
    True

    Returns
    -------
    balanced reactants : dict
    balanced products : dict

    """
    import sympy
    from sympy import (
        MutableDenseMatrix, gcd, zeros, linsolve, numbered_symbols, Wild, Symbol,
        Integer, Tuple, preorder_traversal as pre
    )

    _intersect = set.intersection(*map(set, (reactants, products)))
    if _intersect:
        raise ValueError("Substances on both sides: %s" % str(_intersect))
    if substances is None:
        substances = OrderedDict([(k, substance_factory(k)) for k in chain(reactants, products)])
    if isinstance(substances, str):
        substances = OrderedDict([(k, substance_factory(k)) for k in substances.split()])
    if type(reactants) == set:  # we don't want isinstance since it might be "OrderedSet"
        reactants = sorted(reactants)
    if type(products) == set:
        products = sorted(products)
    subst_keys = list(reactants) + list(products)

    cks = Substance.composition_keys(substances.values())

    if parametric_symbols is None:
        parametric_symbols = numbered_symbols('x', start=1, integer=True, positive=True)

    # ?C2H2 + ?O2 -> ?CO + ?H2O
    # Ax = 0
    #   A:                    x:
    #
    #   C2H2   O2  CO  H2O
    # C -2     0    1   0      x0    =   0
    # H -2     0    0   2      x1        0
    # O  0    -2    1   1      x2        0
    #                          x3

    def _get(ck, sk):
        return substances[sk].composition.get(ck, 0) * (-1 if sk in reactants else 1)

    for ck in cks:  # check that all components are present on reactant & product sides
        for rk in reactants:
            if substances[rk].composition.get(ck, 0) != 0:
                break
        else:
            raise ValueError("Component '%s' not among reactants" % ck)
        for pk in products:
            if substances[pk].composition.get(ck, 0) != 0:
                break
        else:
            raise ValueError("Component '%s' not among products" % ck)

    A = MutableDenseMatrix([[_get(ck, sk) for sk in subst_keys] for ck in cks])
    symbs = list(reversed([next(parametric_symbols) for _ in range(len(subst_keys))]))
    sol, = linsolve((A, zeros(len(cks), 1)), symbs)
    wi = Wild('wi', properties=[lambda k: not k.has(Symbol)])
    cd = reduce(gcd, [1] + [1/m[wi] for m in map(
        lambda n: n.match(symbs[-1]/wi), pre(sol)) if m is not None])
    sol = sol.func(*[arg/cd for arg in sol.args])

    def remove(cont, symb, remaining):
        subsd = dict(zip(remaining/symb, remaining))
        cont = cont.func(*[(arg/symb).expand().subs(subsd) for arg in cont.args])
        if cont.has(symb):
            raise ValueError("Bug, please report an issue at https://github.com/bjodah/chempy")
        return cont

    done = False
    for idx, symb in enumerate(symbs):
        for expr in sol:
            iterable = expr.args if expr.is_Add else [expr]
            for term in iterable:
                if term.is_number:
                    done = True
                    break
            if done:
                break
        if done:
            break
        for expr in sol:
            if (expr/symb).is_number:
                sol = remove(sol, symb, MutableDenseMatrix(symbs[idx+1:]))
                break
    for symb in symbs:
        cd = 1
        for expr in sol:
            iterable = expr.args if expr.is_Add else [expr]
            for term in iterable:
                if term.is_Mul and term.args[0].is_number and term.args[1] == symb:
                    cd = gcd(cd, term.args[0])
        if cd != 1:
            sol = sol.func(*[arg.subs(symb, symb/cd) for arg in sol.args])
    if underdetermined is 1:
        from ._release import __version__
        if int(__version__.split('.')[1]) > 6:
            warnings.warn(  # deprecated because comparison with ``1`` problematic (True==1)
                ("Pass underdetermined == None instead of ``1`` (depreacted since 0.7.0,"
                 " will_be_missing_in='0.9.0')"), ChemPyDeprecationWarning)
        underdetermined = None
    if underdetermined is None:
        sol = Tuple(*[Integer(x) for x in _solve_balancing_ilp_pulp(A)])

    fact = gcd(sol)
    sol = MutableDenseMatrix([e/fact for e in sol]).reshape(len(sol), 1)
    sol /= reduce(gcd, sol)
    if 0 in sol:
        raise ValueError("Superfluous species given.")
    if underdetermined:
        if any(x == sympy.nan for x in sol):
            raise ValueError("Failed to balance reaction")
    else:
        for x in sol:
            if len(x.free_symbols) != 0:
                raise ValueError("The system was under-determined")
        if not all(residual == 0 for residual in A * sol):
            raise ValueError("Failed to balance reaction")

    def _x(k):
        coeff = sol[subst_keys.index(k)]
        return int(coeff) if underdetermined is None else coeff

    return (
        OrderedDict([(k, _x(k)) for k in reactants]),
        OrderedDict([(k, _x(k)) for k in products])
    )


def mass_fractions(stoichiometries, substances=None, substance_factory=Substance.from_formula):
    """ Calculates weight fractions of each substance in a stoichiometric dict

    Parameters
    ----------
    stoichiometries : dict or set
        If a ``set``: all entries are assumed to correspond to unit multiplicity.
    substances: dict or None

    Examples
    --------
    >>> r = mass_fractions({'H2': 1, 'O2': 1})
    >>> mH2, mO2 = 1.008*2, 15.999*2
    >>> abs(r['H2'] - mH2/(mH2+mO2)) < 1e-4
    True
    >>> abs(r['O2'] - mO2/(mH2+mO2)) < 1e-4
    True
    >>> mass_fractions({'H2O2'}) == {'H2O2': 1.0}
    True

    """
    if isinstance(stoichiometries, set):
        stoichiometries = {k: 1 for k in stoichiometries}
    if substances is None:
        substances = OrderedDict([(k, substance_factory(k)) for k in stoichiometries])
    tot_mass = sum([substances[k].mass*v for k, v in stoichiometries.items()])
    return {k: substances[k].mass*v/tot_mass for k, v in stoichiometries.items()}
