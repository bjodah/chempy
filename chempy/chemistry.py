# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict, defaultdict
from functools import reduce
from itertools import chain
from operator import itemgetter, mul
import math
import sys
import warnings

from .util.arithmeticdict import ArithmeticDict
from .util.parsing import (
    formula_to_composition, mass_from_composition, to_reaction,
    formula_to_latex, formula_to_unicode, formula_to_html
)

from .units import to_unitless, default_units
from ._util import intdiv
from .util.pyutil import deprecated, ChemPyDeprecationWarning


class Substance(object):
    """ Class representing a chemical substance

    Parameters
    ----------
    name : str
    charge : int (optional, default: None)
        will be stored in composition[0], prefer composition when possible
    latex_name : str
    unicode_name : str
    html_name : str
    composition : dict or None (default)
        dict (int -> number) e.g. {atomic number: count}, zero has special
        meaning (net charge)
    data : dict
        free form dictionary. Could be simple such as ``{'mp': 0, 'bp': 100}``
        or considerably more involved, e.g.: ``{'diffusion_coefficient': {\
 'water': lambda T: 2.1*m**2/s/K*(T - 273.15*K)}}``

    Attributes
    ----------
    mass
        maps to data['mass'], and when unavailable looks for formula.mass
    attrs
        a tuple of attribute names for serialization
    data
        free form dict
    other_properties
        deprecated alias to :attr:`data`

    Examples
    --------
    >>> ammonium = Substance('NH4+', 1, 'NH_4^+', composition={7: 1, 1: 4},
    ...     data={'mass': 18.0385, 'pKa': 9.24})
    >>> ammonium.name
    'NH4+'
    >>> ammonium.composition  # note that charge was inserted as composition[0]
    {0: 1, 1: 4, 7: 1}
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

    @property
    @deprecated(will_be_missing_in='0.5.0')
    def other_properties(self):
        return self.data

    @other_properties.setter
    @deprecated(will_be_missing_in='0.5.0')
    def other_properties(self, value):
        self.data = value

    def molar_mass(self, units=None):
        """ Returns the molar mass (with units) of the substance

        Examples
        --------
        >>> nh4p = Substance.from_formula('NH4+')  # simpler
        >>> from chempy.units import default_units as u
        >>> nh4p.molar_mass(u)
        array(18.0384511) * g/mol

        """
        if units is None:
            units = default_units
        return self.mass*units.g/units.mol

    def __init__(self, name=None, charge=None, latex_name=None, unicode_name=None,
                 html_name=None, composition=None, data=None, other_properties=None):
        if other_properties is not None:  # will_be_missing_in='0.5.0'
            if data is not None:
                raise ValueError("Cannot take both data and other_properties")
            else:
                data = other_properties
                warnings.warn("use data kwarg instead of other_properties", ChemPyDeprecationWarning)
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
        self.data = data or {}

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

        2 R --> A + P; r = k*A*R*R

    this would be represented as ``Reaction({'A': 1, 'R': 2},
    {'A': 2, 'P': 1}, param=k)``. Some reactions have a larger
    stoichiometric coefficient than what appears in the rate
    expression, e.g.:

        5 A + B --> C; r = k*A*B

    this can be represented as ``Reaction({'C1': 1, 'C2': 1},
    {'B': 1}, inact_reac={'C1': 4}, param=k)``.

    The rate constant information in ``param`` may be a subclass of
    :class:`chempy.kinetics.rates.RateExpr` or carry a :meth:`_as_RateExpr`,
    if neither: `param` will be assumed to be a rate constant for a mass-action
    type of kinetic expression.

    Additional data may be stored in the ``data`` dict.


    Parameters
    ----------
    reac : dict (str -> int)
        if reac is a set multiplicities are assumed to be 1
    prod : dict (str -> int)
        if reac is a set multiplicities are assumed to be 1
    param : float or callable
        Special case (side-effect): if param is a subclass of
        :class:`.kinetics.rates.RateExpr` and its :attr:`rxn`
        is `None` it will be set to `self`.
    inact_reac : dict (optional)
    inact_prod : dict (optional)
    name : str (optional)
    k : deprecated (alias for param)
    ref : object
        reference
    data : dict (optional)
    checks : iterable of str
        raises value error if any method check_%s returns False
        for all %s in checks.

    Attributes
    ----------
    reac: dict
    prod: dict
    param: object
    inact_reac: dict
    inact_prod: dict
    name: str
    ref: str
    data: dict

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

    str_arrow = '->'
    latex_arrow = r'\rightarrow'
    unicode_arrow = u'→'
    html_arrow = '&rarr;'
    param_char = 'k'  # convention

    def __init__(
            self, reac, prod, param=None, inact_reac=None, inact_prod=None,
            name=None, ref=None, data=None,
            checks=('any_effect', 'all_positive', 'all_integral')):
        if isinstance(reac, set):
            reac = {k: 1 for k in reac}
        if isinstance(inact_reac, set):
            inact_reac = {k: 1 for k in inact_reac}
        if isinstance(prod, set):
            prod = {k: 1 for k in prod}
        if isinstance(inact_prod, set):
            inact_prod = {k: 1 for k in inact_prod}
        self.reac = reac
        self.prod = prod
        self.param = param
        self.inact_reac = inact_reac or {}
        self.inact_prod = inact_prod or {}
        self.name = name
        self.ref = ref
        self.data = data or {}

        from .kinetics.rates import RateExpr
        if isinstance(self.param, RateExpr) and self.param.rxn is None:
            self.param.rxn = self  # register instance in rate expression

        for check in checks:
            if not getattr(self, 'check_'+check)():
                raise ValueError("Check failed %s" % check)

    @classmethod
    def from_string(cls, string, substance_keys=None, globals_=None, **kwargs):
        """ Parses a string into a Reaction instance

        Parameters
        ----------
        string : str
            String representation of the reaction.
        substance_keys : iterable of strings or string or None
            Used prevent e.g. misspelling.
            if str: split is invoked, if None: no checking done.
        globals_ : dict (optional)
            Dictionary for eval for (default: None -> {'chempy': chempy})
            If ``False``: no eval will be called (useful for web-apps).
        \*\*kwargs :
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


        Notes
        -----
        :func:`chempy.util.parsing.to_reaction` is used which in turn calls
        :func:`eval` which is a severe security concern for untrusted input.

        """
        if isinstance(substance_keys, str):
            if ' ' in substance_keys:
                substance_keys = substance_keys.split()
        return to_reaction(string, substance_keys, cls.str_arrow, cls, globals_, **kwargs)

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

    def all_prod_stoich(self, substances):
        """ Per substance product stoichiometry tuple (active & inactive) """
        return tuple(self.prod.get(k, 0) + self.inact_prod.get(k, 0) for k in substances)

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
        for s_name in chain(self.reac.keys(), self.prod.keys(), self.inact_reac.keys(), self.inact_prod.keys()):
            if substances[s_name].precipitate:
                return True
        return False

    def _get_str_parts(self, name_attr, arrow_attr, substances, _str=str):
        def not_None(arg, default):
            if arg is None:
                return default
            return arg
        nullstr, space = _str(''), _str(' ')
        reac, prod, i_reac, i_prod = [[
            ((_str(v)+space) if v > 1 else nullstr) + _str(not_None(getattr(substances[k], name_attr, k), k))
            for k, v in filter(itemgetter(1), d.items())
        ] for d in (self.reac, self.prod, self.inact_reac, self.inact_prod)]
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

    def _str_param(self, magnitude_fmt=lambda x: '%.3g' % x, unit_fmt=str, _str=str):
        try:
            magnitude_str = magnitude_fmt(self.param.magnitude)
            unit_str = unit_fmt(self.param.dimensionality)
        except AttributeError:
            try:
                return magnitude_fmt(self.param)
            except TypeError:
                return str(self.param)
        else:
            return magnitude_str + _str(' ') + unit_str

    def string(self, substances=None, with_param=False, magnitude_fmt=lambda m: r'\num{%s}' % m):
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
        if substances is None:
            substances = {
                k: k for k in chain(self.reac.keys(), self.prod.keys(),
                                    self.inact_reac.keys(),
                                    self.inact_prod.keys())
            }
        res = self._get_str('name', 'str_arrow', substances)
        if with_param and self.param is not None:
            res += '; ' + self._str_param()
        return res

    def __str__(self):
        return self.string(with_param=True)

    def latex(self, substances, with_param=False, magnitude_fmt=lambda m: r'\num{%s}' % m):
        r""" Returns a LaTeX representation of the reaction

        Parameters
        ----------
        substances: dict
            mapping substance keys to Substance instances
        with_param: bool
            whether to print the parameter (default: False)
        magnitude_fmt: callback
            how to format the number (default requires siunitx's num)

        Examples
        --------
        >>> keys = 'H2O H+ OH-'.split()
        >>> subst = {k: Substance.from_formula(k) for k in keys}
        >>> r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", subst)
        >>> r.latex(subst) == r'H_{2}O \rightarrow H^{+} + OH^{-}'
        True
        >>> r2 = Reaction.from_string("H2O -> H+ + OH-; 1e-8/molar/second", subst)
        >>> ref = r'H_{2}O \rightarrow H^{+} + OH^{-}; 1\cdot 10^{-8} $\mathrm{\frac{1}{(s{\cdot}M)}}$'
        >>> r2.latex(subst, with_param=True) == ref
        True

        """
        res = self._get_str('latex_name', 'latex_arrow', substances)
        if with_param and self.param is not None:
            from .util.parsing import number_to_scientific_latex as _fmt
            res += '; %s' % self._str_param(magnitude_fmt=_fmt, unit_fmt=lambda dim: dim.latex)
        return res

    def unicode(self, substances, with_param=False):
        u""" Returns a unicode string representation of the reaction

        Examples
        --------
        >>> keys = 'H2O H+ OH-'.split()
        >>> subst = {k: Substance.from_formula(k) for k in keys}
        >>> r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", subst)
        >>> r.unicode(subst) == u'H₂O → H⁺ + OH⁻'
        True
        >>> r2 = Reaction.from_string("H2O -> H+ + OH-; 1e-8/molar/second", subst)
        >>> r2.unicode(subst, with_param=True) == u'H₂O → H⁺ + OH⁻; 1·10⁻⁸ 1/(s·M)'
        True

        """
        res = self._get_str('unicode_name', 'unicode_arrow', substances,
                            _str=str if sys.version_info[0] > 2 else unicode)
        if with_param and self.param is not None:
            from .util.parsing import number_to_scientific_unicode
            res += u'; ' + self._str_param(
                magnitude_fmt=number_to_scientific_unicode,
                unit_fmt=lambda dim: (
                    dim.unicode if sys.version_info[0] > 2
                    else dim.unicode.decode(encoding='utf-8')
                ), _str=str if sys.version_info[0] > 2 else unicode)
        return res

    def html(self, substances, with_param=False):
        """ Returns a HTML representation of the reaction

        Examples
        --------
        >>> keys = 'H2O H+ OH-'.split()
        >>> subst = {k: Substance.from_formula(k) for k in keys}
        >>> r = Reaction.from_string("H2O -> H+ + OH-; 1e-4", subst)
        >>> r.html(subst)
        'H<sub>2</sub>O &rarr; H<sup>+</sup> + OH<sup>-</sup>'
        >>> r2 = Reaction.from_string("H2O -> H+ + OH-; 1e-8/molar/second", subst)
        >>> r2.html(subst, with_param=True)
        'H<sub>2</sub>O &rarr; H<sup>+</sup> + OH<sup>-</sup>&#59; 1&sdot;10<sup>-8</sup> 1/(s*M)'

        """
        res = self._get_str('html_name', 'html_arrow', substances)
        if with_param and self.param is not None:
            from .util.parsing import number_to_scientific_html as _fmt
            res += '&#59; '
            try:
                res += self.param.string(_fmt)
            except AttributeError:
                res += self._str_param(magnitude_fmt=_fmt)
        return res

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
        """
        if composition_keys is None:
            composition_keys = Substance.composition_keys(substances.values())
        net = [0]*len(composition_keys)
        for substance, coeff in zip(substances.values(),
                                    self.net_stoich(substances.keys())):
            for idx, key in enumerate(composition_keys):
                net[idx] += substance.composition.get(key, 0) * coeff
        return net

    def rate_expr(self):
        """ Turns self.param into a RateExpr instance (if not already)

        Examples
        --------
        >>> r = Reaction.from_string('2 A + B -> 3 C; 7')
        >>> ratex = r.rate_expr()
        >>> ratex.args[0] == 7
        True

        """
        from chempy.kinetics.rates import RateExpr, MassAction
        if isinstance(self.param, RateExpr):
            if self.param.rxn is None:
                self.param.rxn = self
            return self.param
        else:
            try:
                convertible = self.param._as_RateExpr
            except AttributeError:
                return MassAction([self.param], rxn=self)
            else:
                return convertible(self)

    def rate(self, variables=None, backend=math, substance_keys=None):
        """ Evaluate the rate of a reaction

        Parameters
        ----------
        variables: dict
        backend: module, optional
        substance_keys: iterable of str, optional

        Examples
        --------
        >>> rxn = Reaction.from_string('2 H2 + O2 -> 2 H2O; 3', None)
        >>> r = 3*5*5*7
        >>> rxn.rate({'H2': 5, 'O2': 7}) == {'H2': -2*r, 'O2': -r, 'H2O': 2*r}
        True

        """
        if variables is None:
            variables = {}
        if substance_keys is None:
            substance_keys = self.keys()
        r = self.rate_expr()(variables, backend)
        return {k: r*v for k, v in zip(substance_keys, self.net_stoich(substance_keys))}


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

    def as_reactions(self, kf=None, kb=None, units=None, variables=None, backend=math):
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
            if kb is None:
                raise ValueError("Exactly one rate needs to be provided")
            kf = kb * self.equilibrium_constant(variables, backend=backend) * c0**(nb - nf)
        elif kb is None:
            kb = kf / (self.equilibrium_constant(variables, backend=backend) * c0**(nb - nf))
        else:
            raise ValueError("Exactly one rate needs to be provided")
        return (
            Reaction(self.reac, self.prod, kf, self.inact_reac,
                     self.inact_prod, ref=self.ref),
            Reaction(self.prod, self.reac, kb, self.inact_prod,
                     self.inact_reac, ref=self.ref)
        )

    def equilibrium_constant(self, variables=None, backend=math):
        """ Return equilibrium constant

        Parameters
        ----------
        variables : dict, optional
        backend : module, optional

        """
        try:
            return self.param(variables, backend=backend)
        except TypeError:
            return self.param

    K = equilibrium_constant

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
    rxns : sequence
         sequence of :py:class:`Reaction` instances
    substances : OrderedDict or string or None
         mapping str -> Substance instances, None => deduced from reactions.
    name : string (optional)
         Name of ReactionSystem (e.g. model name / citation key)
    checks : iterable of str, optional
        raises value error if any method check_%s returns False
        for all %s in checks.
    substance_factory : callback
        e.g. :meth:`Substance.from_formula`

    Attributes
    ----------
    rxns : list of objects
        sequence of :class:`Reaction` instances
    substances : OrderedDict or string or iterable of strings/Substance
        mapping substance name to substance index
    ns : int
        number of substances
    nr : int
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

    _BaseReaction = Reaction
    _BaseSubstance = Substance

    def __init__(self, rxns, substances=None, name=None, checks=('balance', 'substance_keys',
                                                                 'duplicate', 'duplicate_names'),
                 substance_factory=Substance):
        self.rxns = list(rxns)
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

        self.name = name

        for check in checks:
            getattr(self, 'check_'+check)(throw=True)

    def _repr_html_(self):
        def _format(r):
            return r.html(self.substances, with_param=True)
        return '<br>'.join(map(_format, self.rxns))

    def check_duplicate(self, throw=False):
        """ Raies ValueError if there are duplicates in self.rxns """
        for i1, rxn1 in enumerate(self.rxns):
            for i2, rxn2 in enumerate(self.rxns[i1+1:], i1+1):
                if rxn1 == rxn2:
                    if throw:
                        raise ValueError("Duplicate reactions %d & %d" % (i1, i2))
                    else:
                        return False
        return True

    def check_duplicate_names(self, throw=False):
        names_seen = {}
        for idx, rxn in enumerate(self.rxns):
            if rxn.name is None:
                continue
            if rxn.name in names_seen:
                if throw:
                    raise ValueError("Duplicate names at %d: %s" % (idx, rxn.name))
                else:
                    return False
            else:
                names_seen[rxn.name] = idx
        return True

    def check_balance(self, strict=False, throw=False):
        """ Raies ValueError there are unbalanecd reactions in self.rxns """
        for subst in self.substances.values():
            if subst.composition is None:
                if strict:
                    if throw:
                        raise ValueError("No composition for %s" % str(subst))
                    else:
                        return False
                else:
                    return True
        for rxn in self.rxns:
            for net in rxn.composition_violation(self.substances):
                if net != 0:
                    if throw:
                        raise ValueError("Composition violation in %s" % str(rxn))
                    else:
                        return False
        return True

    def check_substance_keys(self, throw=False):
        for rxn in self.rxns:
            for key in chain(rxn.reac, rxn.prod, rxn.inact_reac,
                             rxn.inact_prod):
                if key not in self.substances:
                    if throw:
                        raise ValueError("Unknown key: %s" % key)
                    else:
                        return False
        return True

    @classmethod
    def from_string(cls, s, **kwargs):
        """ Create a reaction system from a string

        Parameters
        ----------
        s : str
            multiline string

        Examples
        --------
        >>> rs = ReactionSystem.from_string('\\n'.join(['2 HNO2 -> H2O + NO + NO2; 3', '2 NO2 -> N2O4; 4']))
        >>> r1, r2 = 5*5*3, 7*7*4
        >>> rs.rates({'HNO2': 5, 'NO2': 7}) == {'HNO2': -2*r1, 'H2O': r1, 'NO': r1, 'NO2': r1 - 2*r2, 'N2O4': r2}
        True

        """
        rxns = [cls._BaseReaction.from_string(r) for r in s.split('\n')]
        return cls(rxns, substance_factory=cls._BaseSubstance.from_formula, **kwargs)

    def __getitem__(self, key):
        candidate = None
        for r in self.rxns:
            if r.name == key:
                if candidate is None:
                    candidate = r
                else:
                    raise ValueError('Multiple reactions with the same name')
        if candidate is None:
            raise KeyError("No reaction with name %s found" % key)
        return candidate

    def __iadd__(self, other):
        try:
            self.substances.update(other.substances)
        except AttributeError:
            self.rxns.extend(other)
        else:
            self.rxns.extend(other.rxns)
        return self

    def __add__(self, other):
        try:
            substances = list(chain(self.substances.items(), other.substances.items()))
        except AttributeError:
            substances = self.substances.copy()
        return self.__class__(chain(self.rxns, getattr(other, 'rxns', other)), substances, checks=())

    def __eq__(self, other):
        if self is other:
            return True
        return self.rxns == other.rxns and self.substances == other.substances

    def substance_names(self):
        """ Returns a tuple of the substances' names """
        return tuple(substance.name for substance in self.substances.values())

    def substance_participation(self, substance_key):
        """ Returns indices of reactions where substance_key occurs

        Parameters
        ----------
        substance_key: str

        Returns
        -------
        List of indices for self.rxns where `substance_key` participates

        """
        return [ri for ri, rxn in enumerate(self.rxns) if substance_key in rxn.keys()]

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

    def as_per_substance_array(self, cont, dtype='float64', unit=None, raise_on_unk=False):
        """ Turns a dict into an ordered array

        Parameters
        ----------
        cont : array_like or dict
        dtype : str or numpy.dtype object
        unit : unit, optional
        raise_on_unk : bool

        """
        import numpy as np
        if unit is not None:
            cont = to_unitless(cont, unit)
        if isinstance(cont, np.ndarray):
            pass
        elif isinstance(cont, dict):
            substance_keys = self.substances.keys()
            if raise_on_unk:
                for k in cont:
                    if k not in substance_keys:
                        raise KeyError("Unkown substance key: %s" % k)
            cont = [cont[k] for k in substance_keys]

        cont = np.atleast_1d(np.asarray(cont, dtype=dtype).squeeze())
        if cont.shape[-1] != self.ns:
            raise ValueError("Incorrect size")
        return cont*(unit if unit is not None else 1)

    def as_per_substance_dict(self, arr):
        return dict(zip(self.substances.keys(), arr))

    def as_substance_index(self, substance_key):
        """ Returns the index of a Substance in the system"""
        if isinstance(substance_key, int):
            return substance_key
        else:
            return list(self.substances.keys()).index(substance_key)

    def per_substance_varied(self, per_substance, varied=None):
        """ Dense nd-array for all combinations of varied levels per substance

        Parameters
        ----------
        per_substance: dict or array
        varied: dict

        Examples
        --------
        >>> rsys = ReactionSystem([], 'A B C')
        >>> arr, keys = rsys.per_substance_varied({'A': 2, 'B': 3, 'C': 5}, {'C': [5, 7, 9, 11]})
        >>> arr.shape, keys
        ((4, 3), ('C',))
        >>> all(arr[1, :] == [2, 3, 7])
        True

        Returns
        -------
        ndarray : with len(varied) + 1 number of axes, and with last axis length == self.ns

        """
        import numpy as np
        varied = varied or {}
        varied_keys = tuple(k for k in self.substances if k in varied)
        n_varied = len(varied)
        shape = tuple(len(varied[k]) for k in self.substances if k in varied)
        result = np.empty(shape + (self.ns,))
        result[..., :] = self.as_per_substance_array(per_substance)
        if varied:
            for k, vals in varied.items():
                varied_axis = varied_keys.index(k)
                for varied_idx, val in enumerate(vals):
                    index = tuple(varied_idx if i == varied_axis else slice(None) for i in range(n_varied))
                    result[index + (self.as_substance_index(k),)] = val
        return result, varied_keys

    def rates(self, variables=None, backend=math, substance_keys=None):
        """ Per substance sums of reaction rates rates.

        Parameters
        ----------
        variables : dict
        backend : module, optional
        substance_keys : iterable of str, optional

        Returns
        -------
        dict
            per substance_key time derivatives of concentrations.

        Examples
        --------
        >>> r = Reaction({'R': 2}, {'P': 1}, 42.0)
        >>> rsys = ReactionSystem([r])
        >>> rates = rsys.rates({'R': 3, 'P': 5})
        >>> abs(rates['P'] - 42*3**2) < 1e-14
        True

        """
        result = {}
        for rxn in self.rxns:
            for k, v in rxn.rate(variables, backend, substance_keys).items():
                if k not in result:
                    result[k] = v
                else:
                    result[k] += v
        return result

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
        ) for idx, eq in enumerate(self.rxns)], dtype=object)

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

    def composition_balance_vectors(self):
        subs = self.substances.values()
        ck = Substance.composition_keys(subs)
        return [[s.composition.get(k, 0) for s in subs] for k in ck], ck

    def _html_table_cell_factory(self, title=True):
        if title:
            def _fmt(r):
                return '<a title="%s">%s</a>' % (r.unicode(self.substances, with_param=True), r.name)
        else:
            def _fmt(r):
                return r.name
        missing = [len(self.substance_participation(k)) == 0 for k in self.substances]

        def cell(A, ri, ci=None):
            args = []
            if ci is not None and ri > ci:
                r = '-'
            else:
                if ci is None:
                    c = A[ri]
                    color_red = missing[ri]
                else:
                    c = A[ri][ci]
                    color_red = missing[ri] or missing[ci]
                if c is None:
                    r = ''
                else:
                    r = ', '.join(_fmt(r) for r in c)
                if color_red:
                    args.append('style="background-color:#faa"')
            return '<td %s>%s</td>' % (' '.join(args), r)
        return cell

    def _unimolecular_reactions(self):
        A = [None]*self.ns
        unconsidered = []
        for r in self.rxns:
            if r.order() == 1:
                keys = list(r.reac.keys())
                if len(keys) == 1:
                    ri = self.as_substance_index(keys[0])
                else:
                    raise NotImplementedError("Need 1 or 2 keys")
                if A[ri] is None:
                    A[ri] = list()
                A[ri].append(r)
            else:
                unconsidered.append(r)
        return A, unconsidered

    def unimolecular_html_table(self, title=True):
        """ Returns a HTML table of unimolecular reactions

        Parameters
        ----------
        title: bool

        Returns
        -------
        string: html representation
        list: reactions not considered
        """
        A, unconsidered = self._unimolecular_reactions()
        _cell = self._html_table_cell_factory(title)
        rows = '\n'.join('<tr><td>%s</td>%s</tr>' % (
            (s.html_name or s.name), _cell(A, ri)
        ) for ri, s in enumerate(self.substances.values()))
        html = '<table>%s</table>' % rows
        return html, unconsidered

    def _bimolecular_reactions(self):
        A = [[None]*self.ns for _ in range(self.ns)]
        unconsidered = []
        for r in self.rxns:
            if r.order() == 2:
                keys = list(r.reac.keys())
                if len(keys) == 1:
                    ri = ci = self.as_substance_index(keys[0])
                elif len(keys) == 2:
                    ri, ci = sorted(map(self.as_substance_index, keys))
                else:
                    raise NotImplementedError("Need 1 or 2 keys")
                if A[ri][ci] is None:
                    A[ri][ci] = list()
                A[ri][ci].append(r)
            else:
                unconsidered.append(r)
        return A, unconsidered

    def bimolecular_html_table(self, title=True):
        """ Returns a HTML table of bimolecular reactions

        Parameters
        ----------
        title: bool

        Returns
        -------
        string: html representation
        list: reactions not considered
        """
        A, unconsidered = self._bimolecular_reactions()
        header = '<th></th>' + ''.join('<th>%s</th>' % (s.html_name or s.name) for s in self.substances.values())
        _cell = self._html_table_cell_factory(title)
        rows = '\n'.join('<tr><td>%s</td>%s</tr>' % (
            (s.html_name or s.name), ''.join(_cell(A, ri, ci) for ci in range(self.ns))
        ) for ri, s in enumerate(self.substances.values()))
        html = '<table>%s</table>' % '\n'.join([header, rows])
        return html, unconsidered


def balance_stoichiometry(reactants, products, substances=None,
                          substance_factory=Substance.from_formula):
    """ Balances stoichiometric coefficients of a reaction

    Parameters
    ----------
    reactants: iterable of reactant keys
    products: iterable of product keys
    substances: OrderedDict or string or None
        mapping reactant/product keys to instances of :class:`Substance`
    substance_factory: callback

    Examples
    --------
    >>> ref = {'C2H2': 2, 'O2': 3}, {'CO': 4, 'H2O': 2}
    >>> balance_stoichiometry({'C2H2', 'O2'}, {'CO', 'H2O'}) == ref
    True
    >>> ref2 = {'H2': 1, 'O2': 1}, {'H2O2': 1}
    >>> balance_stoichiometry('H2 O2'.split(), ['H2O2'], 'H2 O2 H2O2') == ref2
    True


    Returns
    -------
    balanced reactants : dict
    balanced products : dict

    """
    from sympy import Matrix

    _intersect = set.intersection(*map(set, (reactants, products)))
    if _intersect:
        raise ValueError("Substances on both sides: %s" % str(_intersect))
    if substances is None:
        substances = OrderedDict([(k, substance_factory(k)) for k
                                  in chain(reactants, products)])
    if isinstance(substances, str):
        substances = OrderedDict([(k, substance_factory(k)) for k
                                  in substances.split()])
    subst_keys = list(substances.keys())

    cks = Substance.composition_keys(substances.values())
    nsubs = len(substances)

    # ?C2H2 + ?O2 -> ?CO + ?H2O
    # Ax = 0
    # A:                 x:
    #
    #   C2H2   O2  CO  H2O
    # C  2     0    1   0    x0
    # H  2     0    0   2    x1
    # O  0    -2    1   1    x2

    def _get(sk, ck):
        return substances[sk].composition.get(ck, 0) * (-1 if sk in reactants else 1)

    A = Matrix([[_get(sk, ck) for sk in subst_keys] for ck in cks])

    # A2 x = b
    #
    # A2:                 x:   b:
    #
    #      O2  CO  H2O       C2H2
    # C    0    1   0    x0   2
    # H    0    0   2    x1   2
    # O   -2    1   1    x2   0

    A_aug, pivot = A.rref()
    if len(pivot) < nsubs-1:
        raise ValueError("Unsatisfiable system of equations")
    x_aug = Matrix(A_aug[:len(pivot), 1:]).LUsolve(Matrix(-A_aug[:len(pivot), 0]))

    # Reorder to original indices
    x = [1]
    for si in range(1, nsubs):
        ai = si - 1  # augmented index
        if ai in pivot:
            x.append(x_aug[pivot.index(ai)])
        else:
            x.append(None)

    # Now solve for the redundant x:s
    for si in range(1, nsubs):
        elem = x[si]
        if elem is None:
            # solve
            col = A[:, si]
            for ri, cell in enumerate(col):
                if cell == 0:
                    continue
                others = 0
                for ci, comp in enumerate(A[ri, :]):
                    if ci == si:
                        continue
                    if x[ci] is None:
                        raise NotImplementedError("Need a second LU solve")
                    others += comp*x[ci]
                x[si] = -others/cell
                break

    x = Matrix(x)
    while True:
        for idx in range(nsubs):
            elem = x[idx]
            if not elem.is_integer:
                numer, denom = elem.as_numer_denom()
                x *= denom
                break
        else:
            break
    if 0 in x:
        raise ValueError("Unable to balance stoichiometry (did you forget a product?)")

    def _x(k):
        return x[subst_keys.index(k)]

    return (
        {k: _x(k) for k in reactants},
        {k: _x(k) for k in products}
    )


def mass_fractions(stoichiometries, substances=None, substance_factory=Substance.from_formula):
    """ Calculates weight fractions of each substance in a stoichiometric dict

    Parameters
    ----------
    stoichiometries: dict or set
        if a set: all entries are assumed to correspond to unit multiplicity
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
