# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import math

from collections import OrderedDict, defaultdict
from itertools import chain
import warnings

from .chemistry import Reaction, Substance
from .units import to_unitless


class ReactionSystem(object):
    """ Collection of reactions forming a system (model).

    Parameters
    ----------
    rxns : sequence
         Sequence of :py:class:`Reaction` instances.
    substances : OrderedDict or string or None
         Mapping str -> Substance instances, None => deduced from reactions.
         If a set is passed as substances (or a string which is split), then
         ``substance_factory`` will be used to construct substances from the items.
    name : string (optional)
         Name of ReactionSystem (e.g. model name / citation key).
    checks : iterable of str, optional
        Raises ``ValueError`` if any method ``check_%s`` returns False
        for all ``%s`` in ``checks``.
    substance_factory : callback
        Could also be e.g. :meth:`Substance.from_formula`.

    Attributes
    ----------
    rxns : list of objects
        Sequence of :class:`Reaction` instances.
    substances : OrderedDict or string or iterable of strings/Substance
        Mapping substance name to substance index.
    ns : int
        Number of substances.
    nr : int
        Number of reactions.

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
                 substance_factory=Substance, missing_substances_from_keys=False):
        self.rxns = list(rxns)

        if substances is None:
            substances = set.union(*[set(rxn.keys()) for rxn in self.rxns])

        if isinstance(substances, OrderedDict):
            self.substances = substances
        elif isinstance(substances, (str, set)):
            if isinstance(substances, str) and ' ' in substances:
                substances = substances.split()
            self.substances = OrderedDict([
                (s, substance_factory(s)) for s in substances])
        else:
            if all(isinstance(s, Substance) for s in substances):
                self.substances = OrderedDict([(s.name, s) for s in substances])
            elif hasattr(substances, 'values') and all(isinstance(s, Substance) for s in substances.values()):
                self.substances = OrderedDict(substances)
            else:
                self.substances = OrderedDict((k, substance_factory(k)) for k in substances)

        if missing_substances_from_keys:
            for k in set.union(*[set(rxn.keys()) for rxn in self.rxns]) - set(self.substances):
                self.substances[k] = substance_factory(k)

        self.name = name

        for check in checks:
            getattr(self, 'check_'+check)(throw=True)

    def split(self, **kwargs):
        """ Splits the reaction system into multiple disjoint reaction systems. """
        groups = []  # tuples of (list, set) -- list of reactions, set of substance keys
        for i, r in enumerate(self.rxns):
            for gr, gs in groups: # check if reaction is part of group, break out
                rks = r.keys()
                group_found = False
                for k in rks:
                    if k in gs:
                        gr.append(i)
                        gs.update(rks)
                        group_found = True
                        break
                if group_found:
                    break
            else:  # reaction did not fit any group
                groups.append(([i], set(r.keys())))
        # We might have too many groups as this point, we will now recursively fuse groups
        i = 0
        while True:
            for j in range(i+1, len(groups)):
                if groups[i][1] & groups[j][1]:  # do groups share a substance?
                    groups[i][0].extend(groups[j][0])
                    groups[i][1].update(groups[j][1])
                    groups.pop(j)
                    break
            else:
                i += 1
            if i >= len(groups):
                break
        return [self.__class__(
            [self.rxns[i] for i in gr],
            OrderedDict([(k, v) for k, v in self.substances.items() if k in gs]),
            **kwargs
        ) for gr, gs in groups]

    def sinks_sources_disjoint(self, **kwargs):
        """ Returns the keys for sink-, source- and disjoint-substances in the reaction system.

        Sinks are substances that are irreversibly formed (they don't appear on the
        reactant side of any reaction). Sources are eventually completely depleted (they do not
        appear on the product side of any reaction). Disjoint substances don't appear on either
        side and are thus unaffected by the reactionsystem.

        Returns
        -------
        Length-3 tuple of sets with substance keys (strings)

        """
        import numpy as np
        irrev_rxns = []
        for r in self.rxns:
            try:
                irrev_rxns.extend(r.as_reactions())
            except AttributeError:
                irrev_rxns.append(r)
        irrev_rsys = ReactionSystem(irrev_rxns, self.substances, **kwargs)
        all_r = irrev_rsys.all_reac_stoichs()
        all_p = irrev_rsys.all_prod_stoichs()
        if np.any(all_r < 0) or np.any(all_p < 0):
            raise ValueError("Expected positive stoichiometric coefficients")
        result = sinks, sources, disjoint = set(), set(), set()
        for i, sk in enumerate(irrev_rsys.substances.keys()):
            in_r = np.any(all_r[:, i])
            in_p = np.any(all_p[:, i])
            if in_r and in_p:
                pass
            elif in_r:
                sources.add(sk)
            elif in_p:
                sinks.add(sk)
            else:
                disjoint.add(sk)
        return result

    def sort_substances_inplace(self, key=lambda kv: kv[0]):
        """ Sorts the OrderedDict attribute ``substances`` """
        self.substances = OrderedDict(sorted(self.substances.items(), key=key))

    def html(self, with_param=True, checks=(), color_sinks_sources=True):
        """ Returns a string with an HTML representation """
        parts = self.split(checks=checks)
        if len(parts) > 1:
            return '<br><br>'.join(rs.html(with_param) for rs in parts)

        sinks, sources, disjoint = self.sinks_sources_disjoint(checks=checks)
        def _str_formula(s, k):
            if k in sinks:
                fmt = '<span style="color:darkgreen">%s</span>'
            elif k in sources:
                fmt = '<span style="color:darkred">%s</span>'
            else:
                fmt = '%s'
            return fmt % s

        def _format(r):
            return r.html(self.substances, with_param=with_param,
                          str_formula=_str_formula if color_sinks_sources else None)
        return '<br>'.join(map(_format, self.rxns))

    def string(self, with_param=True):
        def _format(r):
            return r.string(self.substances, with_param=with_param)
        return '\n'.join(map(_format, self.rxns))

    def _repr_html_(self):
        return self.html()

    def check_duplicate(self, throw=False):
        """ Raies ValueError if there are duplicates in ``self.rxns`` """
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

    def check_balance(self, strict=False, throw=False):
        """ Checks if all reactions are balanced.

        Parameters
        ----------
        strict : bool
            Puts a requirement on all substances to have their ``composition`` attribute set.
        throw : bool
            Raies ValueError if there are unbalanecd reactions in self.rxns

        """
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

    @classmethod
    def from_string(cls, s, substances=None, rxn_parse_kwargs=None, **kwargs):
        """ Create a reaction system from a string

        Parameters
        ----------
        s : str
            Multiline string.
        substances : convertible to iterable of str
        rxn_parse_kwargs : dict
            Keyword arguments passed on to the Reaction baseclass' method ``from_string``.
        \\*\\*kwargs:
            Keyword arguments passed to the constructor of the class

        Examples
        --------
        >>> rs = ReactionSystem.from_string('\\n'.join(['2 HNO2 -> H2O + NO + NO2; 3', '2 NO2 -> N2O4; 4']))
        >>> r1, r2 = 5*5*3, 7*7*4
        >>> rs.rates({'HNO2': 5, 'NO2': 7}) == {'HNO2': -2*r1, 'H2O': r1, 'NO': r1, 'NO2': r1 - 2*r2, 'N2O4': r2}
        True

        """
        substance_keys = None if kwargs.get('missing_substances_from_keys', False) else substances
        rxns = [cls._BaseReaction.from_string(r, substance_keys, **(rxn_parse_kwargs or {}))
                for r in s.split('\n') if r.strip() != '']
        if 'substance_factory' not in kwargs:
            kwargs['substance_factory'] = cls._BaseSubstance.from_formula
        return cls(rxns, substances, **kwargs)

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
            substances = OrderedDict(chain(self.substances.items(), other.substances.items()))
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
        if isinstance(cont, np.ndarray):
            pass
        elif isinstance(cont, dict):
            substance_keys = self.substances.keys()
            if raise_on_unk:
                for k in cont:
                    if k not in substance_keys:
                        raise KeyError("Unkown substance key: %s" % k)
            cont = [cont[k] for k in substance_keys]
        if unit is not None:
            cont = to_unitless(cont, unit)

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

    def per_reaction_effect_on_substance(self, substance_key):
        result = {}
        for ri, rxn in enumerate(self.rxns):
            n, = rxn.net_stoich((substance_key,))
            if n != 0:
                result[ri] = n
        return result

    def rates(self, variables=None, backend=math, substance_keys=None, ratexs=None, cstr_fr_fc=None):
        """ Per substance sums of reaction rates rates.

        Parameters
        ----------
        variables : dict
        backend : module, optional
        substance_keys : iterable of str, optional
        ratexs : iterable of RateExpr instances
        cstr_fr_fc : tuple (str, tuple of str)
            Continuously stirred tank reactor conditions. Pair of
            flow/volume ratio key (feed-rate/tank-volume) and feed concentration
            keys. (if second item is a string it is taken to be a prefix)
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
        if ratexs is None:
            ratexs = [None]*self.nr
        for rxn, ratex in zip(self.rxns, ratexs):
            for k, v in rxn.rate(variables, backend, substance_keys, ratex=ratex).items():
                if k not in result:
                    result[k] = v
                else:
                    result[k] += v
        if cstr_fr_fc:
            if substance_keys is not None and tuple(substance_keys) != tuple(self.substances.keys()):
                warnings.warn("Only a subset of substances subject to CSTR treatment")
            substance_keys = (substance_keys or tuple(self.substances.keys()))

            fr_key, fc = cstr_fr_fc
            if isinstance(fc, str):
                fc_keys = [fc + k for k in substance_keys]
            elif isinstance(fc, dict):
                fc_keys = [fc[k] for k in substance_keys]
            else:
                fc_keys = fc
            if len(fc) != len(substance_keys):
                raise ValueError("Got incorrect number of feed concentration keys")
            fr = variables[fr_key]  # feed rate / tank volume ratio

            for fck, sk in zip(fc_keys, substance_keys):
                result[sk] += fr*(variables[fck] - variables[sk])
        return result

    def _stoichs(self, attr, keys=None):
        import numpy as np
        if keys is None:
            keys = self.substances.keys()
        # dtype: see https://github.com/sympy/sympy/issues/10295
        return np.array([(getattr(eq, attr)(keys)) for eq in self.rxns], dtype=object)

    def net_stoichs(self, keys=None):
        return self._stoichs('net_stoich', keys)

    def all_reac_stoichs(self, keys=None):
        return self._stoichs('all_reac_stoich', keys)

    def active_reac_stoichs(self, keys=None):
        return self._stoichs('active_reac_stoich', keys)

    def all_prod_stoichs(self, keys=None):
        return self._stoichs('all_prod_stoich', keys)

    def active_prod_stoichs(self, keys=None):
        return self._stoichs('active_prod_stoich', keys)

    def stoichs(self, non_precip_rids=()):  # TODO: rename to cond_stoichs
        """ Conditional stoichiometries depending on precipitation status """
        # dtype: see https://github.com/sympy/sympy/issues/10295
        import numpy as np
        return np.array([(
            -np.array(eq.precipitate_stoich(self.substances)[0]) if idx
            in non_precip_rids else
            eq.non_precipitate_stoich(self.substances)
        ) for idx, eq in enumerate(self.rxns)], dtype=object)

    def composition_balance_vectors(self):
        """ Returns a list of lists with compositions and a list of composition keys.

        The list of lists can be viewed as a matrix with rows corresponding to composition keys
        (which are given as the second item in the returned tuple) and columns corresponding to
        substances. Multiplying the matrix with a vector of concentrations give an equation which
        is an invariant (corresponds to mass & charge conservation).

        Examples
        --------
        >>> s = 'Cu+2 + NH3 -> CuNH3+2'
        >>> import re
        >>> substances = re.split(' \+ | -> ', s)
        >>> rsys = ReactionSystem.from_string(s, substances)
        >>> rsys.composition_balance_vectors()
        ([[2, 0, 2], [0, 3, 3], [0, 1, 1], [1, 0, 1]], [0, 1, 7, 29])

        Returns
        -------
        A: list of lists
        ck: (sorted) tuple of composition keys

        """
        subs = self.substances.values()
        ck = Substance.composition_keys(subs)
        return [[s.composition.get(k, 0) for s in subs] for k in ck], ck

    def upper_conc_bounds(self, init_concs, min_=min, dtype=None, skip_keys=(0,)):
        r""" Calculates upper concentration bounds per substance based on substance composition.

        Parameters
        ----------
        init_concs : dict or array_like
            Per substance initial conidtions.
        min_ : callbable
        dtype : dtype or None
        skip_keys : tuple
            What composition keys to skip.

        Returns
        -------
        numpy.ndarray :
            Per substance upper limit (ordered as :attr:`substances`).

        Notes
        -----
        The function does not take into account wheter there actually exists a
        reaction path leading to a substance. Note also that the upper limit is
        per substance, i.e. the sum of all upper bounds amount to more substance than
        available in ``init_conc``.

        Examples
        --------
        >>> rs = ReactionSystem.from_string('2 HNO2 -> H2O + NO + NO2 \n 2 NO2 -> N2O4')
        >>> from collections import defaultdict
        >>> c0 = defaultdict(float, HNO2=20)
        >>> ref = {'HNO2': 20, 'H2O': 10, 'NO': 20, 'NO2': 20, 'N2O4': 10}
        >>> rs.as_per_substance_dict(rs.upper_conc_bounds(c0)) == ref
        True

        """
        import numpy as np
        # A, composition_keys = self.composition_balance_vectors()
        # composition_amounts = np.dot(A, self.as_per_substance_array(init_conc))
        # return np.min(composition_amounts.reshape((composition_amounts.size, 1))/A, axis=0)
        if dtype is None:
            dtype = np.float64
        init_concs_arr = self.as_per_substance_array(init_concs, dtype=dtype)
        composition_conc = defaultdict(float)
        for conc, s_obj in zip(init_concs_arr, self.substances.values()):
            for comp_nr, coeff in s_obj.composition.items():
                if comp_nr in skip_keys:  # charge may be created (if compensated)
                    continue
                composition_conc[comp_nr] += coeff*conc
        bounds = []
        for s_obj in self.substances.values():
            choose_from = []
            for comp_nr, coeff in s_obj.composition.items():
                if comp_nr == 0:
                    continue
                choose_from.append(composition_conc[comp_nr]/coeff)
            if len(choose_from) == 0:
                bounds.append(float('inf'))
            else:
                bounds.append(min_(choose_from))
        return bounds

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
        title : bool

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

    def identify_equilibria(self):
        eq = []
        for ri1, rxn1 in enumerate(self.rxns):
            for ri2, rxn2 in enumerate(self.rxns[ri1+1:], ri1+1):

                all_eq = (rxn1.all_reac_stoich(self.substances) == rxn2.all_prod_stoich(self.substances) and
                          rxn1.all_prod_stoich(self.substances) == rxn2.all_reac_stoich(self.substances))
                if all_eq:
                    eq.append((ri1, ri2))
                    break
        return eq
