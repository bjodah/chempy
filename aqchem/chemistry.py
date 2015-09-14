# -*- coding: utf-8 -*-

from __future__ import division

from operator import itemgetter
from collections import defaultdict

from aqchem.equilibria import equilibrium_quotient


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
    elemental_composition: dict {atomic number: count}
    """
    def __init__(self, name=None, charge=None, mass=None,
                 latex_name=None, formula=None,
                 elemental_composition=None):
        self.name = name
        self.charge = charge
        self.mass = mass
        self.latex_name = latex_name
        self.formula = formula
        self.elemental_composition = elemental_composition

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
        if elemental_composition is None:
            # try:
            self.elemental_composition = {
                k.number: v for k, v in elements(formula).items()}
            # except AttributeError:
            #     pass


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

    def active_net_stoich(self, substances):
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


class Equilibrium(Reaction):
    """
    Represents equilibrium reaction, subclass of Reatcion.
    Extra attributes:
    - solubility_product: bool (default: False)
    """


    latex_arrow = '\\rightleftharpoons'

    # def __init__(self, *args, **kwargs):
    #     self.solubility_product = kwargs.pop('solubility_product', False)
    #     super(self.__class__, self).__init__(*args, **kwargs)

    def K(self, T=None):
        if callable(self.params):
            if T is None:
                raise ValueError("No T provided")
            return self.params(T)
        else:
            if T is not None:
                raise ValueError("T provided but params not callble")
            return self.params

    def solid_factor(self, substances, concs):
        factor = 1
        for r, n in self.reac.items():
            if r.solid:
                factor *= concs[substances.index(r)]**-n
        for p, n in self.prod.items():
            if p.solid:
                factor *= concs[substances.index(p)]**n
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

class ReactionSystem(object):

    def __init__(self, rxns, substances):
        self.rxns = rxns
        self.substances = substances

def atom_balance(substances, concs, atom_number):
    res = 0
    for s, c in zip(substances, concs):
        res += s.elemental_composition.get(atom_number, 0)*c
    return res

def charge_balance(substances, concs):
    res = 0
    for s, c in zip(substances, concs):
        res += s.charge*c
    return res


def dot(iter_a, iter_b):
    return sum([a*b for a, b in zip(iter_a, iter_b)])

class EqSystem(ReactionSystem):

    def equilibrium_quotients(self, c):
        return [equilibrium_quotient(
            c, eq.active_net_stoich(
                self.substances)) for eq in self.rxns]

    def charge_balance(self, concs):
        return charge_balance(self.substances, concs)

    def charge_balance_vector(self):
        return [s.charge for s in self.substances]

    def atom_balance(self, concs, atom_number):
        # Convenience method
        return atom_balance(self.substances, concs, atom_number)

    def atom_balance_vectors(self):
        atom_numbers = set()
        for s in self.substances:
            for key in s.elemental_composition:
                atom_numbers.add(key)
        vectors = []
        sorted_atom_numbers = sorted(atom_numbers)
        for atm_nr in sorted_atom_numbers:
            vectors.append([s.elemental_composition.get(atm_nr, 0) for
                            s in self.substances])
        return vectors, sorted_atom_numbers

    def rref(self):
        """
        Calculates the Reduced Row Echelon Form of the linear equations
        system for charge and atom balance.

        Returns
        -------
        Length-2 tuple:
           - list of lists with coefficients for the substances
           - pivot array
        """
        from sympy import Matrix
        atom_vecs, atm_nrs = self.atom_balance_vectors()
        M = Matrix([self.charge_balance_vector()] + atom_vecs)
        rref, pivot = M.rref()
        #indep_atm_nrs =
        return rref.tolist()[:len(pivot)], pivot #indep_atm_nrs

    def independent_atoms_from_pivot(self, pivot):
        """ atom number 0 represents charge """
        atom_vecs, atm_nrs = self.atom_balance_vectors()
        return [([0] + atm_nrs)[idx] for idx in pivot]

    def qk(self, concs, scaling=1.0):
        return [q - scaling**rxn.dimensionality()*rxn.params*rxn.solid_factor(
            self.substances, concs) for q, rxn in zip(
                self.equilibrium_quotients(concs), self.rxns)]

    def preserved(self, concs, init_concs):
        rref, pivot = self.rref()
        p0s = [dot(row, init_concs) for row in rref]
        return [dot(row, concs) - p0 for row, p0 in zip(rref, p0s)], pivot

    def f(self, concs, init_concs, scaling=1, reduced=None):
        pres, pivot = self.preserved(concs, init_concs)
        if reduced is None:
            return self.qk(concs, scaling) + pres
        else:
            if reduced:
                subs = []
                for idx, p in zip(pivot, pres)[::-1]:
                    c = concs[idx]
                    subs.append((c, (c-p.subs(subs)).simplify()))
                qk = [expr.subs(subs) for expr in self.qk(concs, scaling)]
                return qk, pivot
            else:
                return self.qk(concs, scaling) + pres, None

    def num_cb_factory(self, init_concs, jac=False, scaling=1.0, logC=False,
                       square=False, reduced=None):
        import sympy as sp
        y = sp.symarray('y', len(self.substances))
        f, elim = self.f(y, init_concs, scaling, reduced=bool(reduced))
        if reduced is not None:
            if reduced:
                y = [y[idx] for idx in range(len(y)) if idx not in elim]
                print(y)
        if square:
            f = [_.subs([(yi, yi*yi) for yi in y]) for _ in f]
        if logC:
            f = [_.subs([(yi, sp.exp(yi)) for yi in y]) for _ in f]
        f_lmbd = sp.lambdify(y, f, modules='numpy')
        def f_cb(arr):
            print(arr) ## DEBUG
            return f_lmbd(*arr)
        j_cb = None
        if jac:
            print(len(f), len(y))
            j = sp.Matrix(1, len(y), lambda _, q: f[q]).jacobian(y)
            j_lmbd = sp.lambdify(y, j, modules='numpy')
            def j_cb(arr):
                return j_lmbd(*arr)
        if reduced is None:
            return f_cb, j_cb
        else:
            return f_cb, j_cb, elim

    def root(self, init_concs, scaling=1.0, logC=False, square=False,
             delta=None, reduced=False, **kwargs):
        import numpy as np
        from scipy.optimize import root
        init_concs = np.asarray([init_concs[k] for k in self.substances] if
                             isinstance(init_concs, dict) else init_concs)
        f, j, elim = self.num_cb_factory(init_concs*scaling, jac=True,
                                         scaling=scaling, logC=logC,
                                         square=square, reduced=reduced)
        if delta is None:
            delta = kwargs.get('tol', 1e-12)
        x0 = init_concs*scaling + delta
        if reduced:
            x0 = np.array([x for idx, x in enumerate(x0) if idx not in elim])
        if square:
            x0 = x0**0.5
        if logC:
            x0 = np.log(x0)
        result = root(f, x0, jac=j)
        if result.success:
            # TODO: back-transform result
            x = result.x
            if logC:
                x = np.exp(x)
            if square:
                x *= x
            x /= scaling
            print('Success: ', x)
        if reduced:
            return result # TODO: reintroduce eliminated...
        else:
            return result
