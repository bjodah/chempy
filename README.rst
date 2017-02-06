ChemPy
======

.. image:: http://hera.physchem.kth.se:9090/api/badges/bjodah/chempy/status.svg
   :target: http://hera.physchem.kth.se:9090/bjodah/chempy
   :alt: Build status
.. image:: https://img.shields.io/pypi/v/chempy.svg
   :target: https://pypi.python.org/pypi/chempy
   :alt: PyPI version
.. image:: https://img.shields.io/badge/python-2.7,3.4,3.5-blue.svg
   :target: https://www.python.org/
   :alt: Python version
.. image:: https://zenodo.org/badge/8840/bjodah/chempy.svg
   :target: https://zenodo.org/badge/latestdoi/8840/bjodah/chempy
.. image:: https://img.shields.io/pypi/l/chempy.svg
   :target: https://github.com/bjodah/chempy/blob/master/LICENSE
   :alt: License
.. image:: http://img.shields.io/badge/benchmarked%20by-asv-green.svg?style=flat
   :target: http://hera.physchem.kth.se/~chempy/benchmarks
   :alt: airspeedvelocity
.. image:: http://hera.physchem.kth.se/~chempy/branches/master/htmlcov/coverage.svg
   :target: http://hera.physchem.kth.se/~chempy/branches/master/htmlcov
   :alt: coverage


.. contents::


About ChemPy
------------
`ChemPy <https://github.com/bjodah/chempy>`_ is a `Python <https://www.python.org>`_ package useful for
chemistry (mainly physical/inorganic/analytical chemistry). Currently it includes:

- Solver for equilibria (including multiphase systems)
- Numerical integration routines for chemical kinetics (ODE solver front-end)
- Integrated rate expressions (and convenience fitting routines)
- Relations in physical chemistry:

  - Debye-Hückel expressions
  - Arrhenius equation
  - Einstein-Smoluchowski equation

- Properties

  - water density as function of temperature
  - water permittivity as function of temperature and pressure
  - water diffusivity as function of temperature
  - sulfuric acid density as function of temperature & weight fraction H₂SO₄


Documentation
-------------
Auto-generated API documentation for the latest stable release is found here:
`<https://bjodah.github.io/chempy/latest>`_
(and the development version for the current master branch is found here:
`<http://hera.physchem.kth.se/~chempy/branches/master/html>`_).


Installation
------------
Simplest way to install ChemPy and its (optional) dependencies is to use the `conda package manager <https://conda.pydata.org/docs/>`_::

   $ conda install -c bjodah chempy pytest
   $ python -m pytest --pyargs chempy  # runs the test-suite

alternatively you may also use `pip`::

   $ python -m pip install --user --upgrade chempy pytest
   $ python -m pytest --pyargs chempy

you can skip the ``--user`` flag if you have got root permissions.
See `setup.py <setup.py>`_ for optional requirements.


Examples
--------
See demonstration scripts in `examples/ <https://github.com/bjodah/chempy/tree/master/examples>`_,
and some rendered `jupyter <https://www.jupyter.org>`_ notebooks here:
`<http://hera.physchem.kth.se/~chempy/branches/master/examples>`_.
You may also browse the documentation for more examples. Below you will find a few code snippets:

Parsing formulae
~~~~~~~~~~~~~~~~
.. code:: python

   >>> from chempy import Substance
   >>> ferricyanide = Substance.from_formula('Fe(CN)6-3')
   >>> ferricyanide.composition == {0: -3, 26: 1, 6: 6, 7: 6}  # 0 for charge
   True
   >>> print(ferricyanide.unicode_name)
   Fe(CN)₆³⁻
   >>> print(ferricyanide.latex_name + ", " + ferricyanide.html_name)
   Fe(CN)_{6}^{3-}, Fe(CN)<sub>6</sub><sup>3-</sup>
   >>> print('%.3f' % ferricyanide.mass)
   211.955


as you see, in composition, the atomic numbers (and 0 for charge) is used as
keys and the count of each kind became respective value.

Balancing stoichiometry of a chemical reaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code:: python

   >>> from chempy import balance_stoichiometry  # Main reaction in NASA's booster rockets:
   >>> reac, prod = balance_stoichiometry({'NH4ClO4', 'Al'}, {'Al2O3', 'HCl', 'H2O', 'N2'})
   >>> from pprint import pprint
   >>> pprint(reac)
   {'Al': 10, 'NH4ClO4': 6}
   >>> pprint(prod)
   {'Al2O3': 5, 'H2O': 9, 'HCl': 6, 'N2': 3}
   >>> from chempy import mass_fractions
   >>> for fractions in map(mass_fractions, [reac, prod]):
   ...     pprint({k: '{0:.3g} wt%'.format(v*100) for k, v in fractions.items()})
   ...
   {'Al': '27.7 wt%', 'NH4ClO4': '72.3 wt%'}
   {'Al2O3': '52.3 wt%', 'H2O': '16.6 wt%', 'HCl': '22.4 wt%', 'N2': '8.62 wt%'}


Balancing reactions
~~~~~~~~~~~~~~~~~~~
.. code:: python

   >>> from chempy import Equilibrium
   >>> from sympy import symbols
   >>> K1, K2, Kw = symbols('K1 K2 Kw')
   >>> e1 = Equilibrium({'MnO4-': 1, 'H+': 8, 'e-': 5}, {'Mn+2': 1, 'H2O': 4}, K1)
   >>> e2 = Equilibrium({'O2': 1, 'H2O': 2, 'e-': 4}, {'OH-': 4}, K2)
   >>> coeff = Equilibrium.eliminate([e1, e2], 'e-')
   >>> coeff
   [4, -5]
   >>> redox = e1*coeff[0] + e2*coeff[1]
   >>> print(redox)
   20 OH- + 32 H+ + 4 MnO4- = 26 H2O + 4 Mn+2 + 5 O2; K1**4/K2**5
   >>> autoprot = Equilibrium({'H2O': 1}, {'H+': 1, 'OH-': 1}, Kw)
   >>> n = redox.cancel(autoprot)
   >>> n
   20
   >>> redox2 = redox + n*autoprot
   >>> print(redox2)
   12 H+ + 4 MnO4- = 4 Mn+2 + 5 O2 + 6 H2O; K1**4*Kw**20/K2**5


Chemical equilibria
~~~~~~~~~~~~~~~~~~~
.. code:: python

   >>> from chempy import Equilibrium
   >>> from chempy.chemistry import Species
   >>> water_autop = Equilibrium({'H2O'}, {'H+', 'OH-'}, 10**-14)  # unit "molar" assumed
   >>> ammonia_prot = Equilibrium({'NH4+'}, {'NH3', 'H+'}, 10**-9.24)  # same here
   >>> from chempy.equilibria import EqSystem
   >>> substances = map(Species.from_formula, 'H2O OH- H+ NH3 NH4+'.split())
   >>> eqsys = EqSystem([water_autop, ammonia_prot], substances)
   >>> print('\n'.join(map(str, eqsys.rxns)))  # "rxns" short for "reactions"
   H2O = H+ + OH-; 1e-14
   NH4+ = H+ + NH3; 5.75e-10
   >>> from collections import defaultdict
   >>> init_conc = defaultdict(float, {'H2O': 1, 'NH3': 0.1})
   >>> x, sol, sane = eqsys.root(init_conc)
   >>> assert sol['success'] and sane
   >>> print(sorted(sol.keys()))  # see package "pyneqsys" for more info
   ['fun', 'intermediate_info', 'internal_x_vecs', 'nfev', 'njev', 'success', 'x', 'x_vecs']
   >>> print(', '.join('%.2g' % v for v in x))
   1, 0.0013, 7.6e-12, 0.099, 0.0013

Please note that the API of the ``chempy.equilibria`` module is not finalized at
the moment.

Ionic strength
~~~~~~~~~~~~~~
.. code:: python

   >>> from chempy.electrolytes import ionic_strength
   >>> ionic_strength({'Fe+3': 0.050, 'ClO4-': 0.150}) == .3
   True

Chemical kinetics (system of ordinary differential equations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code:: python

   >>> from chempy import ReactionSystem  # The rate constants below are arbitrary
   >>> rsys = ReactionSystem.from_string("""2 Fe+2 + H2O2 -> 2 Fe+3 + 2 OH-; 42
   ...     2 Fe+3 + H2O2 -> 2 Fe+2 + O2 + 2 H+; 17
   ...     H+ + OH- -> H2O; 1e10
   ...     H2O -> H+ + OH-; 1e-4
   ...     Fe+3 + 2 H2O -> FeOOH(s) + 3 H+; 1
   ...     FeOOH(s) + 3 H+ -> Fe+3 + 2 H2O; 2.5""")  # "[H2O]" = 1.0 (actually 55.4 at RT)
   >>> from chempy.kinetics.ode import get_odesys
   >>> odesys, extra = get_odesys(rsys)
   >>> from collections import defaultdict
   >>> import numpy as np
   >>> tout = sorted(np.concatenate((np.linspace(0, 23), np.logspace(-8, 1))))
   >>> c0 = defaultdict(float, {'Fe+2': 0.05, 'H2O2': 0.1, 'H2O': 1.0, 'H+': 1e-7, 'OH-': 1e-7})
   >>> result = odesys.integrate(tout, c0, atol=1e-12, rtol=1e-14)
   >>> import matplotlib.pyplot as plt
   >>> _ = plt.subplot(1, 2, 1)
   >>> _ = result.plot(names=[k for k in rsys.substances if k != 'H2O'])
   >>> _ = plt.legend(loc='best', prop={'size': 9}); _ = plt.xlabel('Time'); _ = plt.ylabel('Concentration')
   >>> _ = plt.subplot(1, 2, 2)
   >>> _ = result.plot(names=[k for k in rsys.substances if k != 'H2O'], xscale='log', yscale='log')
   >>> _ = plt.legend(loc='best', prop={'size': 9}); _ = plt.xlabel('Time'); _ = plt.ylabel('Concentration')
   >>> _ = plt.tight_layout()
   >>> plt.show()  # doctest: +SKIP

.. image:: https://raw.githubusercontent.com/bjodah/chempy/master/examples/kinetics.png


License
-------
The source code is Open Source and is released under the very permissive
`"simplified (2-clause) BSD license" <https://opensource.org/licenses/BSD-2-Clause>`_.
See `LICENSE <LICENSE>`_ for further details.

See also
--------
- `SymPy <https://github.com/sympy/sympy>`_
- `pyneqsys <https://github.com/bjodah/pyneqsys>`_
- `pyodesys <https://github.com/bjodah/pyodesys>`_
- `thermo <https://github.com/CalebBell/thermo>`_

Contributing
------------
Contributors are welcome to suggest improvements at https://github.com/bjodah/chempy


Author
------
Björn I. Dahlgren, contact:
 - gmail address: bjodah
 - kth.se address: bda
