v0.7.10
=======
- Added support for newer SymPy version.

v0.7.9
======
- Fixed regression in balance_stoichiometry
- More robust expression handling.
- Fixed SyntaxWarning for upcoming CPython 3.8

v0.7.8
======
- Reaction __init__ got a new kwarg: ``dont_check``
- Fixes to check logic
- Fixes to pretty printing
- ArrheniusParam got a new convenience classmethod: ``from_rateconst_at_T``
- ReactionSystem got a new method: ``concatenate``
- New submodule for template rendering: .util.rendering

v0.7.7
======
- Fixed a bug in ``chempy.kinetics.ode._validate``
- ``balance_stoichiometry`` now is even more allowing with "allow_duplicates"

v0.7.6
======
- ``balance_stoichiometry`` got a new keyword argument "allow_duplicates" (see gh-120)

v0.7.5
======
- Made chempy.units.to_unitless more strict

v0.7.4
======
- Fixed bug with quantities.constants when used with Arrhenius class

v0.7.3
======
- Fixed an issue with PuLP under FreeBSD (gh-109)

v0.7.2
======
- Fixed bug in Reaction.precipitate_stoich
- Changed Python-3-only syntax to support Python 2
- Python < 3.5 is now deprecated for ChemPy <0.7

v0.7.1
======
- Updated requirements on upstream packages

v0.7.0
======
- Drop official support for Python 2.7

v0.6.10
=======
- Fixed bug in Reaction.precipitate_stoich
- Changed Python-3-only syntax to support Python 2

v0.6.8
======
- Fix use of numpy.linalg.lstsq
- Upper limit on pyodesys (<0.12), chempy 0.7.x will use pyodesys>=0.12

v0.6.7
======
- Updated manuscript for JOSS.

v0.6.6
======
- Support for containers in ``chempy.units.unit_of``.
- ``balance_stoichiometry`` now correctly find the canoncial solution when ``underdetermined=None`` is passed.
- ``chempy.kinetics.integrated`` was refactored to have more approachable API.

v0.6.5
======
- Fix assertion firing

v0.6.4
======
- Enhancements for ``.units.is_unitless`` & ``.units.get_physical_quantity``
- New functions in ``chempy.units``: ``compare_equality`` & ``uniform``.

v0.6.3
======
- Fix bug in Reaction.check_consistent_units

v0.6.2
======
- More relaxed tests with respect to 3rd party programs

v0.6.1
======
- Extensive test suite in conda pacakge no longer require graphviz & latex

v0.6.0
======
- ``balance_stoichiometry`` now accepts either of ``True``, ``False``, ``None`` as ``underdetermined``.
- ``NameSpace`` and ``AttributeContainer`` are now public in ``.util.pyutil``.
- New printers in ``chempy.printing``, allows user to subclass printers.
- Jupyter Notebook representation of ``ReactionSystem`` is now interactive (JavaScript/CSS)
- More data from the litterature: water viscosity (``chempy.properties.water_viscosity_korson_1969``).
- New methods for ``ReactionSystem``:
  - ``split``: splits reaction-system into disjoint parts
  - ``categorize_substances``: e.g. "nonparticipating", "unaffected".
- Better documentation throughout.

v0.5.7
======
- New option in ``.kinetics._native.get_native``: conc_roots

v0.5.6
======
- New method: ``ReactionSystem.sort_substances_inplace()``
- New patched NumPy "module": ``.units.patched_numpy``
- Updated ``.util.bkh.integration_with_sliders`` to be compatible with
  latest bokeh.

v0.5.5
======
- Fix non-deterministic ordering of dictionary in ``get_odesys()``.

v0.5.4
======
- Fix to bokeh interface (``chempy.util.bkh``).

v0.5.3
======
- Fixes balance_stoichiometry
- Documentation fixes
- More k_fmt & landscape options in .util.table.rsys2pdf_table

v0.5.2
======
- Fix balance_reacions (non-deterministic ordering could cause endless loop)
- Fix unit scaling of .kinetics.rates.Eyring

v0.5.1
======
- Moved ReactionSystem to .reactionsystem, (import directly from chempy).
- Steady state analysis
- now in default_units: molar, milli-, micro- & nano-
- CSTR kinetics
- Minor fixes, new notebooks

v0.5.0
======
- ``.electrochemistry.nernst_formula`` - thanks to Adel Qalieh (@adelq)
- moved ``.util.parsing.number_to_scientific_*`` to ``.printing(.numbers)``
- Number formating now handles uncertainties.
- ``refereence`` in reimplementations now a dict
- Fixes to ``.kinetics.ode.get_odesys`` (refactored)

v0.4.1
======
- Fixes for enhanced robustness:
  - ``.kinetics.ode.get_odesys``
  - ``.chemistry.as_per_substance_array``
- Minor changes.

v0.4.0
======
- Multiple fixes throughout
- Refactored .equilibria
- .core and .debye_huckel was merged into .electrolytes
- New functions: balance_stoichiometry, mass_fractions
- kwargs one=, exp=, ln= changed throughout to use backend=None (backen=math)
- .chemistry.ArrheniusRate moved (and changed) to .arrhenius.ArrheniusParam
- Equilibrium got a new method: cancel and a new staticmethod: eliminate
- Reaction now raises ValueError if the Reaction has a zero net effect.
- It is now possible to use (parts of) chempy even when only Python stdlib is available
- Substance got a new method: molar_mass, and a two new attributes: unicode_name, html_name
- .util.parsing.to_latex was renamed to formula_to_latex.
- New functions in util.parsing: formula_to_unicode, formula_to_html
- Parsing of crystal water now supported.
- ReactionSystem.__init__ got a new kwarg: substance_factory
- ReactionSystem raises ValueError if it contains duplicate instances of Reaction
- ReactionSystem got new methods:
  - as_per_substance_dict (inverse of as_per_substance_array)
  - unimolecular_html_table
  - bimolecular_html_table
- .kinetics.ode.law_of_mass_action_rates was updated to handle RateExpr
- fix in .properties.sulfuric_acid_density_myhre_1998.density_from_concentration for input with units
- enhancements to .util.deprecation.Deprecation
- .util.stoich.decompose_yields now takes iterable of Reaction instances as second arg.
- .util.table.rsys2tablines now pretty-prints ref={'doi': 'abc123'} too.
- ``chempy.util.stoich.decompose_yields`` now takes reactions instead of
  iterable of dicts (backward incompatible change).

v0.3.5
======
- More robust setup.py

v0.3.3
======
- ``chempy.units.allclose`` now handles iterables with disparate units.

v0.3.2
======
- Substance.from_formula now prefers e.g. Fe+3 over Fe/3+, latter deprecated

v0.3.1
======
- chemistry.Solute deprecated, will be removed in v0.4.0, use chemistry.Species instead
- ReactionSystem now handles "substances" argument more robustely.

v0.3.0
======
- Signature of chempy.chemistry.Substance changed
- New module chempy.util.parsing, (drop dependency on periodictable)
- EqSystem.root and EqSystem.roots got new kwarg: neqsys_type
- chemistry.Equilibrium learned to handle inactive reactants/products
- chemistry.Reaction dropped kwarg 'k' (deprecated since v0.2.0)

v0.2.0
======
- Signature of chempy.equilibria.roots, changed.
- Added two new modules: chempy.util.table, chempy.util.graph
- chempy.einstein_smoluchowski added
- Reaction, ReactionSystems now expects stoichs etc. to be given wrt to Substance names.
- Added chempy.chemistry.ArrheniusRate
- EqSystemLog, EqSystemLin -> EqSystem, (NumSysLog, NumSysLin)
- Support for solid phases in equilibria
- Submodules for water properties moved to chempy.properties
- Moved class ``Equilibrium`` from .equilibria to .chemistry
- Renamed Reaction.params to Reaction.param
- Added method: Reaction.order()
- Added chempy.properties.sulfuric_acid_density_myhre_1998

v0.1.0
======
- Initial release
