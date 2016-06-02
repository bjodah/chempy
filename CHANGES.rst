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
