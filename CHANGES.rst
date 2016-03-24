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
