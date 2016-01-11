v0.2.0
======
- Signature of chempy.equilibria.roots, changed.
- Added two new modules: chempy.util.table, chempy.util.graph
- Reaction, ReactionSystems now expects stoichs etc. to be given wrt to Substance names.
- Added chempy.chemistry.ArrheniusRate
- EqSystemLog, EqSystemLin -> EqSystem, (NumSysLog, NumSysLin)
- Support for solid phases in equilibria
- Submodules for water properties moved to chempy.properties
- Moved class ``Equilibrium`` from .equilibria to .chemistry
- Renamed Reaction.params to Reaction.param

v0.1
====
- Initial release
