---
title: 'ChemPy: A package useful for chemistry written in Python'
tags:
  - chemistry
  - physical chemistry
authors:
 - name: Björn Dahlgren
   orcid: 0000-0003-0596-0222
   affiliation: 1
affiliations:
 - name: KTH Royal Institute of Technology
   index: 1
date: 26 January 2018
bibliography: paper.bib
---

# Summary
ChemPy is a Python library which provides functions and classes for
solving chemistry related problems. It includes classes for
representing substances, reactions, and systems of reactions. It also
includes well established formulae from physical chemistry, as well as
analytic solutions to some differential equations commonly encountered
in chemical kinetics. Last, but not the least, it collects
parametrizations of chemical properties of substances from the
literature.

The class for substances are represented by name and optionally
contain information on their composition, weight, charge etc., as
well as how to pretty print them in e.g. LaTeX. Both the composition
and pretty printing forms can be deduced by ChemPy's
parser. Reactions are represented through their stoichiometry and
thermodynamic/kinetic parameters. If the stoichiometry of a reaction
is unknown, ChemPy can balance it based on the composition of the
substances. The classes for representing systems of reactions provide
methods to analyze e.g. if there are disjoint sets of reactions, or if all
are connected in the same network. The classes also offer
a series of "checks" to be performed at construction, ensuring balanced
reactions with sane coefficients and consistent units.

Systems of reactions can be represented as graphs, tables, systems of
ordinary differential equations (chemical kinetics) or non-linear
equations (chemical equilibria). The latter two forms can be solved
numerically using pyodesys [@dahlgren_pyodesys_2018] and pyneqsys
[@dahlgren_pyneqsys_2018] respectively.

Thanks to the use of SymPy [@Meurer2017], a user can not only solve
stoichiometry problems with a single unique solution, but also
under-determined systems, where the answer then contains a free parameter.
In fact, most equations and parametrizations in ChemPy support, in addition to
NumPy [@vanderWalt2011] arrays, also symbolic input, as well as arrays
with explicit units. The latter allows ChemPy to check that the correct
dimensionality is used based on e.g. reaction order.

# Features
- Pretty printing of chemical formulae and reaction sets.
- Interactive JavaScript enabled widgets in the Jupyter notebook
  [@Kluyver2016].
- Parsing of chemical formulae, reactions and systems
  thereof.
- Functions for expressing systems of reactions as ordinary
  differential equations.
- Functions for expressing systems of equilibria as non-linear
  equations (including multi-phase systems).
- Analytic solutions for a selection of kinetic problems.


# References
