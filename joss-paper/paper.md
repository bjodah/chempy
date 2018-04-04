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
 - name: KTH Royal Institute of Technology, Stockholm, Sweden
   index: 1
date: 26 January 2018
bibliography: paper.bib
---

# Summary
ChemPy is a Python library that provides functions and classes for
solving chemistry related problems. It includes classes for
representing substances, reactions, and systems of reactions. It also
includes well established formulae from physical chemistry, as well as
analytic solutions to some differential equations commonly encountered
in chemical kinetics. Last, but not the least, it collects
parametrizations of chemical properties of substances from the
literature.

Its intended audience is primarily researchers and engineers who need
to perform modelling work. But since the intermediate representations
of, e.g., ODE systems and systems of non-linear equations are available
symbolically, ChemPy may also be used in an educational setting.

Substances are represented by a class that holds their names and, optionally,
information on their composition, weight, charge etc., as
well as how to pretty print them using LaTeX, HTML and unicode. Both the composition
and stylistic representations can be deduced by ChemPy's
parser. Reactions are represented through their stoichiometry and
thermodynamic/kinetic parameters. If the stoichiometry of a reaction
is unknown, ChemPy can balance it based on the composition of the
substances. The classes for representing systems of reactions provide
methods to analyze, e.g., if there are disjoint sets of reactions, or if all
are connected in the same network. The classes also offer
a series of checks performed at initialization, ensuring balanced
reactions with sane coefficients and consistent units.

Systems of reactions can be represented as graphs, tables, systems of
ordinary differential equations (chemical kinetics) or non-linear
equations (chemical equilibria). The latter two forms can be solved
numerically using pyodesys [@dahlgren_pyodesys_2018] and pyneqsys
[@dahlgren_pyneqsys_2018] respectively.

Thanks to the use of SymPy [@Meurer2017], stoichiometry problems with
a single unique solution can be solved analytically, as well as
under-determined systems, where the answer then contains a free parameter.
The under-determined formulation can also be expressed in a canoncial form
with coefficients minimzed using PuLP [@mitchell2011pulp;@lougee2003common].
In fact, most equations and parametrizations in ChemPy support---in addition to
NumPy arrays [@vanderWalt2011]---also symbolic input, as well as arrays
with explicit units. The latter allows ChemPy to check that, e.g., the correct
dimensionality with respect to reaction order is used for rate constants.

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
