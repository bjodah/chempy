======
chempy
======

.. image:: http://hera.physchem.kth.se:8080/github.com/bjodah/chempy/status.svg?branch=master
   :target: http://hera.physchem.kth.se:8080/github.com/bjodah/chempy
   :alt: Build status on hera

Python package useful for (physical) chemistry. Currently includes:

- water density
- water permittivity
- Debye-H\"uckel expressions
- Solver for equilibria
- Integrated rate expressions (and convenience fitting routines)

Optional dependencies
=====================

- numpy
- scipy
- quantities
- matplotlib
- sympy

Units
=====
Use of quantities assumed, but only following attributes are accessed (parenthesis optional):

- Kelvin
- bar
- meter
- kilogram
- mol
- (nanometer)

Tests
=====
Run ``py.test``

License
=======
The source code is Open Source and is released under the very permissive
"simplified (2-clause) BSD license". See ``LICENSE.txt`` for further details.
Contributors are welcome to suggest improvements at https://github.com/bjodah/chempy

Author
======
Björn Dahlgren, contact:
 - gmail adress: bjodah
 - kth.se adress: bda
