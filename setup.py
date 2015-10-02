#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

pkg_name = "aqchem"
exec(open(pkg_name + '/_release.py').read())

with open(pkg_name + '/__init__.py') as f:
    long_description = f.read().split('"""')[1]

submodules = [
    'aqchem.kinetics'
]

tests = [
    'aqchem.tests',
    'aqchem.kinetics.tests',
]

classifiers = [
    "Development Status :: 3 - Alpha",
    'License :: OSI Approved :: BSD License',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
]

setup_kwargs = {
    'name': pkg_name,
    'version': eval('__version__'),  # silence pyflakes
    'description': (
        'Package useful for (physical) chemistry'
    ),
    'long_description': long_description,
    'author': 'Björn Dahlgren',
    'author_email': 'bjodah@DELETEMEgmail.com',
    'license': 'BSD',
    'keywords': ("chemistry", "water properties"),
    'url': 'https://github.com/bjodah/' + pkg_name,
    'packages': [pkg_name] + submodules + tests,
    'classifiers': classifiers,
}

if __name__ == '__main__':
    setup(**setup_kwargs)
