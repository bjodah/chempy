#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import shutil
from setuptools import setup

pkg_name = "chempy"

RELEASE_VERSION = os.environ.get('CHEMPY_RELEASE_VERSION', '')

# http://conda.pydata.org/docs/build.html#environment-variables-set-during-the-build-process
if os.environ.get('CONDA_BUILD', '0') == '1':
    try:
        RELEASE_VERSION = 'v' + io.open(
            '__conda_version__.txt', 'rt', encoding='utf-8'
        ).readline().rstrip()
    except IOError:
        pass


def _path_under_setup(*args):
    return os.path.join(os.path.dirname(__file__), *args)

release_py_path = _path_under_setup(pkg_name, '_release.py')

if (len(RELEASE_VERSION) > 1 and RELEASE_VERSION[0] == 'v'):
    TAGGED_RELEASE = True
    __version__ = RELEASE_VERSION[1:]
else:
    TAGGED_RELEASE = False
    # read __version__ attribute from _release.py:
    exec(io.open(release_py_path, encoding='utf-8').read())


submodules = [
    'chempy.kinetics',
    'chempy.properties',
    'chempy.util',
]

tests = [
    'chempy.tests',
    'chempy.kinetics.tests',
    'chempy.properties.tests',
    'chempy.util.tests',
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
    'Programming Language :: Python :: 3.5',
]

with io.open(_path_under_setup(pkg_name, '__init__.py'), 'rt',
             encoding='utf-8') as f:
    short_description = f.read().split('"""')[1].split('\n')[1]
assert 10 < len(short_description) < 255
long_descr = io.open(_path_under_setup('README.rst'), encoding='utf-8').read()
assert len(long_descr) > 100


setup_kwargs = {
    'name': pkg_name,
    'version': __version__,
    'description': short_description,
    'long_description': long_descr,
    'author': 'Björn Dahlgren',
    'author_email': 'bjodah@DELETEMEgmail.com',
    'license': 'BSD',
    'keywords': ("chemistry", "water properties", "physical chemistry"),
    'url': 'https://github.com/bjodah/' + pkg_name,
    'packages': [pkg_name] + submodules + tests,
    'classifiers': classifiers,
    'install_requires': [
        'numpy>1.7', 'scipy>=0.16.1', 'matplotlib>=1.3.1',
        'sympy>=0.7.6.1', 'quantities>=0.11.1', 'pyneqsys>=0.3.0',
        'pyodesys>=0.5.0', 'pyparsing>=2.0.3'
        # 'dot2tex>=2.9.0'
    ],
    'extras_require': {
        'all': ['argh', 'pycvodes', 'pygslodeiv2', 'pyodeint', 'pykinsol',
                'pytest>=2.8.1', 'pytest-pep8>=1.0.6', 'bokeh>=0.11.1']}
}

if __name__ == '__main__':
    try:
        if TAGGED_RELEASE:
            # Same commit should generate different sdist
            # depending on tagged version (set CHEMPY_RELEASE_VERSION)
            # this will ensure source distributions contain the correct version
            shutil.move(release_py_path, release_py_path+'__temp__')
            open(release_py_path, 'wt').write(
                "__version__ = '{}'\n".format(__version__))
        setup(**setup_kwargs)
    finally:
        if TAGGED_RELEASE:
            shutil.move(release_py_path+'__temp__', release_py_path)
