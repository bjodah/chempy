#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
from itertools import chain
import os
import re
import shutil
import subprocess
import sys
import warnings

from setuptools import setup

pkg_name = "chempy"
url = 'https://github.com/bjodah/' + pkg_name
license = 'BSD'

RELEASE_VERSION = os.environ.get('%s_RELEASE_VERSION' % pkg_name.upper(), '')  # v*


def _path_under_setup(*args):
    return os.path.join(os.path.dirname(__file__), *args)

release_py_path = _path_under_setup(pkg_name, '_release.py')

if len(RELEASE_VERSION) > 0:
    if RELEASE_VERSION[0] == 'v':
        TAGGED_RELEASE = True
        __version__ = RELEASE_VERSION[1:]
    else:
        raise ValueError("Ill formated version")
else:
    TAGGED_RELEASE = False
    # read __version__ attribute from _release.py:
    exec(open(release_py_path).read())
    if __version__.endswith('git'):
        try:
            _git_version = subprocess.check_output(
                ['git', 'describe', '--dirty']).rstrip().decode('utf-8').replace('-dirty', '.dirty')
        except subprocess.CalledProcessError:
            warnings.warn("A git-archive is being installed - version information incomplete.")
        else:
            if 'develop' not in sys.argv:
                warnings.warn("Using git to derive version: dev-branches may compete.")
                _ver_tmplt = r'\1.post\2' if os.environ.get('CONDA_BUILD', '0') == '1' else r'\1.post\2+\3'
                __version__ = re.sub('v([0-9.]+)-(\d+)-(\S+)', _ver_tmplt, _git_version)  # .dev < '' < .post

submodules = [
    'chempy.electrochemistry',
    'chempy.kinetics',
    'chempy.printing',
    'chempy.properties',
    'chempy.thermodynamics',
    'chempy.util',
]

tests = [
    'chempy.tests',
    'chempy.electrochemistry.tests',
    'chempy.kinetics.tests',
    'chempy.printing.tests',
    'chempy.properties.tests',
    'chempy.thermodynamics.tests',
    'chempy.util.tests',
]

classifiers = [
    'Development Status :: 4 - Beta',
    'License :: OSI Approved :: BSD License',
    'Operating System :: OS Independent',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
]

with io.open(_path_under_setup(pkg_name, '__init__.py'), 'rt', encoding='utf-8') as f:
    short_description = f.read().split('"""')[1].split('\n')[1]
if not 10 < len(short_description) < 255:
    warnings.warn("Short description from __init__.py proably not read correctly")
long_descr = io.open(_path_under_setup('README.rst'), encoding='utf-8').read()
if not len(long_descr) > 100:
    warnings.warn("Long description from README.rst probably not read correctly.")
_author, _author_email = open(_path_under_setup('AUTHORS'), 'rt').readline().split('<')

extras_req = {
    'integrators': ['pyodeint>=0.10.1', 'pycvodes>=0.11.9', 'pygslodeiv2>=0.9.1'],
    'solvers': ['pykinsol>=0.1.3'],
    'native': ['pycompilation>=0.4.3', 'pycodeexport>=0.1.2', 'appdirs'],
    'docs': ['Sphinx', 'sphinx_rtd_theme', 'numpydoc'],
    'plotting': ['bokeh>=0.13.0', 'ipywidgets'],
    'testing': ['pytest', 'pytest-cov', 'pytest-flakes', 'pytest-pep8', 'rstcheck']
}
extras_req['all'] = list(chain(extras_req.values()))

setup_kwargs = dict(
    name=pkg_name,
    version=__version__,
    description=short_description,
    long_description=long_descr,
    author=_author.strip(),
    author_email=_author_email.split('>')[0].strip(),
    url=url,
    license=license,
    keywords=("chemistry", "water properties", "physical chemistry"),
    packages=[pkg_name] + submodules + tests,
    classifiers=classifiers,
    install_requires=[
        'numpy>1.11.3', 'scipy>=1.0.1', 'matplotlib>=2.2.3',
        'sympy>=1.1.1,!=1.2', 'quantities>=0.12.1', 'pyneqsys>=0.5.4',
        'pyodesys>=0.12.4', 'pyparsing>=2.0.3', 'sym>=0.3.4', 'jupyter',
        'pulp>=1.6.8',
        # 'dot2tex>=2.9.0'
    ],
    extras_require=extras_req,
    python_requires='>=3.5'
)

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
